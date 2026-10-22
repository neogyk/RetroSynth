"""
build_dataset.py - build a multi-route polymer retrosynthesis dataset for LLM fine-tuning.

Inputs
  * seed_commercial.py      curated commercial polymers, several routes each (incl. enzymatic)
  * polyRETRO-format CSVs   PID, Polymerization_type, Template_1, Polymer_smiles,
                            reactant1, reactant2, Reaction_smiles, functional_group_column ...

Outputs (in --out)
  routes_kb.jsonl     one record per polymer (frame-invariant key) with ALL its routes
  qc_report.csv       every route with QC status (parse / mass balance / by-product)
  sft_{train,val,test}.jsonl   chat-format examples ({"messages": [...], "task": ..., "polymer_key": ...})
                      directly usable by TRL SFTTrainer / Unsloth / Axolotl ("messages" column)
  split_manifest.json which polymer families/keys went to which split

Usage
  python build_dataset.py --csv train_condensation.csv train_addition.csv train_2k.csv \
      --out data_out --n-aug 4 --test-families polyamide fluoropolymer paek --seed 13
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from polychem import (canonical, canonical_ru_key, formula_str, frame_shifts, infer_byproduct,
                      is_valid, mass_balance, mol_formula, randomized_smiles, star_count)
from seed_commercial import SEED

SYSTEM = ("You are a polymer chemist specialising in retrosynthesis. Polymers are given as "
          "repeat units in PSMILES, with [*] marking the two connection points. Answer with "
          "valid JSON only.")

ROUTE_FIELDS = ["rank", "mechanism_class", "route_type", "process", "monomers", "coreactants",
                "byproducts", "catalyst", "conditions", "scale"]

# --------------------------------------------------------------------------- #
# 1. curated seed
# --------------------------------------------------------------------------- #
def load_seed() -> tuple[dict, list]:
    polymers, qc = {}, []
    for p in SEED:
        ru = canonical(p["ru"])
        if star_count(ru) != 2:
            raise ValueError(f"{p['abbr']}: repeat unit must have two [*]")
        key = canonical_ru_key(ru)
        rec = dict(key=key, ru=ru, name=p["name"], abbr=p["abbr"], family=p["family"],
                   aliases=p["aliases"], commercial=True, routes=[])
        for r in p["routes"]:
            reactants = {canonical(s): v for s, v in {**r["monomers"], **r["coreactants"]}.items()}
            byp = {canonical(s): v for s, v in r["byproducts"].items()}
            ok, resid = mass_balance(ru, reactants, byp, r["ru_per_set"])
            pre_ok = all(_balanced_small(st["reactants"], st["products"])
                         for st in r["precursor_steps"])
            route = dict(r, source="curated",
                         monomers={canonical(s): v for s, v in r["monomers"].items()},
                         coreactants={canonical(s): v for s, v in r["coreactants"].items()},
                         byproducts=byp)
            route["qc"] = "ok" if ok and pre_ok else f"UNBALANCED residual={formula_str(resid)}"
            qc.append(dict(source="curated", polymer=p["abbr"], route_id=r["route_id"], ru=ru,
                           monomers=".".join(route["monomers"]), byproducts=json.dumps(byp),
                           scale=r["scale"], qc=route["qc"]))
            rec["routes"].append(route)
        polymers[key] = rec
    return polymers, qc


def _balanced_small(reactants: dict, products: dict) -> bool:
    lhs, rhs = Counter(), Counter()
    for s, v in reactants.items():
        for k, n in mol_formula(s).items():
            lhs[k] += v * n
    for s, v in products.items():
        for k, n in mol_formula(s).items():
            rhs[k] += v * n
    return all(abs(lhs[k] - rhs[k]) < 1e-6 for k in set(lhs) | set(rhs))


# --------------------------------------------------------------------------- #
# 2. polyRETRO-format CSVs
# --------------------------------------------------------------------------- #
CLASS_MAP = {"condensation": "Condensation", "addition": "Addition", "ringopening": "Ringopening",
             "ring-opening": "Ringopening", "ring opening": "Ringopening", "romp": "ROMP"}


def _first(row, *cols):
    for c in cols:
        v = row.get(c)
        if isinstance(v, str) and v.strip() and v.strip().lower() != "nan":
            return v.strip()
    return None


def load_polyretro(paths: list[str]) -> tuple[list, list]:
    routes, qc = [], []
    for path in paths:
        df = pd.read_csv(path, dtype=str, keep_default_na=False)
        for i, row in df.iterrows():
            row = row.to_dict()
            psmi = _first(row, "Polymer_smiles")
            ptype = (_first(row, "Polymerization_type") or "").lower()
            r1 = _first(row, "reactant1", "reactant_1", "Monomer_smiles")
            r2 = _first(row, "reactant2", "reactant_2")
            pid = _first(row, "PID") or f"{Path(path).stem}_row{i}"
            base = dict(source="polyretro", file=Path(path).name, polymer=pid, route_id=f"{pid}#{i}",
                        ru=psmi, monomers=".".join(x for x in (r1, r2) if x))
            if not psmi or not r1 or not is_valid(psmi) or not is_valid(r1) or (r2 and not is_valid(r2)):
                qc.append(dict(base, qc="DROP: missing/invalid SMILES"))
                continue
            if star_count(psmi) != 2:
                qc.append(dict(base, qc="DROP: repeat unit needs exactly two [*]"))
                continue
            ru = canonical(psmi)
            m1 = canonical(r1)
            m2 = canonical(r2) if r2 else None
            # polyRETRO repeats the monomer in reactant2 for A-A/addition/AB cases -> one monomer type
            monomers = {m1: 1.0} if (m2 is None or m2 == m1) else {m1: 1.0, m2: 1.0}
            pclass = CLASS_MAP.get(ptype, ptype.title() or "Unknown")
            n_ru = 1.0
            byp, resid = infer_byproduct(ru, monomers)
            if byp is None and pclass in ("Ringopening", "ROMP", "Addition") and len(monomers) == 1:
                # cyclic oligomer monomers (e.g. cyclophosphazene trimer): try n_ru = 2..6
                for k in range(2, 7):
                    ok, _ = mass_balance(ru, monomers, {}, k)
                    if ok:
                        byp, n_ru, resid = {}, float(k), Counter()
                        break
            status = "ok" if byp is not None else f"UNBALANCED residual={formula_str(resid)}"
            route = dict(route_id=base["route_id"], source="polyretro", pid=pid,
                         polyretro_class=pclass, mechanism_class=None, route_type=None, process=None,
                         monomers=monomers, coreactants={}, byproducts=byp or {}, ru_per_set=n_ru,
                         catalyst="", conditions="", scale="literature", rank=None,
                         template=_first(row, "Template_1"),
                         functional_groups=_first(row, "functional_group_column",
                                                  "original_functional_group_column"),
                         qc=status, ru=ru, key=canonical_ru_key(ru))
            routes.append(route)
            qc.append(dict(base, ru=ru, monomers=".".join(monomers), byproducts=json.dumps(byp),
                           scale="literature", qc=status))
    return routes, qc


# --------------------------------------------------------------------------- #
# 3. merge
# --------------------------------------------------------------------------- #
def merge(seed: dict, pr_routes: list) -> dict:
    kb = {k: dict(v) for k, v in seed.items()}
    for r in pr_routes:
        key = r.pop("key")
        ru = r.pop("ru")
        if key not in kb:
            kb[key] = dict(key=key, ru=ru, name=None, abbr=None, family=None, aliases=[],
                           commercial=False, routes=[])
        rec = kb[key]
        mset = frozenset(r["monomers"])
        if any(frozenset(x["monomers"]) == mset for x in rec["routes"]):
            continue  # same disconnection already present (curated wins)
        rec["routes"].append(r)
    # polymers without a curated family get a pseudo-family = own key (grouped split)
    for rec in kb.values():
        if not rec["family"]:
            rec["family"] = "pr_" + hashlib.md5(rec["key"].encode()).hexdigest()[:8]
    return kb


# --------------------------------------------------------------------------- #
# 4. split (by polymer family -> no leakage of alternative routes across splits)
# --------------------------------------------------------------------------- #
def split(kb: dict, test_families: list[str], val_frac: float, test_frac: float, seed: int):
    fams = sorted({r["family"] for r in kb.values()})
    rng = random.Random(seed)
    assign = {}
    for f in fams:
        if f in test_families:
            assign[f] = "test"
            continue
        if f.startswith("pr_"):
            x = rng.random()
            assign[f] = "test" if x < test_frac else "val" if x < test_frac + val_frac else "train"
        else:
            assign[f] = "train"
    # make sure val has at least one curated family if possible
    curated_train = [f for f in fams if not f.startswith("pr_") and assign[f] == "train"]
    if curated_train and not any(assign[f] == "val" and not f.startswith("pr_") for f in fams):
        assign[rng.choice(curated_train)] = "val"
    return assign


# --------------------------------------------------------------------------- #
# 5. SFT example generation
# --------------------------------------------------------------------------- #
def _route_json(r: dict) -> dict:
    d = {"monomers": list(r["monomers"].keys())}
    stoich = {s: v for s, v in r["monomers"].items() if v != 1}
    if stoich:
        d["monomer_stoichiometry"] = stoich
    if r.get("ru_per_set", 1) != 1:
        d["repeat_units_per_monomer_set"] = r["ru_per_set"]
    for f in ROUTE_FIELDS:
        if f in ("monomers",):
            continue
        v = r.get(f)
        if v in (None, "", {}, []):
            continue
        d[f] = v
    if r.get("polyretro_class"):
        d["polymerization_class"] = r["polyretro_class"]
    if r.get("precursor_steps"):
        d["precursor_steps"] = r["precursor_steps"]
    if r.get("notes"):
        d["notes"] = r["notes"]
    return d


def _good(r):
    return r["qc"] == "ok" and r.get("scale") != "infeasible"


def _bad(r):
    return r["qc"] == "ok" and r.get("scale") == "infeasible"


def _msg(user, assistant):
    return [{"role": "system", "content": SYSTEM},
            {"role": "user", "content": user},
            {"role": "assistant", "content": json.dumps(assistant, ensure_ascii=False)}]


def target_variants(ru: str, n_aug: int, seed: int) -> list[str]:
    if n_aug <= 1:
        return [ru]
    pool = set(frame_shifts(ru))
    for f in list(pool):
        pool.update(randomized_smiles(f, n=2, seed=seed))
    pool.discard(ru)
    rng = random.Random(seed)
    extra = rng.sample(sorted(pool), min(n_aug - 1, len(pool)))
    return [ru] + extra


def examples_for(rec: dict, kb_list: list, n_aug: int, rng: random.Random) -> list[dict]:
    good = sorted([r for r in rec["routes"] if _good(r)],
                  key=lambda r: (r.get("rank") or 50, r["route_id"]))
    bad = [r for r in rec["routes"] if _bad(r)]
    if not good:
        return []
    out = []
    ident = {k: rec[k] for k in ("name", "abbr") if rec.get(k)}
    for t_i, tgt in enumerate(target_variants(rec["ru"], n_aug, rng.randint(0, 10**6))):
        # (a) enumerate all routes, ranked
        ans = {**ident, "n_routes": len(good), "routes": [_route_json(r) for r in good]}
        if bad:
            ans["rejected_disconnections"] = [
                {"monomers": list(r["monomers"]), "route_type": r.get("route_type"),
                 "reason": r.get("notes")} for r in bad]
        out.append(dict(task="all_routes", messages=_msg(
            f"Target polymer repeat unit: {tgt}\nList every established polymerization route to "
            f"this polymer, ranked from most to least industrially relevant.", ans)))

        if t_i > 0:   # the remaining tasks only on the canonical form + 1 variant to limit size
            if t_i > 1:
                continue
        # (b) polymerization class (polyRETRO objective 1)
        classes = sorted({r["polyretro_class"] for r in good if r.get("polyretro_class")})
        out.append(dict(task="class", messages=_msg(
            f"Target polymer repeat unit: {tgt}\nWhich polymerization classes (Condensation, "
            f"Addition, Ringopening, ROMP) can produce it?", {"classes": classes})))

        # (c) route-type conditioned retrosynthesis (e.g. 'enzymatic')
        by_type = defaultdict(list)
        for r in good:
            if r.get("route_type"):
                by_type[r["route_type"]].append(r)
        for rtype, rs in by_type.items():
            out.append(dict(task="route_conditioned", messages=_msg(
                f"Target polymer repeat unit: {tgt}\nPropose a route of type '{rtype}': monomers, "
                f"catalyst/initiator and conditions.", {"routes": [_route_json(r) for r in rs]})))

        # (d) feasibility: positives, curated negatives, and mass-balance hard negatives
        r = rng.choice(good)
        out.append(dict(task="feasibility", messages=_msg(
            f"Target polymer repeat unit: {tgt}\nProposed monomers: {'.'.join(r['monomers'])}"
            f"{' via ' + r['route_type'] if r.get('route_type') else ''}\nIs this a valid route?",
            {"feasible": True, "scale": r.get("scale")})))
        for b in bad:
            out.append(dict(task="feasibility", messages=_msg(
                f"Target polymer repeat unit: {tgt}\nProposed monomers: {'.'.join(b['monomers'])}"
                f" via {b.get('route_type')}\nIs this a valid route?",
                {"feasible": False, "reason": b.get("notes")})))
        decoy = _decoy(rec, kb_list, rng)
        if decoy:
            out.append(dict(task="feasibility", messages=_msg(
                f"Target polymer repeat unit: {tgt}\nProposed monomers: {'.'.join(decoy)}\n"
                f"Is this a valid route?",
                {"feasible": False, "reason": "Monomers cannot give this repeat unit: the atom "
                 "balance (monomers minus repeat unit) is not a plausible small-molecule by-product."})))

    # (e) forward check: monomers + route -> repeat unit (once per route, canonical only)
    for r in good:
        out.append(dict(task="forward", messages=_msg(
            f"Monomers: {'.'.join(r['monomers'])}\nRoute: {r.get('route_type') or r.get('polyretro_class')}"
            f"\nGive the polymer repeat unit in PSMILES.", {"repeat_unit": rec["ru"]})))
    for e in out:
        e["polymer_key"] = rec["key"]
        e["family"] = rec["family"]
    return out


def _decoy(rec, kb_list, rng):
    """Monomer set from another polymer that fails the mass balance against this RU."""
    for _ in range(20):
        other = rng.choice(kb_list)
        if other["key"] == rec["key"]:
            continue
        rs = [r for r in other["routes"] if _good(r)]
        if not rs:
            continue
        r = rng.choice(rs)
        byp, _ = infer_byproduct(rec["ru"], r["monomers"])
        if byp is None:
            return list(r["monomers"])
    return None


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", nargs="*", default=[])
    ap.add_argument("--out", default="data_out")
    ap.add_argument("--n-aug", type=int, default=4, help="target SMILES variants per polymer (train only)")
    ap.add_argument("--test-families", nargs="*", default=["polyamide", "fluoropolymer", "paek"])
    ap.add_argument("--val-frac", type=float, default=0.1)
    ap.add_argument("--test-frac", type=float, default=0.1)
    ap.add_argument("--keep-unbalanced", action="store_true",
                    help="keep polyRETRO routes whose by-product could not be inferred")
    ap.add_argument("--seed", type=int, default=13)
    a = ap.parse_args()

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    seed, qc_seed = load_seed()
    pr, qc_pr = load_polyretro(a.csv)
    if a.keep_unbalanced:
        for r in pr:
            if r["qc"].startswith("UNBALANCED"):
                r["qc"] = "ok"
    kb = merge(seed, pr)

    with open(out / "qc_report.csv", "w", newline="") as fh:
        rows = qc_seed + qc_pr
        cols = sorted({k for r in rows for k in r})
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    with open(out / "routes_kb.jsonl", "w") as fh:
        for rec in kb.values():
            fh.write(json.dumps(rec, ensure_ascii=False) + "\n")

    assign = split(kb, a.test_families, a.val_frac, a.test_frac, a.seed)
    rng = random.Random(a.seed)
    kb_list = list(kb.values())
    buckets = defaultdict(list)
    for rec in kb_list:
        sp = assign[rec["family"]]
        n_aug = a.n_aug if sp == "train" else 1  # never augment eval targets
        buckets[sp].extend(examples_for(rec, kb_list, n_aug, rng))
    for sp, exs in buckets.items():
        rng.shuffle(exs)
        with open(out / f"sft_{sp}.jsonl", "w") as fh:
            for e in exs:
                fh.write(json.dumps(e, ensure_ascii=False) + "\n")

    manifest = defaultdict(list)
    for rec in kb_list:
        manifest[assign[rec["family"]]].append(rec.get("abbr") or rec["ru"])
    json.dump(manifest, open(out / "split_manifest.json", "w"), indent=1)

    # ---- summary
    n_bad = sum(1 for r in qc_seed + qc_pr if r["qc"] != "ok")
    print(f"polymers in KB: {len(kb)}  (curated {len(seed)}, polyRETRO-only {len(kb) - len(seed)})")
    print(f"routes: {sum(len(r['routes']) for r in kb_list)}   QC failures: {n_bad} (see qc_report.csv)")
    multi = sum(1 for r in kb_list if sum(_good(x) for x in r["routes"]) > 1)
    print(f"polymers with >1 valid route: {multi}")
    for sp in ("train", "val", "test"):
        tasks = Counter(e["task"] for e in buckets.get(sp, []))
        print(f"{sp:5s}: {len(buckets.get(sp, []))} examples  {dict(tasks)}")
    print("test polymers:", ", ".join(manifest["test"]))


if __name__ == "__main__":
    main()
