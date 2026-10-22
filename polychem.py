"""
polychem.py - RDKit utilities for polymer repeat units (PSMILES with [*] end points).

Functions
---------
ru_formula(smi)            element Counter of a repeat unit, ignoring [*]
mol_formula(smi)           element Counter of a small molecule
canonical_ru_key(smi)      frame- and direction-invariant key for a repeat unit
                           ([*]CCO[*] == [*]OCC[*] == [*]COC[*] ...)
mass_balance(...)          checks  sum(v_i * monomer_i) == n_ru * RU + sum(v_j * byproduct_j)
infer_byproduct(...)       guesses the small-molecule by-product from the atom difference
randomized_smiles(smi, n)  RDKit random SMILES (works with dummy atoms)
frame_shifts(smi)          all equivalent repeat-unit frames obtained by moving the
                           cut along the acyclic backbone ([*]OCC[*] -> [*]COC[*], [*]CCO[*])
"""
from __future__ import annotations

import random
from collections import Counter

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


# --------------------------------------------------------------------------- #
# basic parsing / formulae
# --------------------------------------------------------------------------- #
def _mol(smi: str) -> Chem.Mol:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise ValueError(f"RDKit could not parse SMILES: {smi}")
    return m


def is_valid(smi: str) -> bool:
    try:
        return Chem.MolFromSmiles(smi) is not None
    except Exception:
        return False


def to_psmiles(s: str) -> str:
    """RDKit writes dummies as '*'; normalise to the bracketed PSMILES form '[*]'."""
    return s.replace("[*]", "*").replace("*", "[*]")


def canonical(smi: str) -> str:
    return to_psmiles(Chem.MolToSmiles(_mol(smi)))


def mol_formula(smi: str) -> Counter:
    m = _mol(smi)
    c = Counter()
    for a in m.GetAtoms():
        if a.GetAtomicNum() == 0:
            continue
        c[a.GetSymbol()] += 1
        c["H"] += a.GetTotalNumHs()
    return Counter({k: v for k, v in c.items() if v})


ru_formula = mol_formula  # dummy atoms are skipped, so the same code works


def formula_str(c: Counter) -> str:
    order = ["C", "H"] + sorted(k for k in c if k not in ("C", "H"))
    out = []
    for k in order:
        v = c.get(k, 0)
        if abs(v) < 1e-9:
            continue
        v = round(v, 3)
        v = int(v) if float(v).is_integer() else v
        out.append(f"{k}{'' if v == 1 else v}")
    return "".join(out) or "(none)"


def star_count(smi: str) -> int:
    return sum(1 for a in _mol(smi).GetAtoms() if a.GetAtomicNum() == 0)


# --------------------------------------------------------------------------- #
# canonical, frame-invariant repeat-unit key
# --------------------------------------------------------------------------- #
def _star_info(m: Chem.Mol):
    stars = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() == 0]
    if len(stars) != 2:
        raise ValueError("repeat unit must contain exactly two [*]")
    nbrs, btypes = [], []
    for s in stars:
        a = m.GetAtomWithIdx(s)
        if a.GetDegree() != 1:
            raise ValueError("each [*] must have exactly one neighbour")
        b = a.GetBonds()[0]
        nbrs.append(b.GetOtherAtomIdx(s))
        btypes.append(b.GetBondType())
    return stars, nbrs, btypes


def canonical_ru_key(smi: str) -> str:
    """Cyclic-dimer key: join two copies head-to-tail into a ring and canonicalise.
    Independent of where the repeat unit was 'cut' and of its writing direction."""
    m = _mol(smi)
    try:
        stars, nbrs, btypes = _star_info(m)
    except ValueError:
        return canonical(smi)  # branched / network / unusual -> fall back
    n = m.GetNumAtoms()
    combo = Chem.RWMol(Chem.CombineMols(m, m))
    h0, t0 = nbrs
    h1, t1 = nbrs[0] + n, nbrs[1] + n
    combo.AddBond(t0, h1, btypes[1])
    combo.AddBond(t1, h0, btypes[0])
    for s in sorted([stars[0], stars[1], stars[0] + n, stars[1] + n], reverse=True):
        combo.RemoveAtom(s)
    out = combo.GetMol()
    try:
        Chem.SanitizeMol(out)
        return "RUKEY:" + Chem.MolToSmiles(out)
    except Exception:
        return canonical(smi)


# --------------------------------------------------------------------------- #
# mass balance
# --------------------------------------------------------------------------- #
def _sum(species: dict[str, float]) -> Counter:
    tot = Counter()
    for smi, coef in (species or {}).items():
        for k, v in mol_formula(smi).items():
            tot[k] += coef * v
    return tot


def mass_balance(ru: str, reactants: dict[str, float], byproducts: dict[str, float] | None,
                 n_ru: float = 1.0) -> tuple[bool, Counter]:
    """Returns (balanced?, residual) where residual = reactants - n_ru*RU - byproducts."""
    lhs = _sum(reactants)
    rhs = Counter({k: n_ru * v for k, v in ru_formula(ru).items()})
    for k, v in _sum(byproducts).items():
        rhs[k] += v
    keys = set(lhs) | set(rhs)
    resid = Counter({k: lhs.get(k, 0) - rhs.get(k, 0) for k in keys})
    resid = Counter({k: v for k, v in resid.items() if abs(v) > 1e-6})
    return (len(resid) == 0), resid


# small molecules that are plausible (or formally acceptable) condensation by-products
KNOWN_BYPRODUCTS = {
    "O": "water", "Cl": "HCl", "F": "HF", "Br": "HBr", "I": "HI", "CO": "methanol",
    "CCO": "ethanol", "CC(=O)O": "acetic acid", "Oc1ccccc1": "phenol", "[HH]": "H2 (formal)",
    "ClCl": "Cl2 (formal, dehalogenative coupling)", "BrBr": "Br2 (formal, dehalogenative coupling)",
    "[Na]Cl": "NaCl", "N": "ammonia", "O=C=O": "CO2", "OCCO": "ethylene glycol",
}


def infer_byproduct(ru: str, reactants: dict[str, float]) -> tuple[dict[str, float] | None, Counter]:
    """Try to explain (reactants - RU) as k * one known by-product (k = 1..4)."""
    _, diff = mass_balance(ru, reactants, None)
    if not diff:
        return {}, diff
    for smi in KNOWN_BYPRODUCTS:
        f = mol_formula(smi)
        if set(f) != set(diff):
            continue
        ks = {diff[e] / f[e] for e in f}
        if len(ks) == 1:
            k = ks.pop()
            if k > 0 and abs(k - round(k)) < 1e-6 and round(k) <= 4:
                return {smi: float(round(k))}, Counter()
    return None, diff


# --------------------------------------------------------------------------- #
# augmentation
# --------------------------------------------------------------------------- #
def randomized_smiles(smi: str, n: int = 5, seed: int = 0) -> list[str]:
    m = _mol(smi)
    rng = random.Random(seed)
    out = set()
    for _ in range(n * 10):
        s = Chem.MolToSmiles(m, doRandom=True, canonical=False)
        # keep it compact: RDKit writes '*' for dummies, convert to [*] for PSMILES style
        s = to_psmiles(s)
        out.add(s)
        if len(out) >= n:
            break
    out = sorted(out)
    rng.shuffle(out)
    return out[:n]


def frame_shifts(smi: str) -> list[str]:
    """Equivalent repeat units obtained by cutting a different acyclic single backbone bond."""
    m = _mol(smi)
    try:
        stars, nbrs, btypes = _star_info(m)
    except ValueError:
        return [canonical(smi)]
    h, t = nbrs
    if h == t or m.GetBondBetweenAtoms(h, t) is not None:
        return [canonical(smi)]
    orig_ring = {b.GetIdx() for b in m.GetBonds() if b.IsInRing()}
    rw = Chem.RWMol(m)
    rw.AddBond(h, t, Chem.BondType.SINGLE)
    for s in sorted(stars, reverse=True):
        rw.RemoveAtom(s)
    closed = rw.GetMol()
    try:
        Chem.SanitizeMol(closed)
    except Exception:
        return [canonical(smi)]
    # map: atoms keep order except the removed stars; recompute ring bonds on the closed mol
    results = set()
    for b in closed.GetBonds():
        if b.GetBondType() != Chem.BondType.SINGLE or not b.IsInRing():
            continue
        # skip bonds that belong to a ring that existed before (e.g. benzene)
        a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        ri = closed.GetRingInfo()
        in_small_ring = any(len(r) < closed.GetNumAtoms() and a1 in r and a2 in r and
                            _ring_was_original(r, closed, m, stars)
                            for r in ri.AtomRings())
        if in_small_ring:
            continue
        cut = Chem.RWMol(closed)
        cut.RemoveBond(a1, a2)
        s1 = cut.AddAtom(Chem.Atom(0))
        s2 = cut.AddAtom(Chem.Atom(0))
        cut.AddBond(a1, s1, Chem.BondType.SINGLE)
        cut.AddBond(a2, s2, Chem.BondType.SINGLE)
        new = cut.GetMol()
        try:
            Chem.SanitizeMol(new)
            if new.GetRingInfo().NumRings() != m.GetRingInfo().NumRings():
                continue
            results.add(to_psmiles(Chem.MolToSmiles(new)))
        except Exception:
            continue
    return sorted(results) or [canonical(smi)]


def _ring_was_original(ring, closed, orig, stars) -> bool:
    """True if this ring (atom indices in `closed`) already existed in the original RU."""
    # index shift: atoms with idx > star idx moved down by one per removed star
    def to_orig(i):
        for s in sorted(stars):
            if i >= s:
                i += 1
        return i
    orig_atoms = {to_orig(i) for i in ring}
    return any(set(r) == orig_atoms for r in orig.GetRingInfo().AtomRings())
