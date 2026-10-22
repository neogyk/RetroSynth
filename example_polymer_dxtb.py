"""Validate one polymerization reaction from a dataset row with the dxtb workflow.

Same reaction as example_polymer.py: reactant1 + reactant2 -> Polymer_smiles, with the by-product
found by atom balance. The energies come from dxtb (PyTorch GFN xTB):
- each molecule is relaxed with dxtb's autograd gradients
- the Hessian comes from dxtb, and it gives the enthalpy H and the free energy G
- the report gives ΔE (electronic), ΔH and ΔG (per repeat unit, kcal/mol)

GFN2 in dxtb needs the libcint library. Use --method GFN1 on a machine without libcint.

Usage:
    python example_polymer_dxtb.py --csv polymer_templates.csv --pid PT00001 --method GFN1 --no-llm
    python example_polymer_dxtb.py --csv polymer_templates.csv --pid PT00001     # GFN2 with the LLM review
"""

import argparse
import dataclasses

from config import load_config
from dataset import find_row, row_reactants
from molecules import TreeNode
from validator import RouteValidator, check_step_gate


def main():
    parser = argparse.ArgumentParser(description="Validate a polymerization reaction with dxtb: dE, dH, dG")
    parser.add_argument("--csv", required=True, help="dataset CSV")
    parser.add_argument("--pid", required=True, help="PID of the row to validate")
    parser.add_argument("--method", default="GFN2", choices=["GFN1", "GFN2"], help="dxtb method (GFN2 needs libcint)")
    parser.add_argument("--no-llm", action="store_true", help="skip the LLM review, energies only")
    args = parser.parse_args()

    settings = load_config().validator
    config = dataclasses.replace(settings, backend="dxtb", dxtb_method=args.method, criterion="de, dh, dg")
    try:
        row = find_row(args.csv, args.pid)
    except KeyError as exc:
        raise SystemExit(str(exc))
    counts = row_reactants(row)
    route = TreeNode(row["Polymer_smiles"], 1.0, [TreeNode(s) for s in counts for _ in range(counts[s])])
    validator = RouteValidator(config)

    print(f"PID {row['PID']}: {row['reactant1_name']} + {row['reactant2_name']} -> {row['common_name'] or row['Polymer_smiles']}")
    print(f"  repeat unit (PSMILES): {row['Polymer_smiles']}")
    print(f"  backend: dxtb {args.method}")

    if args.no_llm:
        gate_ok, _, gate_message = check_step_gate(route.smiles, counts)
        print(f"  gate: {'ok' if gate_ok else 'FAIL'} ({gate_message})")
        if not gate_ok:
            return
        detail = validator.step_deltas(route.smiles, counts)
        for name in ("dE", "dH", "dG"):
            if name in detail:
                print(f"  {name} = {detail[name]:.2f} kcal/mol per repeat unit")
            else:
                print(f"  {name} not computed")
        return

    verdict = validator(route)
    step = verdict.steps[0]
    print(f"  gate: {'ok' if step.gate_ok else 'FAIL'} ({step.message})")
    for name in ("dE", "dH", "dG"):
        if name in step.deltas:
            print(f"  {name} = {step.deltas[name]:.2f} kcal/mol per repeat unit")
        else:
            print(f"  {name} not computed")
    print(f"  validator: valid={verdict.valid}, favourable={verdict.favourable}")
    print(f"  {verdict.rationale}")


if __name__ == "__main__":
    main()
