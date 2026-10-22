"""Example: validate a reaction and report its energies.

Two modes:

1. Default, addition of HCl to propene. Each candidate is a protonation site, and the script
   reports ΔG, the barrier Ea from a restrained scan, and regioselectivity.

2. A row from a polymerization dataset (--csv FILE --pid PID). The reaction is
   reactant1 + reactant2 -> Polymer_smiles, with the by-product found by atom balance.
   The script reports the gate, ΔG and ΔH per repeat unit and the validator verdict. Barrier and
   regioselectivity are not computed: for symmetric step-growth monomers every reacting group
   is equivalent, so there is no site choice to compare.

The chemistry lives in the validator (RouteValidator.evaluate_candidates and validate_row). This
script only picks the example data and prints the results.

Usage:
    python example_reaction.py                                  # propene + HCl, needs MISTRAL_API_KEY
    python example_reaction.py --no-llm                         # energies only
    python example_reaction.py --csv data.csv --pid PT00001     # one dataset row
"""

import argparse

from rich.console import Console
from rich.table import Table
from rich.text import Text

from config import load_config
from dataset import find_row
from validator import RouteValidator, format_deltas

REACTANTS = "C=CC.Cl"  # propene + HCl
DONOR = "Cl"  # the element that carries the transferred proton (HCl)

CANDIDATES = [
    {"name": "protonate C1 (Markovnikov)", "site": 0, "product": "CC(C)Cl"},
    {"name": "protonate C2 (anti-Markovnikov)", "site": 1, "product": "CCCCl"},
]


def print_table(results):
    table = Table(title="Possible reactions: HCl + propene", title_style="bold")
    table.add_column("Candidate", style="bold")
    table.add_column("Product")
    table.add_column("ΔG", justify="right")
    table.add_column("Ea", justify="right")
    table.add_column("TS (Å)", justify="right")
    table.add_column("Ratio", justify="right")
    table.add_column("Valid", justify="center")
    for r in results:
        ea = f"{r['Ea']:.2f}" + (" *" if r["at_end"] else "")
        valid = Text("yes", style="green") if r["valid"] else Text("no", style="red")
        table.add_row(Text(r["name"]), Text(r["product"]), f"{r['dG']:.2f}", ea,
                      f"{r['ts_distance']:.2f}", f"{100 * r['fraction']:.1f}%", valid)
    console = Console()
    console.print(table)
    console.print("ΔG and Ea in kcal/mol. Ratio = Boltzmann fraction of product at 298 K from the barriers.", style="dim")
    if any(r["at_end"] for r in results):
        console.print("* The barrier maximum is at the end of the scan, so the true barrier may be higher.", style="yellow")


def print_row(row, verdict, use_llm):
    print(f"PID {row['PID']}: {row['reactant1_name']} + {row['reactant2_name']} -> {row['common_name'] or row['Polymer_smiles']}")
    print(f"  product (PSMILES): {row['Polymer_smiles']}")
    if use_llm:
        print(f"  gate: {'ok' if verdict['gate_ok'] else 'FAIL'} ({verdict['gate_message']})")
        print(f"  {format_deltas(verdict['deltas'])} per repeat unit" if verdict["deltas"]
              else f"  energies not computed: {verdict['note']}")
        print(f"  validator: valid={verdict['valid']}")
        print(f"  {verdict['rationale']}")
    else:
        favourable = verdict["valid"]
        print(f"  {format_deltas(verdict['deltas'])} per repeat unit ({'favourable' if favourable else 'not favourable'})")
    print("  barrier and regioselectivity: not computed for this step-growth reaction (see module docstring)")


def main():
    parser = argparse.ArgumentParser(description="Validate a reaction: ΔG, and Ea and regioselectivity for propene + HCl")
    parser.add_argument("--no-llm", action="store_true", help="skip the LLM review, energies only")
    parser.add_argument("--csv", help="dataset CSV; with --pid, validates one row instead of the propene example")
    parser.add_argument("--pid", help="PID of the row to validate")
    args = parser.parse_args()

    config = load_config()
    if config.validator.backend != "gfn2":
        raise SystemExit("this example needs backend = gfn2 in [validator]")
    validator = RouteValidator(config.validator)

    if args.csv:
        if not args.pid:
            raise SystemExit("--csv needs --pid")
        try:
            row = find_row(args.csv, args.pid)
        except KeyError as exc:
            raise SystemExit(str(exc))
        print_row(row, validator.validate_row(row, use_llm=not args.no_llm), use_llm=not args.no_llm)
        return
    results = validator.evaluate_candidates(REACTANTS, CANDIDATES, DONOR, use_llm=not args.no_llm)
    print_table(results)


if __name__ == "__main__":
    main()
