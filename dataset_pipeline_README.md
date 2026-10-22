# Multi-route polymer retrosynthesis dataset

Builds a fine-tuning set where each target polymer carries **all** of its known routes
(free-radical, Ziegler-Natta, metallocene, anionic, ROP, polycondensation, interfacial,
enzymatic, biological, post-polymerisation modification), plus rejected disconnections.

```
pip install rdkit pandas
python build_dataset.py --csv train_condensation.csv train_addition.csv train_2k.csv \
       --out data_out --n-aug 4 --test-families polyamide fluoropolymer paek
```

## Files
| file | role |
|---|---|
| `polychem.py` | RDKit helpers: frame-invariant repeat-unit key, mass balance, by-product inference, SMILES + frame-shift augmentation |
| `seed_commercial.py` | curated commercial polymers (42) with ranked multi-route entries; every route is mass-balance checked |
| `build_dataset.py` | ingests polyRETRO-format CSVs, merges with the seed, QC, family-level split, writes chat JSONL |

## Design decisions
1. **One record per polymer, many routes.** Rows are grouped by a cyclic-dimer key, so
   `[*]CCO[*]`, `[*]OCC[*]` and `[*]COC[*]` are recognised as the same polymer.
2. **Mass balance is the gatekeeper.** `Σ ν·monomers = n·RU + Σ ν·by-products` must hold
   (by-products in formal form: H2O, HCl, HF, MeOH, PhOH, NaCl...). polyRETRO rows get
   their by-product inferred; rows that cannot be balanced are flagged in `qc_report.csv`.
3. **Negatives are explicit.** `scale: "infeasible"` routes (radical PP, vinyl alcohol,
   EG dehydration to PEO) and mass-balance decoys teach the model to reject routes.
4. **Split by polymer family**, never by row, so alternative routes of a test polymer
   never leak into training. Augmentation (random SMILES + repeat-unit frame shifts) is
   applied to train only.
5. **Tasks** (field `task`): `all_routes`, `route_conditioned`, `class`, `feasibility`, `forward`.

## Before release
- Fill `refs` (DOI/patent) for every curated route; conditions are handbook-typical ranges.
- No yields or rates are included. Add them only from sourced data.
