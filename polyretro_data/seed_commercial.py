"""
seed_commercial.py - hand-curated, multi-route knowledge base for commercial polymers.

Every route is machine-checked by build_dataset.py:
  * all SMILES must parse in RDKit
  * the route must be mass balanced:
        sum(v_i * monomer_i) + sum(v_j * coreactant_j) == ru_per_set * RU + sum(v_k * byproduct_k)
    (by-products such as HCl / HF / NaCl are written in their *formal* form)
  * precursor_steps (e.g. enzymatic monomer synthesis) are balanced as small-molecule reactions

Field meanings
--------------
mechanism_class : chain_growth | step_growth_condensation | step_growth_polyaddition |
                  ring_opening | romp | oxidative_coupling | biological
polyretro_class : the 4-class label used by polyRETRO (Condensation | Addition |
                  Ringopening | ROMP) so this seed can be merged with polyRETRO-format data
route_type      : the "how" a chemist/engineer would name it - free_radical, ziegler_natta,
                  metallocene, phillips, anionic, cationic, coordination, melt_polycondensation,
                  interfacial, solution_polycondensation, snar_polycondensation, enzymatic,
                  enzymatic_initiation, chemoenzymatic, biological, post_polymerization_modification
process         : bulk | solution | suspension | emulsion | slurry | gas_phase | melt |
                  interfacial | in_vivo | high_pressure
scale           : industrial | pilot | lab | infeasible
rank            : 1 = dominant industrial route for this polymer (ties allowed)

IMPORTANT: conditions here are typical literature/handbook ranges, written as text.
No yields or rates are given - add them only from sourced data (put the DOI in `refs`).
`refs` are deliberately left for you to fill (DOI / patent) before release.
"""

H2O, HCl, HF, MEOH, ETOH, PHOH, O2, H2O2 = "O", "Cl", "F", "CO", "CCO", "Oc1ccccc1", "O=O", "OO"


def R(route_id, mechanism_class, polyretro_class, route_type, process, monomers, *,
      coreactants=None, byproducts=None, ru_per_set=1, catalyst="", conditions="",
      scale="industrial", rank=1, notes="", precursor_steps=None, refs=None):
    return dict(route_id=route_id, mechanism_class=mechanism_class,
                polyretro_class=polyretro_class, route_type=route_type, process=process,
                monomers=monomers, coreactants=coreactants or {}, byproducts=byproducts or {},
                ru_per_set=ru_per_set, catalyst=catalyst, conditions=conditions, scale=scale,
                rank=rank, notes=notes, precursor_steps=precursor_steps or [], refs=refs or [])


def P(name, abbr, family, ru, routes, aliases=()):
    return dict(name=name, abbr=abbr, family=family, ru=ru, routes=routes, aliases=list(aliases))


ETHYLENE, PROPYLENE, STYRENE = "C=C", "C=CC", "C=Cc1ccccc1"

SEED = [
    # ------------------------------------------------------------------ polyolefins
    P("Polyethylene", "PE", "polyolefin", "[*]CC[*]", [
        R("PE-1", "chain_growth", "Addition", "ziegler_natta", "slurry",
          {ETHYLENE: 1}, catalyst="TiCl4/MgCl2 + AlEt3 (HDPE)",
          conditions="70-110 °C, 10-40 bar, H2 for MW control; slurry or gas-phase"),
        R("PE-2", "chain_growth", "Addition", "free_radical", "high_pressure",
          {ETHYLENE: 1}, catalyst="O2 / organic peroxides (LDPE)",
          conditions="150-300 °C, 1500-3000 bar, tubular or autoclave reactor",
          notes="Gives long- and short-chain branched LDPE"),
        R("PE-3", "chain_growth", "Addition", "phillips", "slurry",
          {ETHYLENE: 1}, catalyst="CrOx on silica (Phillips)",
          conditions="85-110 °C, loop slurry reactor in isobutane"),
        R("PE-4", "chain_growth", "Addition", "metallocene", "gas_phase",
          {ETHYLENE: 1}, catalyst="Cp2ZrCl2-type metallocene / MAO", rank=2,
          conditions="70-100 °C, gas-phase or solution; usually with 1-butene/1-hexene/1-octene for LLDPE",
          notes="Single-site catalyst, narrow MWD"),
    ], aliases=["HDPE", "LDPE", "polyethene"]),

    P("Polypropylene (isotactic)", "PP", "polyolefin", "[*]CC([*])C", [
        R("PP-1", "chain_growth", "Addition", "ziegler_natta", "bulk",
          {PROPYLENE: 1}, catalyst="TiCl4/MgCl2 + AlEt3 + internal/external donors",
          conditions="60-80 °C, liquid-propylene bulk (loop) or gas-phase, H2 for MW control"),
        R("PP-2", "chain_growth", "Addition", "metallocene", "bulk",
          {PROPYLENE: 1}, catalyst="rac-ansa-zirconocene (e.g. rac-Me2Si(2-Me-4-Ph-Ind)2ZrCl2) / MAO",
          conditions="50-70 °C", rank=2),
        R("PP-X", "chain_growth", "Addition", "free_radical", "bulk",
          {PROPYLENE: 1}, scale="infeasible", rank=99,
          notes="Radical polymerisation of propylene gives only oligomers: degradative chain "
                "transfer to the allylic C-H stabilises the radical."),
    ], aliases=["polypropene", "iPP"]),

    P("Polystyrene", "PS", "styrenic", "[*]CC([*])c1ccccc1", [
        R("PS-1", "chain_growth", "Addition", "free_radical", "bulk",
          {STYRENE: 1}, catalyst="thermal self-initiation (Mayo) or peroxide",
          conditions="100-180 °C, continuous mass process with ethylbenzene diluent"),
        R("PS-2", "chain_growth", "Addition", "free_radical", "suspension",
          {STYRENE: 1}, catalyst="benzoyl peroxide, PVA/TCP suspending agent",
          conditions="80-130 °C, aqueous suspension (used for EPS beads)", rank=2),
        R("PS-3", "chain_growth", "Addition", "anionic", "solution",
          {STYRENE: 1}, catalyst="sec-BuLi or n-BuLi", scale="industrial", rank=3,
          conditions="40-80 °C in cyclohexane; living, narrow MWD (mainly for SBS block copolymers)"),
        R("PS-4", "chain_growth", "Addition", "metallocene", "bulk",
          {STYRENE: 1}, catalyst="CpTiCl3 or Cp*Ti(OMe)3 / MAO", rank=3,
          conditions="50-90 °C", notes="Gives syndiotactic PS (sPS)"),
    ]),

    P("Poly(vinyl chloride)", "PVC", "vinyl", "[*]CC([*])Cl", [
        R("PVC-1", "chain_growth", "Addition", "free_radical", "suspension",
          {"C=CCl": 1}, catalyst="peroxydicarbonates / peresters; PVA suspending agent",
          conditions="50-70 °C (T sets K-value), aqueous suspension, ~80-90 % conversion"),
        R("PVC-2", "chain_growth", "Addition", "free_radical", "emulsion",
          {"C=CCl": 1}, catalyst="persulfate, anionic surfactant", rank=2,
          conditions="40-60 °C, emulsion/microsuspension (paste PVC)"),
    ]),

    P("Poly(methyl methacrylate)", "PMMA", "acrylic", "[*]CC([*])(C)C(=O)OC", [
        R("PMMA-1", "chain_growth", "Addition", "free_radical", "bulk",
          {"C=C(C)C(=O)OC": 1}, catalyst="AIBN or peroxides",
          conditions="cell-cast or continuous bulk, 40-160 °C"),
        R("PMMA-2", "chain_growth", "Addition", "free_radical", "suspension",
          {"C=C(C)C(=O)OC": 1}, catalyst="peroxide, PVA stabiliser", rank=2,
          conditions="70-90 °C, bead (moulding) grades"),
        R("PMMA-3", "chain_growth", "Addition", "anionic", "solution",
          {"C=C(C)C(=O)OC": 1}, catalyst="sec-BuLi/LiCl or t-BuMgBr", scale="lab", rank=4,
          conditions="-78 °C in THF; living, tacticity control"),
        R("PMMA-4", "chain_growth", "Addition", "enzymatic_initiation", "solution",
          {"C=C(C)C(=O)OC": 1}, catalyst="horseradish peroxidase (HRP) / H2O2 / 2,4-pentanedione",
          scale="lab", rank=5, conditions="aqueous/organic, ~25-60 °C",
          notes="Enzyme generates the initiating radical; propagation is ordinary radical chemistry"),
    ], aliases=["acrylic glass", "Plexiglas"]),

    P("Polyacrylonitrile", "PAN", "acrylic", "[*]CC([*])C#N", [
        R("PAN-1", "chain_growth", "Addition", "free_radical", "suspension",
          {"C=CC#N": 1}, catalyst="persulfate/bisulfite redox",
          conditions="40-60 °C aqueous slurry (precipitation) polymerisation"),
        R("PAN-2", "chain_growth", "Addition", "free_radical", "solution",
          {"C=CC#N": 1}, catalyst="AIBN", rank=1,
          conditions="DMSO, DMF or aqueous NaSCN, 50-70 °C; dope spun directly to fibre"),
    ]),

    P("Polytetrafluoroethylene", "PTFE", "fluoropolymer", "[*]C(F)(F)C([*])(F)F", [
        R("PTFE-1", "chain_growth", "Addition", "free_radical", "emulsion",
          {"FC(F)=C(F)F": 1}, catalyst="ammonium persulfate, fluorosurfactant",
          conditions="10-30 bar, 60-90 °C, aqueous dispersion"),
        R("PTFE-2", "chain_growth", "Addition", "free_radical", "suspension",
          {"FC(F)=C(F)F": 1}, catalyst="persulfate, little/no surfactant", rank=1,
          conditions="granular PTFE, vigorous agitation"),
    ], aliases=["Teflon"]),

    P("Poly(vinylidene fluoride)", "PVDF", "fluoropolymer", "[*]CC([*])(F)F", [
        R("PVDF-1", "chain_growth", "Addition", "free_radical", "emulsion",
          {"C=C(F)F": 1}, catalyst="persulfate or di-tert-butyl peroxide",
          conditions="75-90 °C, 20-35 bar"),
        R("PVDF-2", "chain_growth", "Addition", "free_radical", "suspension",
          {"C=C(F)F": 1}, catalyst="peroxydicarbonate", conditions="30-60 °C", rank=1),
    ]),

    P("Poly(vinyl acetate)", "PVAc", "vinyl", "[*]CC([*])OC(C)=O", [
        R("PVAc-1", "chain_growth", "Addition", "free_radical", "emulsion",
          {"C=COC(C)=O": 1}, catalyst="persulfate, PVA protective colloid",
          conditions="60-80 °C, aqueous emulsion (adhesives)"),
        R("PVAc-2", "chain_growth", "Addition", "free_radical", "solution",
          {"C=COC(C)=O": 1}, catalyst="AIBN or peroxide", conditions="methanol, 60-70 °C",
          notes="Precursor for PVA"),
    ]),

    P("Poly(vinyl alcohol)", "PVA", "vinyl", "[*]CC([*])O", [
        R("PVA-1", "chain_growth", "Addition", "post_polymerization_modification", "solution",
          {"C=COC(C)=O": 1}, coreactants={MEOH: 1}, byproducts={"COC(C)=O": 1},
          catalyst="step 1: radical (AIBN) in MeOH; step 2: NaOH-catalysed methanolysis",
          conditions="step 1: 60-70 °C; step 2: 20-60 °C, degree of hydrolysis set by NaOH/time",
          notes="Two steps: vinyl acetate -> PVAc -> PVA (+ methyl acetate)"),
        R("PVA-X", "chain_growth", "Addition", "free_radical", "bulk",
          {"C=CO": 1}, scale="infeasible", rank=99,
          notes="Vinyl alcohol is not an isolable monomer: it tautomerises to acetaldehyde."),
    ]),

    P("Polyacrylamide", "PAM", "acrylic", "[*]CC([*])C(N)=O", [
        R("PAM-1", "chain_growth", "Addition", "free_radical", "solution",
          {"C=CC(N)=O": 1}, catalyst="persulfate/bisulfite redox or azo initiator",
          conditions="aqueous gel or inverse-emulsion polymerisation, 10-60 °C"),
        R("PAM-2", "chain_growth", "Addition", "chemoenzymatic", "solution",
          {"C=CC(N)=O": 1}, catalyst="step 1: nitrile hydratase (Rhodococcus) whole cells; step 2: redox radical",
          conditions="step 1: ~0-15 °C aqueous; step 2 as PAM-1",
          notes="Industrial biocatalytic acrylamide (Mitsubishi/Nitto process) feeding radical polymerisation",
          precursor_steps=[dict(reactants={"C=CC#N": 1, H2O: 1}, products={"C=CC(N)=O": 1},
                                catalyst="nitrile hydratase", type="enzymatic")]),
    ]),

    P("Poly(acrylic acid)", "PAA", "acrylic", "[*]CC([*])C(=O)O", [
        R("PAA-1", "chain_growth", "Addition", "free_radical", "solution",
          {"C=CC(=O)O": 1}, catalyst="persulfate", conditions="aqueous, often partially neutralised, 50-90 °C"),
    ]),

    # ------------------------------------------------------------------ polyethers / acetals
    P("Polyoxymethylene", "POM", "polyether", "[*]OC[*]", [
        R("POM-1", "ring_opening", "Ringopening", "cationic", "bulk",
          {"C1OCOCO1": 1}, ru_per_set=3, catalyst="BF3·OEt2 (+ 1,3-dioxolane comonomer for copolymer)",
          conditions="60-100 °C, bulk; end-capping/stabilisation afterwards"),
        R("POM-2", "chain_growth", "Addition", "anionic", "slurry",
          {"C=O": 1}, catalyst="tertiary amines / alkoxides", rank=2,
          conditions="anhydrous formaldehyde in hydrocarbon slurry; acetic-anhydride end-capping"),
    ], aliases=["acetal", "polyacetal", "Delrin"]),

    P("Poly(ethylene oxide)", "PEO", "polyether", "[*]CCO[*]", [
        R("PEO-1", "ring_opening", "Ringopening", "anionic", "bulk",
          {"C1CO1": 1}, catalyst="KOH or alkoxide initiator (diol/alcohol starter)",
          conditions="100-150 °C, semi-batch EO feed (PEG grades)"),
        R("PEO-2", "ring_opening", "Ringopening", "coordination", "slurry",
          {"C1CO1": 1}, catalyst="Ca amide/alkoxide or Zn/Al-based catalysts", rank=2,
          conditions="hydrocarbon slurry; high-MW PEO (>10^5 g/mol)"),
        R("PEO-X", "step_growth_condensation", "Condensation", "melt_polycondensation", "melt",
          {"OCCO": 1}, byproducts={H2O: 1}, scale="infeasible", rank=99,
          notes="Acid-catalysed dehydration of ethylene glycol gives diethylene glycol and 1,4-dioxane, not high polymer"),
    ], aliases=["PEG", "polyethylene glycol"]),

    # ------------------------------------------------------------------ aliphatic polyesters
    P("Poly(lactic acid)", "PLA", "aliphatic_polyester", "[*]OC(C)C(=O)[*]", [
        R("PLA-1", "ring_opening", "Ringopening", "coordination", "melt",
          {"CC1OC(=O)C(C)OC1=O": 1}, ru_per_set=2, catalyst="Sn(Oct)2 (+ alcohol initiator)",
          conditions="180-210 °C melt, continuous; lactide made by back-biting depolymerisation of PLA oligomer"),
        R("PLA-2", "step_growth_condensation", "Condensation", "solution_polycondensation", "solution",
          {"CC(O)C(=O)O": 1}, byproducts={H2O: 1}, catalyst="Sn powder / SnCl2", rank=2, scale="pilot",
          conditions="azeotropic dehydration in diphenyl ether, ~130-150 °C"),
        R("PLA-3", "ring_opening", "Ringopening", "enzymatic", "bulk",
          {"CC1OC(=O)C(C)OC1=O": 1}, ru_per_set=2, catalyst="lipase (Burkholderia cepacia lipase PS or CALB)",
          scale="lab", rank=4, conditions="80-130 °C bulk", notes="Low to moderate MW; L-lactide is a poor lipase substrate"),
    ], aliases=["polylactide", "PLLA"]),

    P("Poly(ε-caprolactone)", "PCL", "aliphatic_polyester", "[*]OCCCCCC(=O)[*]", [
        R("PCL-1", "ring_opening", "Ringopening", "coordination", "bulk",
          {"O=C1CCCCCO1": 1}, catalyst="Sn(Oct)2 + alcohol/diol initiator",
          conditions="110-160 °C bulk"),
        R("PCL-2", "ring_opening", "Ringopening", "enzymatic", "solution",
          {"O=C1CCCCCO1": 1}, catalyst="Candida antarctica lipase B, immobilised (Novozym 435)",
          scale="lab", rank=2, conditions="toluene, 60-90 °C, dry enzyme",
          notes="Best-established enzymatic ROP; metal-free"),
        R("PCL-3", "step_growth_condensation", "Condensation", "enzymatic", "solution",
          {"OCCCCCC(=O)O": 1}, byproducts={H2O: 1}, catalyst="CALB (Novozym 435)",
          scale="lab", rank=4, conditions="toluene or diphenyl ether, vacuum/molecular sieves"),
    ]),

    P("Poly(glycolic acid)", "PGA", "aliphatic_polyester", "[*]OCC(=O)[*]", [
        R("PGA-1", "ring_opening", "Ringopening", "coordination", "melt",
          {"O=C1COC(=O)CO1": 1}, ru_per_set=2, catalyst="Sn(Oct)2", conditions="220-240 °C melt"),
        R("PGA-2", "step_growth_condensation", "Condensation", "melt_polycondensation", "melt",
          {"OCC(=O)O": 1}, byproducts={H2O: 1}, scale="lab", rank=3,
          conditions="melt/solid-state, gives low-MW PGA"),
    ], aliases=["polyglycolide"]),

    P("Poly(3-hydroxybutyrate)", "PHB", "aliphatic_polyester", "[*]OC(C)CC(=O)[*]", [
        R("PHB-1", "biological", "Condensation", "biological", "in_vivo",
          {"C[C@@H](O)CC(=O)O": 1}, byproducts={H2O: 1},
          catalyst="PHA synthase (PhaC) in Cupriavidus necator or recombinant E. coli",
          conditions="fed-batch fermentation on sugars/oils under nutrient limitation, ~30 °C",
          notes="Formal monomer shown; the real substrate is (R)-3-hydroxybutyryl-CoA (CoA released)"),
        R("PHB-2", "ring_opening", "Ringopening", "coordination", "solution",
          {"CC1CC(=O)O1": 1}, catalyst="Zn/Y/Al complexes (β-butyrolactone ROP)",
          scale="lab", rank=3, conditions="toluene, 20-110 °C; tacticity depends on catalyst"),
    ], aliases=["P3HB", "PHA"]),

    P("Poly(butylene succinate)", "PBS", "aliphatic_polyester", "[*]OCCCCOC(=O)CCC(=O)[*]", [
        R("PBS-1", "step_growth_condensation", "Condensation", "melt_polycondensation", "melt",
          {"O=C(O)CCC(=O)O": 1, "OCCCCO": 1}, byproducts={H2O: 2}, catalyst="Ti(OBu)4",
          conditions="esterification 150-200 °C, then polycondensation 220-240 °C under vacuum"),
        R("PBS-2", "step_growth_condensation", "Condensation", "enzymatic", "solution",
          {"CCOC(=O)CCC(=O)OCC": 1, "OCCCCO": 1}, byproducts={ETOH: 2},
          catalyst="CALB (Novozym 435)", scale="lab", rank=3,
          conditions="diphenyl ether, 80-95 °C, two-stage (oligomerisation then vacuum)"),
    ]),

    # ------------------------------------------------------------------ aromatic polyesters
    P("Poly(ethylene terephthalate)", "PET", "aromatic_polyester",
      "[*]OCCOC(=O)c1ccc(C(=O)[*])cc1", [
        R("PET-1", "step_growth_condensation", "Condensation", "melt_polycondensation", "melt",
          {"OC(=O)c1ccc(C(=O)O)cc1": 1, "OCCO": 1}, byproducts={H2O: 2},
          catalyst="Sb2O3 (or Ti / Ge)",
          conditions="direct esterification 240-260 °C, melt polycondensation 270-290 °C under vacuum, then SSP"),
        R("PET-2", "step_growth_condensation", "Condensation", "melt_polycondensation", "melt",
          {"COC(=O)c1ccc(C(=O)OC)cc1": 1, "OCCO": 1}, byproducts={MEOH: 2}, rank=2,
          catalyst="Zn/Mn acetate (transesterification) then Sb2O3",
          conditions="transesterification 150-220 °C via BHET, then as PET-1"),
    ], aliases=["polyester", "PETE"]),

    P("Poly(butylene terephthalate)", "PBT", "aromatic_polyester",
      "[*]OCCCCOC(=O)c1ccc(C(=O)[*])cc1", [
        R("PBT-1", "step_growth_condensation", "Condensation", "melt_polycondensation", "melt",
          {"OC(=O)c1ccc(C(=O)O)cc1": 1, "OCCCCO": 1}, byproducts={H2O: 2}, catalyst="Ti(OBu)4",
          conditions="230-260 °C; THF formed as side product from BDO"),
        R("PBT-2", "step_growth_condensation", "Condensation", "melt_polycondensation", "melt",
          {"COC(=O)c1ccc(C(=O)OC)cc1": 1, "OCCCCO": 1}, byproducts={MEOH: 2}, catalyst="Ti(OBu)4",
          conditions="transesterification 150-220 °C, polycondensation ~250 °C", rank=1),
    ]),

    P("Poly(ethylene furanoate)", "PEF", "aromatic_polyester", "[*]OCCOC(=O)c1ccc(C(=O)[*])o1", [
        R("PEF-1", "step_growth_condensation", "Condensation", "melt_polycondensation", "melt",
          {"O=C(O)c1ccc(C(=O)O)o1": 1, "OCCO": 1}, byproducts={H2O: 2}, scale="pilot",
          catalyst="Ti or Sn catalysts", conditions="melt polycondensation then SSP (emerging, bio-based)"),
        R("PEF-2", "step_growth_condensation", "Condensation", "melt_polycondensation", "melt",
          {"COC(=O)c1ccc(C(=O)OC)o1": 1, "OCCO": 1}, byproducts={MEOH: 2}, scale="pilot", rank=2,
          catalyst="Ti(OiPr)4", conditions="transesterification from dimethyl FDCA"),
    ]),

    P("Bisphenol A polycarbonate", "PC", "polycarbonate",
      "[*]Oc1ccc(C(C)(C)c2ccc(OC([*])=O)cc2)cc1", [
        R("PC-1", "step_growth_condensation", "Condensation", "interfacial", "interfacial",
          {"CC(C)(c1ccc(O)cc1)c1ccc(O)cc1": 1, "O=C(Cl)Cl": 1}, byproducts={HCl: 2},
          catalyst="NaOH(aq) / CH2Cl2, triethylamine, phenol end-capper",
          conditions="20-40 °C interfacial phosgenation"),
        R("PC-2", "step_growth_condensation", "Condensation", "melt_polycondensation", "melt",
          {"CC(C)(c1ccc(O)cc1)c1ccc(O)cc1": 1, "O=C(Oc1ccccc1)Oc1ccccc1": 1}, byproducts={PHOH: 2},
          catalyst="basic catalysts (e.g. NaOH, tetraalkylammonium/phosphonium salts)",
          conditions="melt transesterification 180-300 °C under falling pressure (phosgene-free)", rank=1),
    ]),

    # ------------------------------------------------------------------ polyamides / aramids / PI
    P("Polyamide 6", "PA6", "polyamide", "[*]NCCCCCC([*])=O", [
        R("PA6-1", "ring_opening", "Ringopening", "hydrolytic", "melt",
          {"O=C1CCCCCN1": 1}, catalyst="water (initiator), acid regulator",
          conditions="240-280 °C, VK tube reactor"),
        R("PA6-2", "ring_opening", "Ringopening", "anionic", "bulk",
          {"O=C1CCCCCN1": 1}, catalyst="sodium caprolactamate + N-acyl-caprolactam activator", rank=2,
          conditions="140-180 °C, minutes; cast nylon / reactive moulding"),
        R("PA6-3", "step_growth_condensation", "Condensation", "melt_polycondensation", "melt",
          {"NCCCCCC(=O)O": 1}, byproducts={H2O: 1}, scale="lab", rank=4,
          conditions="melt polycondensation of 6-aminohexanoic acid"),
    ], aliases=["nylon 6", "polycaprolactam"]),

    P("Polyamide 66", "PA66", "polyamide", "[*]NCCCCCCNC(=O)CCCCC([*])=O", [
        R("PA66-1", "step_growth_condensation", "Condensation", "melt_polycondensation", "melt",
          {"NCCCCCCN": 1, "O=C(O)CCCCC(=O)O": 1}, byproducts={H2O: 2},
          catalyst="none (nylon-6,6 salt, acetic acid end-capper)",
          conditions="aqueous salt concentrated, then 270-290 °C melt"),
        R("PA66-2", "step_growth_condensation", "Condensation", "interfacial", "interfacial",
          {"NCCCCCCN": 1, "O=C(Cl)CCCCC(=O)Cl": 1}, byproducts={HCl: 2}, scale="lab", rank=3,
          catalyst="NaOH(aq) acid acceptor", conditions="room temperature, 'nylon rope trick'"),
    ], aliases=["nylon 66", "nylon 6,6"]),

    P("Polyamide 11", "PA11", "polyamide", "[*]NCCCCCCCCCCC([*])=O", [
        R("PA11-1", "step_growth_condensation", "Condensation", "melt_polycondensation", "melt",
          {"NCCCCCCCCCCC(=O)O": 1}, byproducts={H2O: 1},
          conditions="melt ~220-250 °C; monomer from castor oil (bio-based)"),
    ], aliases=["nylon 11", "Rilsan"]),

    P("Polyamide 12", "PA12", "polyamide", "[*]NCCCCCCCCCCCC([*])=O", [
        R("PA12-1", "ring_opening", "Ringopening", "hydrolytic", "melt",
          {"O=C1CCCCCCCCCCCN1": 1}, catalyst="water / H3PO4", conditions="~280-300 °C under pressure"),
        R("PA12-2", "step_growth_condensation", "Condensation", "melt_polycondensation", "melt",
          {"NCCCCCCCCCCCC(=O)O": 1}, byproducts={H2O: 1}, rank=2,
          conditions="melt polycondensation of 12-aminododecanoic acid"),
    ], aliases=["nylon 12"]),

    P("Poly(p-phenylene terephthalamide)", "PPTA", "aramid",
      "[*]Nc1ccc(NC(=O)c2ccc(C([*])=O)cc2)cc1", [
        R("PPTA-1", "step_growth_condensation", "Condensation", "solution_polycondensation", "solution",
          {"Nc1ccc(N)cc1": 1, "O=C(Cl)c1ccc(C(=O)Cl)cc1": 1}, byproducts={HCl: 2},
          catalyst="none; NMP/CaCl2 solvent", conditions="low temperature (0-20 °C) solution polycondensation"),
        R("PPTA-2", "step_growth_condensation", "Condensation", "solution_polycondensation", "solution",
          {"Nc1ccc(N)cc1": 1, "OC(=O)c1ccc(C(=O)O)cc1": 1}, byproducts={H2O: 2}, scale="lab", rank=3,
          catalyst="triphenyl phosphite / pyridine (Yamazaki-Higashi)", conditions="NMP/LiCl, ~100 °C"),
    ], aliases=["Kevlar", "Twaron"]),

    P("Polyimide (PMDA-ODA)", "PI", "polyimide",
      "[*]c1ccc(Oc2ccc(N3C(=O)c4cc5c(cc4C3=O)C(=O)N([*])C5=O)cc2)cc1", [
        R("PI-1", "step_growth_condensation", "Condensation", "solution_polycondensation", "solution",
          {"O=C1OC(=O)c2cc3c(cc21)C(=O)OC3=O": 1, "Nc1ccc(Oc2ccc(N)cc2)cc1": 1}, byproducts={H2O: 2},
          catalyst="none (step 2: thermal or Ac2O/pyridine imidisation)",
          conditions="step 1: poly(amic acid) in DMAc/NMP at 0-25 °C; step 2: film cast, cured stepwise to 300-350 °C"),
    ], aliases=["Kapton"]),

    # ------------------------------------------------------------------ polyurethanes
    P("Thermoplastic polyurethane hard segment (MDI-BDO)", "TPU", "polyurethane",
      "[*]OCCCCOC(=O)Nc1ccc(Cc2ccc(NC([*])=O)cc2)cc1", [
        R("TPU-1", "step_growth_polyaddition", "Addition", "polyaddition", "melt",
          {"O=C=Nc1ccc(Cc2ccc(N=C=O)cc2)cc1": 1, "OCCCCO": 1},
          catalyst="dibutyltin dilaurate or bismuth carboxylate (optional)",
          conditions="reactive extrusion 180-240 °C (one-shot) or prepolymer route with a polyol soft segment"),
    ]),

    # ------------------------------------------------------------------ high-performance aromatics
    P("Poly(ether ether ketone)", "PEEK", "paek", "[*]Oc1ccc(Oc2ccc(C(=O)c3ccc([*])cc3)cc2)cc1", [
        R("PEEK-1", "step_growth_condensation", "Condensation", "snar_polycondensation", "solution",
          {"Oc1ccc(O)cc1": 1, "O=C(c1ccc(F)cc1)c1ccc(F)cc1": 1}, byproducts={HF: 2},
          catalyst="Na2CO3/K2CO3 (HF formally removed as KF + CO2 + H2O)",
          conditions="diphenyl sulfone, staged heating to 300-320 °C"),
    ]),

    P("Poly(ether ketone)", "PEK", "paek", "[*]Oc1ccc(C(=O)c2ccc([*])cc2)cc1", [
        R("PEK-1", "step_growth_condensation", "Condensation", "snar_polycondensation", "solution",
          {"O=C(c1ccc(O)cc1)c1ccc(F)cc1": 1}, byproducts={HF: 1}, catalyst="K2CO3",
          conditions="diphenyl sulfone, ~300-340 °C"),
        R("PEK-2", "step_growth_condensation", "Condensation", "friedel_crafts", "solution",
          {"O=C(Cl)c1ccc(Oc2ccccc2)cc1": 1}, byproducts={HCl: 1}, rank=2,
          catalyst="AlCl3 (with Lewis base) or HF/BF3", conditions="electrophilic acylation, -10 to 30 °C"),
    ]),

    P("Polysulfone (bisphenol A)", "PSU", "polysulfone",
      "[*]Oc1ccc(C(C)(C)c2ccc(Oc3ccc(S(=O)(=O)c4ccc([*])cc4)cc3)cc2)cc1", [
        R("PSU-1", "step_growth_condensation", "Condensation", "snar_polycondensation", "solution",
          {"CC(C)(c1ccc(O)cc1)c1ccc(O)cc1": 1, "O=S(=O)(c1ccc(Cl)cc1)c1ccc(Cl)cc1": 1},
          byproducts={HCl: 2}, catalyst="K2CO3 (or NaOH to form the bisphenolate)",
          conditions="DMAc or NMP with toluene azeotrope, 150-190 °C"),
    ], aliases=["Udel"]),

    P("Poly(p-phenylene sulfide)", "PPS", "polysulfide", "[*]Sc1ccc([*])cc1", [
        R("PPS-1", "step_growth_condensation", "Condensation", "solution_polycondensation", "solution",
          {"Clc1ccc(Cl)cc1": 1, "[Na]S[Na]": 1}, byproducts={"[Na]Cl": 2},
          catalyst="none (Phillips process)", conditions="NMP, 220-280 °C under pressure"),
    ], aliases=["Ryton"]),

    P("Poly(2,6-dimethyl-1,4-phenylene oxide)", "PPO", "polyether", "[*]Oc1c(C)cc([*])cc1C", [
        R("PPO-1", "oxidative_coupling", "Condensation", "oxidative_coupling", "solution",
          {"Cc1cccc(C)c1O": 1}, coreactants={O2: 0.5}, byproducts={H2O: 1},
          catalyst="CuBr / amine (e.g. dibutylamine)", conditions="toluene, O2, 25-50 °C"),
        R("PPO-2", "oxidative_coupling", "Condensation", "enzymatic", "solution",
          {"Cc1cccc(C)c1O": 1}, coreactants={H2O2: 1}, byproducts={H2O: 2},
          catalyst="horseradish peroxidase / H2O2", scale="lab", rank=3,
          conditions="aqueous organic mixture, room temperature"),
    ], aliases=["PPE", "polyphenylene ether", "Noryl base resin"]),

    P("Polyaniline (simplified leucoemeraldine repeat unit)", "PANI", "conducting", "[*]Nc1ccc([*])cc1", [
        R("PANI-1", "oxidative_coupling", "Condensation", "oxidative_coupling", "solution",
          {"Nc1ccccc1": 1}, coreactants={O2: 0.5}, byproducts={H2O: 1},
          catalyst="ammonium persulfate in 1 M HCl (O2 written as the formal oxidant)",
          conditions="0-5 °C aqueous acid"),
        R("PANI-2", "oxidative_coupling", "Condensation", "enzymatic", "solution",
          {"Nc1ccccc1": 1}, coreactants={H2O2: 1}, byproducts={H2O: 2},
          catalyst="horseradish peroxidase / H2O2 with polyanionic template (e.g. sulfonated polystyrene)",
          scale="lab", rank=3, conditions="aqueous buffer pH ~4, room temperature",
          notes="Template suppresses branched/ortho coupling"),
    ]),

    # ------------------------------------------------------------------ silicones
    P("Polydimethylsiloxane", "PDMS", "silicone", "[*]O[Si](C)(C)[*]", [
        R("PDMS-1", "step_growth_condensation", "Condensation", "hydrolytic_polycondensation", "solution",
          {"C[Si](C)(Cl)Cl": 1}, coreactants={H2O: 1}, byproducts={HCl: 2},
          conditions="hydrolysis of dimethyldichlorosilane, then acid/base condensation"),
        R("PDMS-2", "ring_opening", "Ringopening", "anionic", "bulk",
          {"C[Si]1(C)O[Si](C)(C)O[Si](C)(C)O[Si](C)(C)O1": 1}, ru_per_set=4,
          catalyst="KOH or tetramethylammonium silanolate (equilibration)", rank=1,
          conditions="120-160 °C; end-blocker (MM) sets chain length"),
        R("PDMS-3", "ring_opening", "Ringopening", "anionic", "solution",
          {"C[Si]1(C)O[Si](C)(C)O[Si](C)(C)O1": 1}, ru_per_set=3, catalyst="n-BuLi / THF promoter",
          scale="lab", rank=3, conditions="living kinetically controlled ROP of D3, 0-25 °C"),
    ], aliases=["silicone", "dimethicone"]),

    # ------------------------------------------------------------------ elastomers
    P("cis-1,4-Polyisoprene", "IR", "diene_elastomer", "[*]CC(C)=CC[*]", [
        R("IR-1", "chain_growth", "Addition", "ziegler_natta", "solution",
          {"C=CC(C)=C": 1}, catalyst="Nd carboxylate / Al alkyl / halide or TiCl4/AlR3",
          conditions="hexane, 30-60 °C; >97 % cis-1,4"),
        R("IR-2", "chain_growth", "Addition", "anionic", "solution",
          {"C=CC(C)=C": 1}, catalyst="n-BuLi or sec-BuLi", rank=2,
          conditions="hydrocarbon solvent, 40-70 °C; ~90-93 % cis-1,4"),
        R("IR-3", "biological", "Addition", "biological", "in_vivo",
          {"C=C(C)CCOP(=O)(O)OP(=O)(O)O": 1}, byproducts={"O=P(O)(O)OP(=O)(O)O": 1},
          catalyst="cis-prenyltransferase (rubber transferase) in Hevea brasiliensis latex",
          conditions="in vivo; natural rubber is harvested as latex",
          notes="Monomer is isopentenyl diphosphate (IPP); natural rubber, ~100 % cis-1,4"),
    ], aliases=["natural rubber", "NR", "polyisoprene"]),

    P("Polybutadiene", "BR", "diene_elastomer", "[*]CC=CC[*]", [
        R("BR-1", "chain_growth", "Addition", "ziegler_natta", "solution",
          {"C=CC=C": 1}, catalyst="Nd-, Co- or Ni-based Ziegler-Natta",
          conditions="hexane/cyclohexane, 50-80 °C; high-cis (>96 %)"),
        R("BR-2", "chain_growth", "Addition", "anionic", "solution",
          {"C=CC=C": 1}, catalyst="n-BuLi", rank=2, conditions="hydrocarbon, 50-100 °C; ~35-40 % cis"),
        R("BR-3", "chain_growth", "Addition", "free_radical", "emulsion",
          {"C=CC=C": 1}, catalyst="redox (hydroperoxide/Fe2+) cold emulsion", rank=3,
          conditions="5-10 °C; mostly used for SBR copolymers"),
    ]),

    P("Polychloroprene", "CR", "diene_elastomer", "[*]CC(Cl)=CC[*]", [
        R("CR-1", "chain_growth", "Addition", "free_radical", "emulsion",
          {"C=CC(Cl)=C": 1}, catalyst="persulfate, rosin soap, thiol/sulfur modifier",
          conditions="10-45 °C aqueous emulsion"),
    ], aliases=["Neoprene"]),

    # ------------------------------------------------------------------ ROMP polymers
    P("Polynorbornene", "PNB", "romp", "[*]C=CC1CCC([*])C1", [
        R("PNB-1", "romp", "ROMP", "romp", "solution",
          {"C1=CC2CCC1C2": 1}, catalyst="RuCl3 / HCl in butanol (Norsorex)",
          conditions="industrial ROMP; Grubbs/Schrock catalysts in the lab"),
    ], aliases=["Norsorex"]),

    P("Polyoctenamer", "TOR", "romp", "[*]C=CCCCCCC[*]", [
        R("TOR-1", "romp", "ROMP", "romp", "solution",
          {"C1=CCCCCCC1": 1}, catalyst="WCl6 / EtAlCl2", conditions="hexane, 0-40 °C"),
    ], aliases=["Vestenamer"]),
]
