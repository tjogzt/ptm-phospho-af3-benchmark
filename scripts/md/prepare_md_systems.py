#!/usr/bin/env python3
"""
prepare_md_systems.py — GROMACS MD system preparation
======================================================
Prepares MD simulation input files for phosphorylation benchmark proteins.
Generates system-specific GROMACS input files and job submission scripts.

Run LOCALLY on MacBook (data prep), then transfer to AutoDL for execution.

Workflow for each protein:
  1. Select protein (from benchmark dataset)
  2. Generate topology with phosphorylated residues (AMBER ff19SB + phos params)
  3. Solvate + add ions
  4. Generate MDP files (minimization, equilibration, production)
  5. Generate AutoDL job submission script

Usage:
  python prepare_md_systems.py --protein CDK2 [--all]
"""

import sys
import json
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "data"))
from benchmark_dataset import BENCHMARK_PROTEINS

# ================================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results" / "md"

# Target proteins for MD (priority subset: diverse effect types, manageable size)
MD_TARGETS = [
    "CDK2",       # local rearrangement, kinase, paired structures
    "ERK2",       # dual phosphorylation, local rearrangement
    "PKA_Calpha", # allosteric effect, paired structures
    "CheY",       # allosteric (uses BeF3 mimic, note special handling)
    "Aurora_A",   # local rearrangement
    "IRF3",       # allosteric + dimerization
    "Ubiquitin",  # negative control, minimal effect
]

# ================================================================

def generate_mdp_files(output_dir: Path, protein_name: str, nsteps_prod: int = 250_000_000):
    """Generate GROMACS MDP parameter files."""
    
    # --- Energy Minimization ---
    em_mdp = f"""; Energy Minimization
; Generated for {protein_name}
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
nstlist         = 10
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
vdwtype         = Cut-off
rvdw            = 1.0
pbc             = xyz
"""
    with open(output_dir / "em.mdp", 'w') as f:
        f.write(em_mdp)
    
    # --- NVT Equilibration ---
    nvt_mdp = f"""; NVT Equilibration (100 ps)
; Generated for {protein_name}
integrator      = md
dt              = 0.001
nsteps          = 100000
nstxout-compressed = 5000
nstlog          = 5000
nstenergy       = 5000
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
vdwtype         = Cut-off
rvdw            = 1.0
DispCorr        = EnerPres
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau-t           = 0.1 0.1
ref-t           = 310 310
pcoupl          = no
pbc             = xyz
constraints     = h-bonds
constraint-algorithm = LINCS
gen-vel         = yes
gen-temp        = 310
gen-seed        = -1
"""
    with open(output_dir / "nvt.mdp", 'w') as f:
        f.write(nvt_mdp)
    
    # --- NPT Equilibration ---
    npt_mdp = f"""; NPT Equilibration (200 ps)
; Generated for {protein_name}
integrator      = md
dt              = 0.002
nsteps          = 100000
nstxout-compressed = 5000
nstlog          = 5000
nstenergy       = 5000
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
vdwtype         = Cut-off
rvdw            = 1.0
DispCorr        = EnerPres
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau-t           = 0.1 0.1
ref-t           = 310 310
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau-p           = 2.0
ref-p           = 1.0
compressibility = 4.5e-5
pbc             = xyz
constraints     = h-bonds
constraint-algorithm = LINCS
"""
    with open(output_dir / "npt.mdp", 'w') as f:
        f.write(npt_mdp)
    
    # --- Production MD ---
    prod_mdp = f"""; Production MD (500 ns)
; Generated for {protein_name}
integrator      = md
dt              = 0.002
nsteps          = {nsteps_prod}
nstxout-compressed = 5000
nstlog          = 50000
nstenergy       = 50000
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
vdwtype         = Cut-off
rvdw            = 1.0
DispCorr        = EnerPres
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau-t           = 0.1 0.1
ref-t           = 310 310
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau-p           = 2.0
ref-p           = 1.0
compressibility = 4.5e-5
pbc             = xyz
constraints     = h-bonds
constraint-algorithm = LINCS
"""
    with open(output_dir / "prod.mdp", 'w') as f:
        f.write(prod_mdp)


def generate_gromacs_run_script(output_dir: Path, protein_name: str, state: str, 
                                 pdb_input: str, chain: str = "A"):
    """
    Generate a GROMACS run script for a specific protein/state.
    
    This script handles:
      1. pdb2gmx (topology generation with phosphorylated residue params)
      2. Solvation + ion addition
      3. Energy minimization
      4. NVT + NPT equilibration
      5. Production MD (500 ns)
    """
    
    script = f"""#!/bin/bash
# ================================================================
# GROMACS MD Simulation: {protein_name} ({state})
# ================================================================
# Auto-generated by prepare_md_systems.py
# Run on AutoDL with RTX 4090

set -e

PROTEIN="{protein_name}"
STATE="{state}"
PDB_INPUT="{pdb_input}"
CHAIN="{chain}"
NPUS=8  # CPU threads

echo "=== MD Simulation: $PROTEIN ($STATE) ==="

# Step 1: Prepare topology
echo "[1/7] Generating topology..."
# Note: phosphorylated residues need custom residue types.
# For AMBER ff19SB, phosphorylated residues are available as SEP, TPO, PTR.
# If using CHARMM36m, use the CHARMM-GUI generated topology instead.

gmx pdb2gmx -f ../data/pdb/$PDB_INPUT.pdb \\
    -o ${PROTEIN}_${STATE}_processed.gro \\
    -p ${PROTEIN}_${STATE}.top \\
    -ff amber99sb-ildn \\
    -water tip3p \\
    -ignh \\
    -chainsep interactive <<EOF
$CHAIN
EOF

# Step 2: Define simulation box
echo "[2/7] Defining simulation box..."
gmx editconf -f ${PROTEIN}_${STATE}_processed.gro \\
    -o ${PROTEIN}_${STATE}_box.gro \\
    -c -d 1.2 -bt cubic

# Step 3: Solvate
echo "[3/7] Solvating..."
gmx solvate -cp ${PROTEIN}_${STATE}_box.gro \\
    -cs spc216.gro \\
    -o ${PROTEIN}_${STATE}_solv.gro \\
    -p ${PROTEIN}_${STATE}.top

# Step 4: Add ions (0.15 M NaCl)
echo "[4/7] Adding ions..."
gmx grompp -f em.mdp \\
    -c ${PROTEIN}_${STATE}_solv.gro \\
    -p ${PROTEIN}_${STATE}.top \\
    -o ${PROTEIN}_${STATE}_ions.tpr \\
    -maxwarn 5

echo "SOL" | gmx genion -s ${PROTEIN}_${STATE}_ions.tpr \\
    -o ${PROTEIN}_${STATE}_neutral.gro \\
    -p ${PROTEIN}_${STATE}.top \\
    -pname NA -nname CL \\
    -neutral -conc 0.15

# Step 5: Energy minimization
echo "[5/7] Energy minimization..."
gmx grompp -f em.mdp \\
    -c ${PROTEIN}_${STATE}_neutral.gro \\
    -p ${PROTEIN}_${STATE}.top \\
    -o ${PROTEIN}_${STATE}_em.tpr \\
    -maxwarn 5

gmx mdrun -deffnm ${PROTEIN}_${STATE}_em \\
    -v -nt $NPUS -nb gpu

# Step 6: Equilibration (NVT + NPT)
echo "[6/7] Equilibration..."

# NVT
gmx grompp -f nvt.mdp \\
    -c ${PROTEIN}_${STATE}_em.gro \\
    -p ${PROTEIN}_${STATE}.top \\
    -o ${PROTEIN}_${STATE}_nvt.tpr \\
    -maxwarn 5

gmx mdrun -deffnm ${PROTEIN}_${STATE}_nvt \\
    -v -nt $NPUS -nb gpu

# NPT
gmx grompp -f npt.mdp \\
    -c ${PROTEIN}_${STATE}_nvt.gro \\
    -p ${PROTEIN}_${STATE}.top \\
    -o ${PROTEIN}_${STATE}_npt.tpr \\
    -maxwarn 5

gmx mdrun -deffnm ${PROTEIN}_${STATE}_npt \\
    -v -nt $NPUS -nb gpu

# Step 7: Production MD (500 ns)
echo "[7/7] Production MD (500 ns)..."
gmx grompp -f prod.mdp \\
    -c ${PROTEIN}_${STATE}_npt.gro \\
    -p ${PROTEIN}_${STATE}.top \\
    -o ${PROTEIN}_${STATE}_prod.tpr \\
    -maxwarn 5

gmx mdrun -deffnm ${PROTEIN}_${STATE}_prod \\
    -v -nt $NPUS -nb gpu -pin on

echo "=== MD Complete: $PROTEIN ($STATE) ==="

# Post-processing: strip water, center protein
echo "SOL" | gmx trjconv -f ${PROTEIN}_${STATE}_prod.xtc \\
    -s ${PROTEIN}_${STATE}_prod.tpr \\
    -o ${PROTEIN}_${STATE}_prod_processed.xtc \\
    -pbc mol -center -ur compact

echo "Output files:"
echo "  Trajectory: ${PROTEIN}_${STATE}_prod_processed.xtc"
echo "  Structure:  ${PROTEIN}_${STATE}_prod.gro"
echo "  Energy:     ${PROTEIN}_${STATE}_prod.edr"
"""
    
    script_path = output_dir / f"run_{protein_name}_{state}.sh"
    with open(script_path, 'w') as f:
        f.write(script)
    script_path.chmod(0o755)
    print(f"  [OK] MD script: {script_path}")


def main():
    parser = argparse.ArgumentParser(description="Prepare MD simulation systems")
    parser.add_argument("--protein", type=str, default="all",
                       help="Target protein name or 'all'")
    parser.add_argument("--output-dir", type=str, default="results/md",
                       help="Output directory")
    args = parser.parse_args()
    
    # Determine target proteins
    if args.protein == "all":
        targets = MD_TARGETS
    elif args.protein in BENCHMARK_PROTEINS:
        targets = [args.protein]
    else:
        print(f"Unknown protein: {args.protein}")
        print(f"Available: {list(BENCHMARK_PROTEINS.keys())}")
        sys.exit(1)
    
    output_root = Path(args.output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    
    # Save metadata
    md_metadata = {}
    
    for protein_name in targets:
        info = BENCHMARK_PROTEINS[protein_name]
        print(f"\n{'='*60}")
        print(f"Preparing MD system: {protein_name} ({info['gene']})")
        print(f"  Sites: {info['phos_sites']} | Type: {info['effect_type']}")
        
        # Create output dir
        protein_dir = output_root / protein_name
        protein_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate MDP files
        generate_mdp_files(protein_dir, protein_name)
        
        # Generate run scripts for each state
        pdbs = info["pdbs"]
        chain = info.get("chain", "A")
        
        # Unphos state
        if pdbs.get("unphos"):
            for pdb_id in pdbs["unphos"]:
                unphos_dir = protein_dir / "unphos"
                unphos_dir.mkdir(exist_ok=True)
                for f in protein_dir.glob("*.mdp"):
                    import shutil
                    shutil.copy(f, unphos_dir / f.name)
                generate_gromacs_run_script(unphos_dir, protein_name, "unphos", pdb_id, chain)
        
        # Phos state
        if pdbs.get("phos"):
            for pdb_id in pdbs["phos"]:
                phos_dir = protein_dir / "phos"
                phos_dir.mkdir(exist_ok=True)
                for f in protein_dir.glob("*.mdp"):
                    import shutil
                    shutil.copy(f, phos_dir / f.name)
                generate_gromacs_run_script(phos_dir, protein_name, "phos", pdb_id, chain)
        
        md_metadata[protein_name] = {
            "gene": info["gene"],
            "phos_sites": info["phos_sites"],
            "effect_type": info.get("effect_type"),
            "length": info.get("length"),
            "blind_test": info.get("blind_test", False),
        }
    
    # Save MD metadata
    import json
    meta_file = output_root / "md_metadata.json"
    with open(meta_file, 'w') as f:
        json.dump(md_metadata, f, indent=2)
    
    print(f"\n{'='*60}")
    print(f"MD preparation complete!")
    print(f"Target proteins: {len(targets)}")
    print(f"Output directory: {output_root}")
    print(f"Metadata: {meta_file}")
    print(f"\nNext step: Transfer to AutoDL and execute run scripts.")

if __name__ == "__main__":
    import shutil
    main()
