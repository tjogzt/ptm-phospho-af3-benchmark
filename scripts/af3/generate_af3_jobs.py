#!/usr/bin/env python3
"""
generate_af3_jobs.py — Generate AlphaFold3 batch JSON for submission
====================================================================
Creates AF3-compatible JSON files for batch upload to alphafoldserver.com.

For each benchmark protein, generates:
  - P-group: Phosphorylated structure prediction (5 seeds)
  - U-group: Unphosphorylated baseline (5 seeds)  
  - T-group: Phosphorylated with custom template (5 seeds, optional)
  - A-group: Phosphorylated WITHOUT template (ablation, for blind test)

AF3 JSON format (v3):
  - sequences[].proteinChain.sequence: standard AA codes
  - sequences[].proteinChain.modifications: [{ptmType, ptmPosition}]
  - sequences[].proteinChain.count: copy number
  - name: job name
  - modelSeeds: array of seed integers

Phosphorylation CCD codes:
  SEP = phosphoserine, TPO = phosphothreonine, PTR = phosphotyrosine, NEP = phosphohistidine

Usage:
  python generate_af3_jobs.py [--output af3_jobs.json] [--limit N]
"""

import sys
import json
import argparse
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent.parent  # scripts/af3/ -> project root
sys.path.insert(0, str(_ROOT / "data"))
from benchmark_dataset import BENCHMARK_PROTEINS

# ================================================================

PROJECT_ROOT = _ROOT
DATA_DIR = PROJECT_ROOT / "data"
FASTA_DIR = DATA_DIR / "fasta"

# AF3 PTM CCD codes
PTM_CCD_CODES = {
    "SEP": "SEP",   # Phosphoserine
    "TPO": "TPO",   # Phosphothreonine
    "PTR": "PTR",   # Phosphotyrosine
    "NEP": "NEP",   # Phosphohistidine
}

def get_sequence_from_fasta(gene, uniprot, state, pdb_id):
    """Read protein sequence from previously saved FASTA file."""
    fasta_file = FASTA_DIR / f"{gene}_{state}_{pdb_id}_{uniprot}.fasta"
    if not fasta_file.exists():
        print(f"  [WARN] FASTA not found: {fasta_file}")
        return None
    
    with open(fasta_file) as f:
        lines = f.readlines()
        if len(lines) >= 2:
            seq = ''.join(line.strip() for line in lines[1:])
            return seq
    return None

def build_modified_sequence(sequence, phos_sites, phos_type):
    """
    Build sequence string with modification annotations.
    
    AF3 format: standard sequence + modifications array.
    Each modification: {"ptmType": "SEP", "ptmPosition": 160}
    """
    modifications = []
    
    for i, site_str in enumerate(phos_sites):
        # Parse site: "Thr160" -> AA="T", pos=160
        aa = site_str[0]  # First letter of residue type
        pos = int(''.join(c for c in site_str if c.isdigit()))
        
        # Determine CCD code based on phos_type or residue
        if phos_type == "SEP":
            code = "SEP"
        elif phos_type == "TPO":
            code = "TPO"
        elif phos_type == "PTR":
            code = "PTR"
        elif phos_type == "TPO_PTR":
            # Dual phosphorylation: Thr183->TPO, Tyr185->PTR
            if aa in ['T']:
                code = "TPO"
            elif aa in ['Y']:
                code = "PTR"
            else:
                code = "SEP"
        else:
            # Infer from site name
            if aa == 'S':
                code = "SEP"
            elif aa == 'T':
                code = "TPO"
            elif aa == 'Y':
                code = "PTR"
            else:
                code = "SEP"
        
        modifications.append({
            "ptmType": code,
            "ptmPosition": pos
        })
    
    return sequence, modifications

def create_af3_job(name, sequence, modifications=None, seeds=None):
    """Create a single AF3 job JSON object."""
    if seeds is None:
        seeds = [1, 2, 3, 4, 5]
    
    protein_chain = {
        "sequence": sequence,
        "count": 1
    }
    
    if modifications:
        protein_chain["modifications"] = modifications
    
    job = {
        "name": name,
        "modelSeeds": seeds,
        "sequences": [
            {"proteinChain": protein_chain}
        ]
    }
    
    return job

def generate_jobs(limit=0):
    """Generate all AF3 jobs for the benchmark."""
    jobs = []
    stats = {"p_group": 0, "u_group": 0, "t_group": 0}
    
    proteins = list(BENCHMARK_PROTEINS.items())
    if limit > 0:
        proteins = proteins[:limit]
    
    for name, info in proteins:
        gene = info["gene"]
        uniprot = info["uniprot"]
        phos_sites = info["phos_sites"]
        phos_type = info["phos_type"]
        pdbs = info["pdbs"]
        chain = info.get("chain", "A")
        blind_test = info.get("blind_test", False)
        
        # Skip CheY (Asp phosphorylation not directly supported)
        if phos_type is None:
            print(f"  [SKIP] {name}: Asp phosphorylation, not supported by AF3")
            continue
        
        # Get unphosphorylated sequence
        unphos_seq = None
        if pdbs.get("unphos"):
            for pdb_id in pdbs["unphos"]:
                unphos_seq = get_sequence_from_fasta(gene, uniprot, "unphos", pdb_id)
                if unphos_seq:
                    break
        
        # Get phosphorylated sequence
        phos_seq = None
        if pdbs.get("phos"):
            for pdb_id in pdbs["phos"]:
                phos_seq = get_sequence_from_fasta(gene, uniprot, "phos", pdb_id)
                if phos_seq:
                    break
        
        # If no phos FASTA, use unphos sequence
        if not phos_seq and unphos_seq:
            phos_seq = unphos_seq
        if not unphos_seq and phos_seq:
            unphos_seq = phos_seq
        if not phos_seq:
            print(f"  [WARN] {name}: no sequence available")
            continue
        
        # Build modified sequence with PTM annotations
        mod_seq, modifications = build_modified_sequence(phos_seq, phos_sites, phos_type)
        
        # U-group: Unphosphorylated baseline
        job_name_u = f"{gene}_unphos"
        jobs.append(create_af3_job(job_name_u, unphos_seq))
        stats["u_group"] += 1
        
        # P-group: Phosphorylated
        job_name_p = f"{gene}_phos"
        jobs.append(create_af3_job(job_name_p, mod_seq, modifications))
        stats["p_group"] += 1
        
        # Additional seeds for key proteins (dual phos, allosteric)
        if phos_type == "TPO_PTR" or info.get("effect_type") == "allosteric":
            job_name_p10 = f"{gene}_phos_seeds6-10"
            jobs.append(create_af3_job(job_name_p10, mod_seq, modifications, 
                                       seeds=[6, 7, 8, 9, 10]))
            stats["p_group"] += 1
        
        print(f"  [OK] {gene}: U+P jobs generated (sites: {phos_sites})")
    
    print(f"\n=== Job Summary ===")
    print(f"Total jobs: {len(jobs)}")
    print(f"  P-group: {stats['p_group']}")
    print(f"  U-group: {stats['u_group']}")
    print(f"Estimated runtime: {len(jobs) * 5 / 30:.1f} days (at 30 jobs/day quota)")
    
    return jobs

def main():
    parser = argparse.ArgumentParser(description="Generate AF3 batch JSON jobs")
    parser.add_argument("--output", type=str, default="af3_batch_jobs.json",
                       help="Output JSON file path")
    parser.add_argument("--limit", type=int, default=0, help="Limit to N proteins")
    parser.add_argument("--split", type=int, default=0, 
                       help="Split into N files (for large batches)")
    args = parser.parse_args()
    
    print("=== AF3 Batch Job Generator ===\n")
    print("Analyzing benchmark dataset...")
    
    jobs = generate_jobs(limit=args.limit)
    
    if args.split > 0:
        # Split into chunks (AF3 accepts up to 100 jobs per upload)
        chunk_size = min(100, max(1, len(jobs) // args.split))
        for i in range(0, len(jobs), chunk_size):
            chunk = jobs[i:i+chunk_size]
            out_path = PROJECT_ROOT / f"af3_jobs_part{i//chunk_size + 1}.json"
            with open(out_path, 'w') as f:
                json.dump(chunk, f, indent=2)
            print(f"\n[OK] Written {len(chunk)} jobs to {out_path}")
    else:
        out_path = PROJECT_ROOT / args.output
        with open(out_path, 'w') as f:
            json.dump(jobs, f, indent=2)
        print(f"\n[OK] Written {len(jobs)} jobs to {out_path}")
    
    # Print submission reminder
    print(f"""
=== SUBMISSION INSTRUCTIONS ===
1. Go to https://alphafoldserver.com
2. Sign in with Google account
3. Click 'Upload JSON' 
4. Select the generated JSON file
5. Wait for all jobs to complete (check quota: 30 jobs/day)

Tips:
- AF3 Server quota: 30 jobs/day
- Each job produces 5 predictions (unless custom seeds)
- Download results as ZIP for each completed job
- Save ZIP files to: results/af3/
""")

if __name__ == "__main__":
    main()
