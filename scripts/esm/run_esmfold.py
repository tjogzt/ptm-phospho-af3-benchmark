#!/usr/bin/env python3
"""
run_esmfold.py — ESMFold batch inference for phosphorylation benchmark
======================================================================
Runs ESMFold v1 on all benchmark protein sequences.
ESMFold uses single-sequence input (no MSA) → ~60x faster than AlphaFold2.

Requirements (AutoDL A100 40G):
  pip install fair-esm torch biopython

Memory: ~35-40 GB VRAM for esmfold_v1 (3B model)
Runtime: ~2-3 min per protein (on A100)

Usage:
  python run_esmfold.py --fasta-dir data/fasta --output-dir results/esm
"""

import os
import sys
import json
import argparse
import time
from pathlib import Path
from typing import List, Dict

import torch
import esm
from esm import pretrained

# ================================================================

def parse_fasta(fasta_path: str) -> Dict[str, str]:
    """Parse FASTA files, return {header: sequence}."""
    sequences = {}
    with open(fasta_path) as f:
        current_header = None
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences[current_header] = ''.join(current_seq)
    return sequences

def run_esmfold(sequence: str, model, device: str = "cuda:0"):
    """Run ESMFold on a single sequence."""
    with torch.no_grad():
        output = model.infer([("protein", sequence)], num_recycles=4)
    return output

def extract_metrics(output, protein_name: str) -> dict:
    """Extract confidence metrics from ESMFold output."""
    plddt = output["plddt"][0].cpu().numpy()  # per-residue pLDDT
    ptm = float(output["ptm"][0].cpu().numpy())  # predicted TM-score
    
    metrics = {
        "name": protein_name,
        "ptm": ptm,
        "mean_plddt": float(plddt.mean()),
        "plddt_per_residue": plddt.tolist(),
    }
    return metrics

def main():
    parser = argparse.ArgumentParser(description="ESMFold batch inference")
    parser.add_argument("--fasta-dir", type=str, default="data/fasta",
                       help="Directory with FASTA files")
    parser.add_argument("--output-dir", type=str, default="results/esm",
                       help="Output directory for predictions")
    parser.add_argument("--gpu", type=int, default=0, help="GPU device ID")
    parser.add_argument("--chunk-size", type=int, default=128,
                       help="Chunk size for long sequences (reduce if OOM)")
    args = parser.parse_args()
    
    fasta_dir = Path(args.fasta_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pdb_dir = output_dir / "pdb"
    pdb_dir.mkdir(exist_ok=True)
    
    device = f"cuda:{args.gpu}" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")
    
    # Load ESMFold model
    print("Loading ESMFold v1 model...")
    model = esm.pretrained.esmfold_v1()
    model = model.eval().to(device)
    model.set_chunk_size(args.chunk_size)
    print(f"Model loaded. VRAM used: {torch.cuda.memory_allocated()/1e9:.1f} GB")
    
    # Find all unique proteins (deduplicate by unphos vs phos states)
    all_metrics = []
    
    for fasta_file in sorted(fasta_dir.glob("*.fasta")):
        print(f"\nProcessing: {fasta_file.name}")
        sequences = parse_fasta(str(fasta_file))
        
        for header, seq in sequences.items():
            # Parse header: gene|uniprot|chain_X or gene_state_pdb_uniprot
            parts = header.split('|')
            if len(parts) >= 3:
                protein_name = parts[0]
            else:
                protein_name = header.split('_')[0]
            
            print(f"  {protein_name}: {len(seq)} residues")
            
            try:
                output = run_esmfold(seq, model, device)
                metrics = extract_metrics(output, protein_name)
                all_metrics.append(metrics)
                
                # Save PDB file
                pdb_content = esm.pretrained.esmfold_v1().output_to_pdb(output)[0]
                pdb_file = pdb_dir / f"{protein_name}_esmfold.pdb"
                with open(pdb_file, 'w') as f:
                    f.write(pdb_content)
                
                print(f"    pTM={metrics['ptm']:.3f}, mean_pLDDT={metrics['mean_plddt']:.1f}")
                
            except torch.cuda.OutOfMemoryError:
                print(f"    [OOM] Skipping {protein_name} (too large for GPU memory)")
                # Try with larger chunk size next time
                torch.cuda.empty_cache()
            except Exception as e:
                print(f"    [ERROR] {protein_name}: {e}")
    
    # Save metrics summary
    metrics_file = output_dir / "esmfold_metrics.json"
    with open(metrics_file, 'w') as f:
        json.dump(all_metrics, f, indent=2)
    
    print(f"\n=== Complete ===")
    print(f"Processed {len(all_metrics)} sequences")
    print(f"Metrics saved to: {metrics_file}")
    print(f"PDB files saved to: {pdb_dir}")

if __name__ == "__main__":
    main()
