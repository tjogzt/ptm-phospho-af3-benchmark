#!/usr/bin/env python3
"""
download_pdbs.py — Batch PDB structure download
=================================================
Downloads all PDB structures from the benchmark dataset.
Also extracts FASTA sequences for AF3/ESMFold input.

Usage: python download_pdbs.py [--limit N]
"""

import os
import sys
import time
import json
import argparse
from pathlib import Path
from urllib.request import urlretrieve
from io import StringIO

import requests
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.SeqUtils import seq1

# Add parent dir to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "data"))
from benchmark_dataset import BENCHMARK_PROTEINS

# === Config ===
PROJECT_ROOT = Path(__file__).parent.parent
PDB_DIR = PROJECT_ROOT / "data" / "pdb"
FASTA_DIR = PROJECT_ROOT / "data" / "fasta"
METADATA_DIR = PROJECT_ROOT / "data" / "metadata"

PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
PDB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"

# Rate limiting: be nice to RCSB servers
DOWNLOAD_DELAY = 0.5  # seconds between downloads

# ================================================================

def ensure_dirs():
    """Create output directories."""
    for d in [PDB_DIR, FASTA_DIR, METADATA_DIR]:
        d.mkdir(parents=True, exist_ok=True)

def download_pdb(pdb_id, file_format="pdb"):
    """Download a PDB structure file from RCSB."""
    if file_format == "pdb":
        url = PDB_DOWNLOAD_URL.format(pdb_id=pdb_id.upper())
        ext = ".pdb"
    else:
        url = PDB_CIF_URL.format(pdb_id=pdb_id.upper())
        ext = ".cif"
    
    dest = PDB_DIR / f"{pdb_id.upper()}{ext}"
    
    if dest.exists():
        print(f"  [SKIP] {pdb_id} already exists")
        return True
    
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            with open(dest, 'w') as f:
                f.write(resp.text)
            print(f"  [OK] {pdb_id} downloaded ({len(resp.text)} bytes)")
            time.sleep(DOWNLOAD_DELAY)
            return True
        else:
            print(f"  [FAIL] {pdb_id}: HTTP {resp.status_code}")
            return False
    except Exception as e:
        print(f"  [ERROR] {pdb_id}: {e}")
        return False

def extract_sequence_from_pdb(pdb_id):
    """Extract amino acid sequence from PDB file using BioPython."""
    pdb_file = PDB_DIR / f"{pdb_id.upper()}.pdb"
    if not pdb_file.exists():
        print(f"  [WARN] PDB file not found: {pdb_file}")
        return None
    
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, str(pdb_file))
        model = structure[0]
        
        sequences = {}
        for chain in model:
            residues = []
            for res in chain:
                if res.get_id()[0] == ' ':  # Standard residues only
                    try:
                        residues.append(seq1(res.get_resname()))
                    except:
                        residues.append('X')
            if residues:
                sequences[chain.get_id()] = ''.join(residues)
        
        return sequences
    except Exception as e:
        print(f"  [ERROR] parsing {pdb_id}: {e}")
        return None

def save_fasta(uniprot_id, gene_name, sequences):
    """Save extracted sequences as FASTA files."""
    fasta_file = FASTA_DIR / f"{gene_name}_{uniprot_id}.fasta"
    with open(fasta_file, 'w') as f:
        for chain_id, seq in sequences.items():
            f.write(f">{gene_name}|{uniprot_id}|chain_{chain_id}\n")
            f.write(f"{seq}\n")
    print(f"  [FASTA] Written to {fasta_file}")

def build_metadata():
    """Compile benchmark metadata JSON for downstream analysis."""
    metadata = {}
    for name, info in BENCHMARK_PROTEINS.items():
        all_pdbs = []
        for state, ids in info.get("pdbs", {}).items():
            all_pdbs.extend(ids)
        
        metadata[name] = {
            "gene": info["gene"],
            "uniprot": info["uniprot"],
            "description": info["description"],
            "phos_sites": info["phos_sites"],
            "phos_type": info.get("phos_type", "SEP"),
            "effect_type": info.get("effect_type", "unknown"),
            "all_pdbs": all_pdbs,
            "chain": info.get("chain", "A"),
            "length": info.get("length", 0),
            "blind_test": info.get("blind_test", False),
            "note": info.get("note", ""),
        }
    
    meta_file = METADATA_DIR / "benchmark_metadata.json"
    with open(meta_file, 'w') as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)
    print(f"\n[OK] Metadata written to {meta_file}")

def main():
    parser = argparse.ArgumentParser(description="Download PDB structures for PTM benchmark")
    parser.add_argument("--limit", type=int, default=0, help="Limit to first N proteins")
    parser.add_argument("--skip-download", action="store_true", help="Skip PDB download, only build metadata/FASTA")
    args = parser.parse_args()
    
    ensure_dirs()
    
    # Collect all unique PDB IDs
    all_pdb_ids = set()
    proteins = list(BENCHMARK_PROTEINS.items())
    if args.limit > 0:
        proteins = proteins[:args.limit]
    
    for name, info in proteins:
        for state, ids in info.get("pdbs", {}).items():
            all_pdb_ids.update(ids)
    
    print(f"=== PTM Benchmark PDB Download ===")
    print(f"Proteins: {len(proteins)} | Unique PDBs: {len(all_pdb_ids)}")
    print(f"Target dir: {PDB_DIR}\n")
    
    # Download PDBs
    if not args.skip_download:
        success = 0
        for pdb_id in sorted(all_pdb_ids):
            if download_pdb(pdb_id):
                success += 1
        
        print(f"\nDownload complete: {success}/{len(all_pdb_ids)} successful")
    
    # Extract FASTA sequences
    print(f"\n=== Extracting FASTA sequences ===")
    for name, info in proteins:
        print(f"\n[{name}] ({info['gene']})")
        for state, ids in info.get("pdbs", {}).items():
            for pdb_id in ids:
                seqs = extract_sequence_from_pdb(pdb_id)
                if seqs:
                    # Use the primary chain
                    chain = info.get("chain", "A")
                    if chain in seqs:
                        target_seq = {chain: seqs[chain]}
                        save_fasta(info["uniprot"], f"{info['gene']}_{state}_{pdb_id}", target_seq)
                    else:
                        # Chain not found, take first chain
                        first_chain = list(seqs.keys())[0]
                        print(f"  [INFO] Chain {chain} not found, using {first_chain}")
                        save_fasta(info["uniprot"], f"{info['gene']}_{state}_{pdb_id}", 
                                  {first_chain: seqs[first_chain]})
    
    # Build metadata
    build_metadata()
    print("\n=== Done ===")

if __name__ == "__main__":
    main()
