#!/bin/bash
# ================================================================
# setup_autodl_esm.sh — AutoDL ESMFold environment setup
# ================================================================
# Run this on an AutoDL instance with A100 40G
# Recommended image: PyTorch 2.3.0 + CUDA 12.1 + Python 3.10
#
# Usage: bash setup_autodl_esm.sh
# ================================================================

set -e

echo "=== ESMFold AutoDL Environment Setup ==="

# 1. Check CUDA
echo "[1/5] Checking CUDA..."
nvidia-smi
python -c "import torch; print(f'PyTorch {torch.__version__}, CUDA available: {torch.cuda.is_available()}')"

# 2. Install dependencies
echo "[2/5] Installing dependencies..."
pip install --upgrade pip
pip install fair-esm  # ESM-2 + ESMFold
pip install biopython pandas tqdm
pip install mdanalysis mdtraj  # for structure analysis

# 3. Verify ESMFold installation
echo "[3/5] Verifying ESMFold..."
python -c "
import esm
print('ESM version:', esm.__version__ if hasattr(esm, '__version__') else 'installed')
# Test loading a small model
model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
print('ESM-2 8M model loaded successfully')
"

# 4. Create directory structure
echo "[4/5] Creating directories..."
mkdir -p /root/ptm_project/{data/{fasta,pdb},results/esm,scripts}

# 5. Download benchmark data from local machine
echo "[5/5] Ready for data transfer"
echo ""
echo "=== Setup Complete ==="
echo ""
echo "Next steps:"
echo "  1. Transfer FASTA files from MacBook:"
echo "     scp data/fasta/*.fasta root@<autoDL-IP>:/root/ptm_project/data/fasta/"
echo ""
echo "  2. Transfer and run ESMFold inference:"
echo "     scp scripts/esm/run_esmfold.py root@<autoDL-IP>:/root/ptm_project/scripts/"
echo "     ssh root@<autoDL-IP> 'cd /root/ptm_project && python scripts/run_esmfold.py'"
echo ""
echo "  3. Download results back:"
echo "     scp root@<autoDL-IP>:/root/ptm_project/results/esm/* results/esm/"
