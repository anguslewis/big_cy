#!/bin/bash
# =============================================================================
# setup_klrep_env_local.sh — create the conda env for the KL-model replication
# (big_cy task Q0). Designed to run LOCALLY on a Mac (Apple Silicon / arm64).
#
# The replication is tensor-native (PyTorch), device-agnostic, GPU-ready. Note:
# Apple MPS is float32-only, so the float64 model solve runs on CPU locally; the
# same code targets CUDA/H100 later (float64-capable) via the Sherlock env.
#
# Usage:  bash klreplication/setup_klrep_env_local.sh
# Then:   conda activate big_cy_klrep_env
# =============================================================================
set -e  # Exit on error

ENV_NAME="big_cy_klrep_env"
PYTHON_VERSION="3.11"

echo "=========================================="
echo "Setting up conda environment: $ENV_NAME (local)"
echo "=========================================="

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda is not available. Please install conda/miniconda first."
    echo "You can install miniconda from: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Initialize conda for bash (if not already initialized)
if ! conda info --envs &> /dev/null; then
    echo "Initializing conda..."
    eval "$(conda shell.bash hook)"
fi

# Remove existing environment if it exists
if conda env list | grep -q "^${ENV_NAME}\s"; then
    echo "Removing existing environment: $ENV_NAME..."
    conda env remove -n $ENV_NAME -y
fi

# Create new conda environment with Python
echo "Creating conda environment with Python $PYTHON_VERSION..."
conda create -n "$ENV_NAME" python=$PYTHON_VERSION -y

# Activate the environment
echo "Activating environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$ENV_NAME"

# Install everything via conda (conda first; pip only for genuine conda gaps).
# PyTorch comes from the pytorch channel — on Apple Silicon this is the CPU + MPS
# build. The float64 solve uses CPU (MPS lacks float64); the code is device-
# agnostic so the CUDA build (Sherlock script) runs the same source.
echo "Installing core packages + PyTorch via conda..."
conda install -y -c pytorch -c conda-forge \
    numpy \
    scipy \
    pandas \
    matplotlib \
    pytest \
    pytorch

# Verify installation
echo ""
echo "=========================================="
echo "Verifying installation..."
echo "=========================================="
python -c "import numpy, scipy, pandas, matplotlib, pytest; print('Core packages OK')"
python -c "import torch; print('torch', torch.__version__, '| mps', torch.backends.mps.is_available(), '| cuda', torch.cuda.is_available())"
python -c "import torch; x=torch.ones(3, dtype=torch.float64); print('float64 on CPU OK:', (x.sum().item()==3.0))"

echo ""
echo "=========================================="
echo "Environment setup complete!"
echo "=========================================="
echo ""
echo "To activate the environment, run:"
echo "  conda activate $ENV_NAME"
echo ""
echo "To run the replication tests:"
echo "  conda activate $ENV_NAME && pytest klreplication/tests -q"
echo ""
