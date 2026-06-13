#!/bin/bash
# =============================================================================
# setup_klrep_env_sherlock.sh — create the conda env for the KL-model
# replication (big_cy task Q0) on the Sherlock HPC cluster.
#
# REFERENCE / FUTURE USE ONLY. Per the big_cy standing policy all work currently
# runs locally (no Sherlock). This script exists so the GPU (H100) move is a
# one-command setup when the time comes. It is NOT run as part of the local task.
#
# The replication is tensor-native (PyTorch) and float64. H100s are float64-
# capable, so the model solve runs on GPU there with the same source as local.
#
# Before running on a compute node WITH INTERNET + a CUDA module loaded:
#   - ensure $big_cy points at the project oak root (see project_strings.py)
#   - load the cluster's CUDA toolkit module so the torch CUDA wheel matches the
#     driver (adjust the --index-url cuXXX tag below to the cluster's CUDA).
#   bash setup_klrep_env_sherlock.sh
# =============================================================================
set -e  # Exit on error

ENV_NAME="big_cy_klrep_env"
PYTHON_VERSION="3.11"

# torch CUDA version for the conda pytorch-cuda meta-package — MATCH to the
# cluster's CUDA toolkit (e.g. 12.1, 12.4). H100 needs CUDA 12.x. Override by
# exporting PYTORCH_CUDA before running.
PYTORCH_CUDA="${PYTORCH_CUDA:-12.1}"

echo "=========================================="
echo "Setting up conda environment: $ENV_NAME (Sherlock/HPC)"
echo "=========================================="

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda is not available. Please load conda module first."
    echo "Example: module load conda   (or source the group miniconda profile)"
    exit 1
fi

# Project root on oak (data root); fall back sensibly. Env + pkgs cache live
# here / on scratch to save home-directory space (mirrors rier/lx1 convention).
BIG_CY_ROOT="${big_cy:-/oak/stanford/groups/maggiori/Lab/lewis/big_cy}"

USER_CONDA_PKGS="${SCRATCH:-$BIG_CY_ROOT/.conda}/conda_pkgs_klrep"
USER_CONDA_ENVS="$BIG_CY_ROOT/.conda/envs"
mkdir -p "$USER_CONDA_PKGS" "$USER_CONDA_ENVS"
export CONDA_PKGS_DIRS="$USER_CONDA_PKGS"
export CONDA_ENVS_DIRS="$USER_CONDA_ENVS"

echo "Using conda directories:"
echo "  Packages cache: $CONDA_PKGS_DIRS"
echo "  Environments:   $CONDA_ENVS_DIRS"
echo ""

ENV_PREFIX="$USER_CONDA_ENVS/$ENV_NAME"

# Remove existing environment if it exists
if [ -d "$ENV_PREFIX" ]; then
    echo "Environment at $ENV_PREFIX already exists. Removing it..."
    conda env remove -p "$ENV_PREFIX" -y 2>/dev/null || rm -rf "$ENV_PREFIX"
fi
if conda env list | grep -q "${ENV_NAME}"; then
    conda env remove -n "$ENV_NAME" -y 2>/dev/null || true
fi

# Create new conda environment with Python at the user-writable prefix
echo "Creating conda environment with Python $PYTHON_VERSION at $ENV_PREFIX..."
conda create --prefix "$ENV_PREFIX" python=$PYTHON_VERSION -y

# Register the envs directory so the env can be activated by name
conda config --append envs_dirs "$USER_CONDA_ENVS" 2>/dev/null || {
    echo "Note: Could not update conda config. Activate using the full path."
}

# Activate the environment
echo "Activating environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$ENV_NAME" 2>/dev/null || conda activate "$ENV_PREFIX"

# Install everything via conda (conda first; pip only for genuine conda gaps).
# PyTorch + CUDA come from the pytorch/nvidia channels (float64-capable on H100).
# Match pytorch-cuda to the cluster's CUDA toolkit (see PYTORCH_CUDA above).
echo "Installing core packages + PyTorch (CUDA $PYTORCH_CUDA) via conda..."
conda install -y -c pytorch -c nvidia -c conda-forge \
    numpy \
    scipy \
    pandas \
    matplotlib \
    pytest \
    pytorch \
    "pytorch-cuda=${PYTORCH_CUDA}"

# Verify installation
echo ""
echo "=========================================="
echo "Verifying installation..."
echo "=========================================="
python -c "import numpy, scipy, pandas, matplotlib, pytest; print('Core packages OK')"
python -c "import torch; print('torch', torch.__version__, '| cuda', torch.cuda.is_available())"
python -c "import torch; x=torch.ones(3, dtype=torch.float64); print('float64 OK:', (x.sum().item()==3.0))"

echo ""
echo "=========================================="
echo "Environment setup complete!"
echo "=========================================="
echo ""
echo "Environment location: $ENV_PREFIX"
echo "Activate with:  conda activate $ENV_NAME   (or the full prefix path)"
echo ""
