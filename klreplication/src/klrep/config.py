"""Global tensor backend config: device + dtype.

float64 is mandatory — it matches the Fortran `dp` (= IEEE double) and the model
solve loses accuracy over thousands of dampened iterations in float32. Apple MPS
is float32-only, so locally the float64 solve runs on CPU; the same device-
agnostic code targets CUDA (H100, float64-capable) later with no source change.

Override the device with the env var KLREP_DEVICE (e.g. "cpu", "cuda", "cuda:0").
"""
import os

import torch

# The one true dtype for the whole replication.
DTYPE = torch.float64


def get_device() -> torch.device:
    """Resolve the compute device. Honors KLREP_DEVICE; else CUDA if present;
    else CPU. MPS is deliberately NOT auto-selected (no float64 support)."""
    override = os.environ.get("KLREP_DEVICE", "").strip()
    if override:
        return torch.device(override)
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


DEVICE = get_device()
