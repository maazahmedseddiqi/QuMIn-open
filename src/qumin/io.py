
from __future__ import annotations
import numpy as np
from typing import Tuple
import czifile

def load_czi(path: str) -> np.ndarray:
    """Load a CZI file as a numpy array with channels last (Y,X,C) if possible.
    Strategy: read with czifile.imread; squeeze singleton dims; if dims > 3,
    we max-project across non-channel dims. Tries axes ordering heuristics.
    """
    arr = czifile.imread(path)
    arr = np.asarray(arr)
    # move channel axis to last if we can guess it
    # Heuristic: the smallest dim among last 4 is often C; but we avoid fantasy—
    # we check for dims <= 5 and assume C <= 5.
    axes = arr.shape
    # Flatten singleton dims
    while arr.ndim > 3:
        # Project over the first non-channel axis
        # We'll assume channel index is the last dim with size <= 5
        chan_idx = None
        for i, sz in enumerate(arr.shape[::-1]):
            if sz <= 5:
                chan_idx = arr.ndim - 1 - i
                break
        if chan_idx is None:
            # Fallback: assume last dim is channels if small; else project last axis
            chan_idx = arr.ndim - 1 if arr.shape[-1] <= 5 else None
        if arr.ndim == 4:
            # choose a projection axis that's not channel
            axes_order = list(range(arr.ndim))
            proj_axis = 0 if chan_idx != 0 else 1
            arr = arr.max(axis=proj_axis)
        else:
            # nd>4: project the first axis
            arr = arr.max(axis=0)
    # At this point arr should be (Y,X) or (Y,X,C) or (C,Y,X)
    if arr.ndim == 2:
        # single-channel
        arr = np.stack([arr], axis=-1)
    elif arr.shape[0] <= 5 and arr.ndim == 3:
        # (C,Y,X) -> (Y,X,C)
        arr = np.moveaxis(arr, 0, -1)
    return arr.astype(float)

def max_projection(arr: np.ndarray) -> np.ndarray:
    """Ensure array is (Y,X,C) float; if extra dims linger, max-project them."""
    a = np.asarray(arr, dtype=float)
    while a.ndim > 3:
        a = a.max(axis=0)
    if a.ndim == 2:
        a = np.stack([a], axis=-1)
    if a.shape[0] <= 5 and a.ndim == 3:
        a = np.moveaxis(a, 0, -1)
    return a
