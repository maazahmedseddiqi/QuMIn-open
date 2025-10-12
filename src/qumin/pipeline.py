
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from typing import Dict, Tuple, Optional
from scipy.ndimage import gaussian_filter, binary_dilation
from skimage import morphology

@dataclass
class Settings:
    lower_k: float = 0.25
    upper_k: float = 0.20
    blur_sigma: float = 0.0  # in pixels; 0 disables
    dilate_size: int = 0     # odd integer kernel; 0 disables

def compute_thresholds(image: np.ndarray, lower_k: float, upper_k: float) -> Tuple[float, float]:
    """Return (lower_threshold, upper_threshold) as mean + k*std."""
    img = image.astype(float)
    std = img.std()
    mean = img.mean()
    lower = mean + lower_k * std
    upper = mean + upper_k * std
    return float(lower), float(upper)

def _prep(image: np.ndarray, s: Settings) -> np.ndarray:
    img = image.astype(float)
    if s.blur_sigma and s.blur_sigma > 0:
        img = gaussian_filter(img, sigma=s.blur_sigma)
    return img

def make_masks(main_channel_img: np.ndarray, s: Settings) -> Dict[str, np.ndarray]:
    """Create background-removed and aggregates masks.
    Returns a dict with keys: 'bg', 'agg' (boolean masks).
    """
    img = _prep(main_channel_img, s)
    lower, upper = compute_thresholds(img, s.lower_k, s.upper_k)
    bg = img > lower
    agg = img > upper
    if s.dilate_size and s.dilate_size >= 3 and s.dilate_size % 2 == 1:
        selem = morphology.square(s.dilate_size)
        bg = morphology.binary_dilation(bg, selem)
        agg = morphology.binary_dilation(agg, selem)
    return {"bg": bg.astype(bool), "agg": agg.astype(bool), "lower": lower, "upper": upper}

def _nonzero_mean(values: np.ndarray, mask: Optional[np.ndarray]) -> float:
    if mask is None:
        v = values
    else:
        v = values[mask]
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float('nan')
    return float(v.mean())

def measure_intensity(main_img: np.ndarray, red_img: np.ndarray, masks: Dict[str, np.ndarray]) -> Dict[str, float]:
    """Compute mean intensities for stages and PC."""
    bg = masks["bg"]; agg = masks["agg"]
    # Stages on main channel
    unedited = _nonzero_mean(main_img, None)
    bg_mean = _nonzero_mean(main_img, bg)
    agg_mean = _nonzero_mean(main_img, agg)
    # Red channel on regions
    red_agg = _nonzero_mean(red_img, agg)
    cytoplasm_mask = np.logical_and(bg, ~agg)
    red_cyto = _nonzero_mean(red_img, cytoplasm_mask)
    return {
        "Unedited (Intensity)": unedited,
        "Background Removed (Intensity)": bg_mean,
        "Aggregates (Intensity)": agg_mean,
        "Red Aggregates (Intensity)": red_agg,
        "Red Cytoplasm (Intensity)": red_cyto,
    }

def partition_coefficient(red_agg: float, red_cyto: float) -> float:
    if red_cyto is None or not np.isfinite(red_cyto) or red_cyto == 0:
        return float('nan')
    return float(red_agg / red_cyto)
