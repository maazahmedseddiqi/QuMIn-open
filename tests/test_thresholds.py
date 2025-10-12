
import numpy as np
from qumin.pipeline import compute_thresholds, make_masks, Settings, measure_intensity, partition_coefficient

def test_threshold_order():
    img = np.zeros((10,10), float)
    img[:,:5] = 0.0
    img[:,5:] = 10.0
    lower, upper = compute_thresholds(img, 0.0, 1.0)  # upper > lower
    assert upper > lower

def test_masks_and_metrics():
    main = np.zeros((10,10), float)
    main[:,5:] = 10.0
    red = main.copy()

    s = Settings(lower_k=0.1, upper_k=0.9, blur_sigma=0, dilate_size=0)
    masks = make_masks(main, s)
    m = measure_intensity(main, red, masks)

    assert m["Aggregates (Intensity)"] >= m["Background Removed (Intensity)"]
    pc = partition_coefficient(m["Red Aggregates (Intensity)"], m["Red Cytoplasm (Intensity)"])
    assert np.isnan(pc) or pc >= 0.0
