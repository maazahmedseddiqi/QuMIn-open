
# QUMIN — Quantitative Microscopy Intensity (CZI → masks → metrics)

**What it does (in plain English)**  
QUMIN loads **.czi** microscope images, creates two masks using thresholds from the image mean and standard deviation, and measures intensities:
- **Background-removed** mask: pixels > (mean + *k1*·std)
- **Aggregates** mask: pixels > (mean + *k2*·std)
It then measures:
- Mean intensities on the **main channel** for each mask (including unedited)
- Mean intensities on the **red channel** for **Aggregates** and **Cytoplasm** (= background-removed minus aggregates)
- A **Partition Coefficient (PC)** = Red Aggregates / Red Cytoplasm (guarded against divide-by-zero)

It also saves **stage PNGs** (Unedited, BackgroundRemoved, Aggregates, RedAggregates, RedCytoplasm), a **Results.xlsx** table, a **boxplot** from the table, and (optionally) a **Sholl analysis** of the red mask.

This repo is a clean, admissions-friendly package built from your original scripts (`1_qumin_aggregates.py`, `cool_figs.py`, `qumin_aggregates.py`, `rstitch_image.py`, `scholl.py`) while **preserving functionality**.

---

## 📦 Quickstart (local)
```bash
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt

# Process all .czi files in the current folder
python -m qumin quantify --main-channel 0 --red-channel 1 --lower-k 0.25 --upper-k 0.20 --save-stages

# Make a boxplot from Results.xlsx
python -m qumin plot --y PC --group-digits 2 --title Figure_1

# Stitch stage PNGs for each source file
python -m qumin stitch

# Optional: Sholl analysis for one file (requires 'skan')
python -m qumin sholl --file sample.czi --red-channel 1
```

**Tip:** If your CZI has Z/T/S dimensions, QUMIN uses a **max projection** over those dims before processing.

---

## 🔧 CLI reference

### `quantify`
```
python -m qumin quantify [--dir .] --main-channel 0 --red-channel 1 \
  [--lower-k 0.25] [--upper-k 0.20] [--blur 0] [--dilate 0] [--save-stages] [--xlsx Results.xlsx]
```
- `--lower-k` = coefficient for background mask threshold (mean + k·std)
- `--upper-k` = coefficient for aggregates threshold (mean + k·std)
- `--blur`    = Gaussian sigma (pixels); 0 to disable
- `--dilate`  = dilation size (odd kernel side length, e.g. 3, 5); 0 to disable

Writes `Results.xlsx` and stage PNGs if `--save-stages` is used.

### `plot`
```
python -m qumin plot [--xlsx Results.xlsx] [--group-digits 2] --y PC --title Figure_1
```
Adds a `Group` column based on the first `N` characters in filenames and saves a Seaborn boxplot.

### `stitch`
```
python -m qumin stitch [--out-dir stitched]
```
Builds a montage for each source using the stage PNGs.

### `sholl` (optional, requires `skan`)
```
python -m qumin sholl --file sample.czi --red-channel 1 [--lower-k 0.25] [--blur 1.0] [--dilate 3]
```

---

## 🧪 Tests & CI
- Unit tests run on array inputs (no microscope files required) so CI is fast and green.
- GitHub Actions workflow: `.github/workflows/python-tests.yml`

---

## 📄 License
Apache-2.0 (see LICENSE)

© 2025-10-12
