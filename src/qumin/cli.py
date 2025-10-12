
from __future__ import annotations
import argparse, os, sys, pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict

from .io import load_czi, max_projection
from .pipeline import Settings, make_masks, measure_intensity, partition_coefficient

STAGES = [
    "Unedited",
    "BackgroundRemoved",
    "Aggregates",
    "RedAggregates",
    "RedCytoplasm",
]

def _save_stage_png(out_dir: str, stem: str, label: str, img: np.ndarray, mask=None):
    import matplotlib.pyplot as plt
    os.makedirs(out_dir, exist_ok=True)
    plt.figure(figsize=(4,4))
    if mask is None:
        plt.imshow(img, cmap='gray')
    else:
        m = np.ma.masked_where(~mask, img)
        plt.imshow(img, cmap='gray')
        plt.imshow(m, cmap='magma', alpha=0.5)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{stem}_{label}.png"), dpi=200)
    plt.close()

def cmd_quantify(args):
    directory = pathlib.Path(args.dir)
    files = sorted([p for p in directory.iterdir() if p.suffix.lower()=='.czi'])
    if not files:
        print("No .czi files found in", directory)
        sys.exit(1)

    s = Settings(lower_k=args.lower_k, upper_k=args.upper_k, blur_sigma=args.blur, dilate_size=args.dilate)

    rows = []
    out_dir = directory / "qumin_output"
    out_dir.mkdir(exist_ok=True)
    for p in files:
        print(f"Processing {p.name}...")
        arr = load_czi(str(p))
        # arr shape: (Y,X,C)
        if args.main_channel >= arr.shape[-1] or args.red_channel >= arr.shape[-1]:
            print(f"ERROR: channel index out of range for {p.name}. Has {arr.shape[-1]} channels.")
            continue
        main_img = arr[..., args.main_channel]
        red_img = arr[..., args.red_channel]

        masks = make_masks(main_img, s)
        metrics = measure_intensity(main_img, red_img, masks)
        pc = partition_coefficient(metrics["Red Aggregates (Intensity)"], metrics["Red Cytoplasm (Intensity)"])

        row = {"Filename": p.name}
        row.update(metrics)
        row["PC"] = pc
        rows.append(row)

        if args.save_stages:
            stem = p.stem
            _save_stage_png(str(out_dir), stem, "Unedited", main_img, None)
            _save_stage_png(str(out_dir), stem, "BackgroundRemoved", main_img, masks["bg"])
            _save_stage_png(str(out_dir), stem, "Aggregates", main_img, masks["agg"])
            _save_stage_png(str(out_dir), stem, "RedAggregates", red_img, masks["agg"])
            cyto = np.logical_and(masks["bg"], ~masks["agg"])
            _save_stage_png(str(out_dir), stem, "RedCytoplasm", red_img, cyto)

    df = pd.DataFrame(rows)
    xlsx = directory / args.xlsx
    df.to_excel(xlsx, index=False)
    print("Saved:", xlsx)

def cmd_plot(args):
    import pandas as pd, seaborn as sns, matplotlib.pyplot as plt
    df = pd.read_excel(args.xlsx)
    if args.group_digits and args.group_digits > 0:
        df["Group"] = df["Filename"].str.slice(0, args.group_digits)
        x = "Group"
    else:
        x = "Filename"
    y = args.y
    plt.figure(figsize=(8,5))
    sns.boxplot(x=x, y=y, data=df)
    plt.title(args.title or y)
    plt.xticks(rotation=30, ha='right')
    out = args.title + ".png" if args.title else f"boxplot_{y}.png"
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    print("Saved:", out)
    if args.show:
        plt.show()
    plt.close()

def _collect_stage_pngs(directory: pathlib.Path) -> Dict[str, List[pathlib.Path]]:
    # Group by stem prefix before underscore (original stem) and map to stage files
    out = {}
    for p in (directory / "qumin_output").glob("*.png"):
        # Expect names like <stem>_<Stage>.png
        if "_" not in p.stem:
            continue
        stem, stage = p.stem.rsplit("_", 1)
        out.setdefault(stem, []).append(p)
    return out

def cmd_stitch(args):
    from PIL import Image
    directory = pathlib.Path(args.dir)
    groups = _collect_stage_pngs(directory)
    if not groups:
        print("No stage PNGs found in qumin_output. Run quantify with --save-stages first.")
        return
    out_dir = directory / (args.out_dir or "qumin_stitched")
    out_dir.mkdir(exist_ok=True)

    order = ["Unedited", "BackgroundRemoved", "Aggregates", "RedAggregates", "RedCytoplasm"]
    for stem, files in groups.items():
        # Map stage -> file
        mapping = {p.stem.rsplit("_",1)[1]: p for p in files}
        imgs = []
        for st in order:
            p = mapping.get(st)
            if p is None:
                continue
            imgs.append(Image.open(p).convert("RGB"))
        if not imgs:
            continue
        # Make a horizontal strip
        widths, heights = zip(*(im.size for im in imgs))
        canvas = Image.new("RGB", (sum(widths), max(heights)), "white")
        x = 0
        for im in imgs:
            canvas.paste(im, (x, 0))
            x += im.size[0]
        out_path = out_dir / f"stitched_{stem}.png"
        canvas.save(out_path)
        print("Saved:", out_path)

def cmd_sholl(args):
    try:
        from skimage.morphology import skeletonize
        from skan import Skeleton, sholl_analysis
    except Exception:
        print("Sholl requires 'skan' and 'scikit-image'. Install with: pip install skan scikit-image")
        sys.exit(1)

    arr = load_czi(args.file)
    red = arr[..., args.red_channel].astype(float)
    # simple mask from lower-k on red
    from .pipeline import compute_thresholds, Settings, make_masks
    s = Settings(lower_k=args.lower_k, upper_k=999, blur_sigma=args.blur, dilate_size=args.dilate)
    masks = make_masks(red, s)
    mask = masks["bg"]
    skel = skeletonize(mask)
    # pick center as image center for radii
    center = (red.shape[0]//2, red.shape[1]//2)
    sk = Skeleton(skel)
    radii = np.arange(5, max(red.shape)//2, 5)
    df = sholl_analysis(sk, center=center, radii=radii)
    print(df.head())
    out = os.path.splitext(args.file)[0] + "_sholl.csv"
    df.to_csv(out, index=False)
    print("Saved:", out)

def build_parser():
    p = argparse.ArgumentParser(prog="qumin", description="CZI → masks → metrics → figures")
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("quantify", help="Process .czi files and write Results.xlsx (+stage PNGs)")
    q.add_argument("--dir", default=".", help="Folder with .czi files")
    q.add_argument("--main-channel", type=int, required=True, dest="main_channel")
    q.add_argument("--red-channel", type=int, required=True, dest="red_channel")
    q.add_argument("--lower-k", type=float, default=0.25, dest="lower_k")
    q.add_argument("--upper-k", type=float, default=0.20, dest="upper_k")
    q.add_argument("--blur", type=float, default=0.0, help="Gaussian sigma; 0 disables")
    q.add_argument("--dilate", type=int, default=0, help="Dilation kernel size (odd int); 0 disables")
    q.add_argument("--save-stages", action="store_true")
    q.add_argument("--xlsx", default="Results.xlsx")
    q.set_defaults(func=cmd_quantify)

    b = sub.add_parser("plot", help="Make a boxplot from Results.xlsx")
    b.add_argument("--xlsx", default="Results.xlsx")
    b.add_argument("--group-digits", type=int, default=0)
    b.add_argument("--y", required=True, help="Column to plot (e.g., PC)")
    b.add_argument("--title", default="")
    b.add_argument("--show", action="store_true")
    b.set_defaults(func=cmd_plot)

    s = sub.add_parser("stitch", help="Stitch stage PNGs horizontally per image")
    s.add_argument("--dir", default=".")
    s.add_argument("--out-dir", default="qumin_stitched")
    s.set_defaults(func=cmd_stitch)

    sh = sub.add_parser("sholl", help="Run Sholl analysis on one file (optional skan)")
    sh.add_argument("--file", required=True)
    sh.add_argument("--red-channel", type=int, required=True, dest="red_channel")
    sh.add_argument("--lower-k", type=float, default=0.25)
    sh.add_argument("--blur", type=float, default=1.0)
    sh.add_argument("--dilate", type=int, default=3)
    sh.set_defaults(func=cmd_sholl)

    return p

def main(argv=None):
    p = build_parser()
    args = p.parse_args(argv)
    args.func(args)

if __name__ == "__main__":
    main()
