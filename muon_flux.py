import numpy as np
import ROOT, argparse, os, glob
from pathlib import Path
from datetime import datetime


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run', type=int, default=7080)
    parser.add_argument('-s', '--selection', type=str, default="*f10*")
    parser.add_argument('--xmin_sf', type=float, default=-42.)
    parser.add_argument('--xmax_sf', type=float, default=-10.)
    parser.add_argument('--xmin_ds', type=float, default=-42.)
    parser.add_argument('--xmax_ds', type=float, default=-10.)
    parser.add_argument('--ymin_sf', type=float, default=19.)
    parser.add_argument('--ymax_sf', type=float, default=48.)
    parser.add_argument('--ymin_ds', type=float, default=19.)
    parser.add_argument('--ymax_ds', type=float, default=48.)

    args = parser.parse_args()
    run = args.run
    selection = args.selection
    xy = {
        "sf": {
            "min": {"x": args.xmin_sf, "y": args.ymin_sf},
            "max": {"x": args.xmax_sf, "y": args.ymax_sf}
        },
        "ds": {
            "min": {"x": args.xmin_ds, "y": args.ymin_ds},
            "max": {"x": args.xmax_ds, "y": args.ymax_ds}
        }
    }
    A = {sys: abs(xy[sys]["max"]["x"] - y[sys]["min"]["x"]) * abs(xy[sys]["max"]["y"] - y[sys]["min"]["y"]) for sys in ("sf", "ds")}


if __name__=="__main__":
    main()
