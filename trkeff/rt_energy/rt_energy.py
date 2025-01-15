import glob
import pandas as pd
import uproot
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import argparse

from array import array
from ddf.pyfuncs import get_cl_sigma, get_current_dmy
from ddf.root.to_numpy import get_as_numpy
from ddf.root.to_pandas import get_as_pandas
from ddf.root.to_dict import get_as_dict
from ddf.root.misc import save_to_root, save_dict_to_root
from ddf.root.teff import get_teff_dict, get_graph_from_teff_1D
from ddf.snd.snd import get_N, print_status_with_time
from ddf.snd.mc import there_is_a_muon, get_intersection_mc, get_angle_xz, get_angle_yz, should_be_a_track
from ddf.snd.trk import sys, sys_name, alg, alg_name, n_name

from getters.histograms import *



def get_trees(
    track_types = (1, 11, 3, 13)
) -> dict:
    eff = {}

    efficiency = {}
    deff_up    = {}
    deff_low   = {}
    E          = {}
    dE_up      = {}
    dE_low     = {}

    for tt in track_types:
        eff[tt] = ROOT.TTree(f"eff_e_{tt}",  f"{sys_name(tt)} {alg_name(tt)} efficiency dependence on energy")

        efficiency[tt] = array('f', [ 0. ])
        deff_up[tt]    = array('f', [ 0. ])
        deff_low[tt]   = array('f', [ 0. ])
        E[tt]          = array('f', [ 0. ])
        dE_up[tt]      = array('f', [ 0. ])
        dE_low[tt]     = array('f', [ 0. ])

        eff[tt].Branch("efficiency", efficiency[tt], "efficiency/F")
        eff[tt].Branch("deff_up",    deff_up[tt],    "deff_up/F")
        eff[tt].Branch("deff_low",   deff_low[tt],   "deff_low/F")
        eff[tt].Branch("E",          E[tt],          "E/F")
        eff[tt].Branch("dE_up",      dE_up[tt],      "dE_up/F")
        eff[tt].Branch("dE_low_",    dE_low[tt],     "dE_low/F")

    return eff


int main():
    plt.style.use("root")

    parser = argparse.ArgumentParser()
    parser.add_argument('-z', '--z_ref', nargs="+", type=float, default=[430., 450., 430., 450.])
    parser.add_argument('-xz', '--xz', type=float, default=20.)
    parser.add_argument('-yz', '--yz', type=float, default=20.)
    parser.add_argument('-o', '--out', type=str, default="")
    parser.add_argument('-i', '--in', type=str, default="/eos/user/i/idioniso/mfout/muGun/reco")
    args = parser.parse_args()

    track_types=(1, 11, 3, 13)
    cl = get_cl_sigma(1)

    xy_range = {
        'min': {'x': -70., 'y': 0. },
        'max': {'x':  10., 'y': 80.}
    }
    xy_eff_range = {
        'min': {'x': -42., 'y': 19.},
        'max': {'x': -10., 'y': 48.}
    }

    track_types = (1, 11, 3, 13)
    z_ref = {1: 430., 11: 450., 3: 430., 13: 450.}

    e_arr = [10, 12.5, 15, 17.5, 20, 25, 30, 35, 55, 100, 200, 300, 450, 600, 1010]
    np_e_arr = np.array(e_arr, dtype=np.float64)

    d, m, y = get_current_dmy()
    if not args.out:
        outfile = f"/eos/user/i/idioniso/mfout/trkeff_muGun.rt_{d:02d}.{m:02d}.{y:04d}.root"
    else:
        if not args.out.endswith(".root"):
            outfile = f"{args.out}.root"
        else:
            outfile = args.out
    fout = ROOT.TFile(f"{args.out}/trkeff_mugun.rt_{d:02d}.{m:02d}.{y:04d}.root", "recreate")




if __name__ == "__main__":
    main()
