import argparse
from SndlhcGeo import GeoInterface
from ROOT import TH1F, TFile, TChain

from sndUtils import SndData, DdfTrack

from ddf.root.misc import save_to_root
from ddf.root.teff import get_teff_dict
from ddf.pyfuncs import get_current_dmy, get_cl_sigma
from hists import get_hists, create_h
from process import loop, get_fit_eq

def main():
    parser = argparse.ArgumentParser(description="Script for calculating the tracking efficiency.")

    parser.add_argument('-r', '--run', type=int, default=7080)
    parser.add_argument('-s', '--selection', type=str, default='*')
    parser.add_argument('-z', '--z_ref', nargs="+", type=float, default=[430., 450., 430., 450.])
    parser.add_argument('-xz', '--xz', type=float, default=20.)
    parser.add_argument('-yz', '--yz', type=float, default=20.)
    parser.add_argument('-o', '--fout', type=str, default="")

    args = parser.parse_args()
    run = args.run
    track_types = (1, 11, 3, 13)
    selection = args.selection
    z_ref = {1: args.z_ref[0], 11: args.z_ref[1], 3: args.z_ref[2], 13: args.z_ref[3]}
    xz_min = {1: -abs(args.xz)/1e3, 11: -abs(args.xz)/1e3, 3: -abs(args.xz)/1e3, 13: -abs(args.xz)/1e3}
    xz_max = {1:  abs(args.xz)/1e3, 11:  abs(args.xz)/1e3, 3:  abs(args.xz)/1e3, 13:  abs(args.xz)/1e3}
    yz_min = {1: -abs(args.yz)/1e3, 11: -abs(args.yz)/1e3, 3: -abs(args.yz)/1e3, 13: -abs(args.yz)/1e3}
    yz_max = {1:  abs(args.yz)/1e3, 11:  abs(args.yz)/1e3, 3:  abs(args.yz)/1e3, 13:  abs(args.yz)/1e3}

    if args.fout:
        fout_name = args.fout
    else:
        fout_name = f"trkeff_Run{run}_tc.root"
    mfout = "/eos/user/i/idioniso/mfout"
    fout_name = f"{mfout}/{fout_name}"

    cl = get_cl_sigma(1)
    print("\n", fout_name)
    print(f"Zref: {z_ref}")

    data = SndData(Run=run, Selection=selection)
    data.Print()

    cbmsim = TChain("cbmsim")
    input_files = f"/eos/user/i/idioniso/1_Data/Tracks/run_{run:06d}/{selection}.root"
    print(f"Input files:\t{input_files}")
    cbmsim.Add(input_files)
    print(f"Total events:\t{cbmsim.GetEntries():,}")

    fout = TFile(fout_name, "recreate")

    geo   = GeoInterface("$SND_DATA/2023/geofile_sndlhc_TI18_V4_2023.root")
    scifi = geo.modules['Scifi']
    mufi  = geo.modules['MuFilter']

    h = create_h(run)


    loop(
        hists = h,
        tree = cbmsim,
        scifi = scifi,
        mufi = mufi,
        run = run,
        track_types = track_types,
        z_ref = z_ref,
        xz_min = xz_min,
        xz_max = xz_max,
        yz_min = yz_min,
        yz_max = yz_max
    )

    teff = get_teff_dict(h, stat_option='kfcp', confidence_level=cl, suffix="tc")


    eq = get_fit_eq(teff, run, track_types)


    save_to_root(teff, tfile=fout)

    for tt in (1, 11, 3, 13):
        h[tt]["dxRef"].Write()
        h[tt]["dyRef"].Write()
        h[tt]["dxz"].Write()
        h[tt]["dyz"].Write()

    fout.Close()


if __name__ == "__main__":
    main()
