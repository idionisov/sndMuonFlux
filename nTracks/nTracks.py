import os
from time import time
import ROOT
import numpy as np
import pandas as pd
import uproot


from sndUtils import SndData, DdfTrack
from ddfUtils import printStatus
from helpers import getEosDir, getInputDir, getOutput, saveToCsv, xy_eff_range



def updateNtracks(rootFile: str = "/eos/user/i/idioniso/mfout/MuonFlux.root"):
    df = pd.read_csv("/eos/user/i/idioniso/mfout/nTracks.csv")

    with uproot.update(rootFile) as f:
        f["nTracks"] = df


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Script for getting the number of tracks from IP1.")

    parser.add_argument('-r', '--run', type=int, default=7080)
    parser.add_argument('-f', '--files', type=str, default='*')
    parser.add_argument('-i', '--input-dir', type=str, default='')
    parser.add_argument('-o', '--fout', type=str, default='')
    parser.add_argument('-z', '--z-ref', nargs="+", type=float, default=[450., 430., 450., 430.])
    parser.add_argument('-xz', '--xz', type=float, default=80.)
    parser.add_argument('-yz', '--yz', type=float, default=80.)
    parser.add_argument('--remote-eos', type=bool, default=True)
    parser.add_argument('--scale', type=int, default=1)

    args = parser.parse_args()

    run = args.run
    scale = args.scale
    input = getInputDir(args.input_dir, args.remote_eos)
    fout = getOutput(args.fout, args.remote_eos)


    TTs = (1, 11, 3, 13)
    z_ref = {1: args.z_ref[0], 11: args.z_ref[1], 3: args.z_ref[2], 13: args.z_ref[3]}
    xz_min = -abs(args.xz)
    xz_max =  abs(args.xz)
    yz_min = -abs(args.yz)
    yz_max =  abs(args.yz)

    data = SndData(Run=run, InputDir=input, Files=args.files)
    data.Print()

    start_time = time()
    count = 0
    tree = data.Tree
    nEntries = tree.GetEntries()


    nTracks = {
        'IP1': {tt: 0 for tt in TTs},
        'all': {tt: 0 for tt in TTs}
    }

    for i_event, event in enumerate(tree):
        count = printStatus(i_event, nEntries, start_time, count)

        for i_trk, trk in enumerate(event.Reco_MuonTracks):
            if trk.getTrackMom().Z()==0:
                continue
            trk = DdfTrack(Track=trk, Event=event, IP1_Angle=xz_max)

            if not trk.IsWithinAref(
                Zref=z_ref[trk.tt],
                xmin = xy_eff_range["min"]["x"],
                xmax = xy_eff_range["max"]["x"],
                ymin = xy_eff_range["min"]["y"],
                ymax = xy_eff_range["max"]["y"]
            ): continue

            tt = trk.tt
            xz = trk.XZ
            yz = trk.YZ

            nTracks['all'][tt] += 1

            if (
                trk.IsIP1() and
                trk.XZ >= xz_min/1e3 and
                trk.XZ <= xz_max/1e3 and
                trk.YZ >= yz_min/1e3 and
                trk.YZ <= yz_max/1e3
            ):
                nTracks['IP1'][tt] += 1


    header = ["run"]
    data = [[run]]
    for tt in TTs:
        for a in "IP1", "all":
            if a == "IP1":
                header.append(f"nTracks{tt}")
            else:
                header.append(f"nTracks{tt}Total")
            data[0].append(nTracks[a][tt])
    header.append("scale")
    data[0].append(scale)

    saveToCsv(header, data, fout)
    #updateNtracks("/eos/user/i/idioniso/mfout/MuonFlux.root")


if __name__=="__main__":
    main()
