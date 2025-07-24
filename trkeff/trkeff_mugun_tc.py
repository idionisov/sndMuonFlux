import glob, os
import numpy as np
from time import time
from array import array
from SndlhcGeo import GeoInterface
from ROOT import TH1F, TFile, TChain, TTree

from sndUtils import SndMCData, DdfTrack, system, algorithm, att, getN
from ddfUtils import printStatus, getEffWithError
from ddfRoot import getTEffDict, saveToRoot

from helpers.hists import getHists, createHists, fillHistsTC, prpts_data_all
from helpers.trkSelection import isSelected
from helpers.misc import getFitEq, saveEffsMugun
from helpers.mugun import xy_full_range, xy_eff_range, getTrees, getRoundedE


def main():
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-e', '--energy', type=float, default=200)
    parser.add_argument('-z', '--z_ref', nargs="+", type=float, default=[430., 430., 450., 450.])
    parser.add_argument('-xz', '--xz', type=float, default=80.)
    parser.add_argument('-yz', '--yz', type=float, default=80.)
    parser.add_argument('--xz-min', type=float, default=0.)
    parser.add_argument('--xz-max', type=float, default=0.)
    parser.add_argument('--yz-min', type=float, default=0.)
    parser.add_argument('--yz-max', type=float, default=0.)
    parser.add_argument('-o', '--fout', type=str, default="")
    parser.add_argument('--remote-eos', type=bool, default=True)

    args = parser.parse_args()
    e = getRoundedE(args.energy)
    track_types = (1, 11, 3, 13)
    z_ref = {1: args.z_ref[0], 11: args.z_ref[1], 3: args.z_ref[2], 13: args.z_ref[3]}

    if args.xz_min > 1e-3:
        xz_min = -abs(args.xz_min)
    else:
        xz_min = -abs(args.xz)

    if args.xz_max > 1e-3:
        xz_max = abs(args.xz_max)
    else:
        xz_max = abs(args.xz)

    if args.yz_min > 1e-3:
        yz_min = -abs(args.yz_min)
    else:
        yz_min = -abs(args.yz)

    if args.yz_max > 1e-3:
        yz_max = abs(args.yz_max)
    else:
        yz_max = abs(args.yz)

    print(f"XZ min:\t{xz_min}")
    print(f"XZ max:\t{xz_max}")
    print(f"YZ min:\t{yz_min}")
    print(f"YZ max:\t{yz_max}")
    print(xy_eff_range)

    if args.fout:
        fout_name = args.fout
    else:
        fout_name = f"trkeff_muGun.{e}GeV_tc.root"

    if args.remote_eos:
        eos = "/eos/user/i/idioniso"
    else:
        eos = "/EOS/user/i/idioniso"

    mfout = f"{eos}/mfout"
    fout_name = f"{mfout}/trkeff/trkeff-E/{fout_name}"


    inputDir = f"{eos}/1_Data/Monte_Carlo/pGun/muons/{e}.GeV"
    geofile = f"{inputDir}/geofile_full.PG_13-TGeant4.root"

    data = SndMCData(InputDir=inputDir, Files=f"reco*{e}GeV.root", Geofile=geofile)
    data.InitGeo()
    data.Print()

    fout = TFile(fout_name, "recreate")


    h = createHists("muGun", prpts=prpts_data_all)
    passed = {tt: 0 for tt in (1, 11, 3, 13)}
    total  = {tt: 0 for tt in (1, 11, 3, 13)}

    start_time = time()
    count = 0
    tree = data.Tree
    nEntries = tree.GetEntries()
    for i_event, event in enumerate(tree):
        count = printStatus(i_event, nEntries, start_time, count)
        if not event.EventHeader.isIP1():
            continue

        for tt in track_types:
            for tag_trk in event.Reco_MuonTracks:
                if tag_trk.getTrackType() != att(tt):
                    continue



                tag_trk = DdfTrack(Track=tag_trk, Event=event, IP1_Angle=1e12)
                if not (
                    tag_trk.IsIP1() and
                    tag_trk.XZ >= xz_min/1e3 and tag_trk.XZ <= xz_max/1e3 and
                    tag_trk.YZ >= yz_min/1e3 and tag_trk.YZ <= yz_max/1e3
                ):
                    continue

                if tag_trk.tt==1 or tag_trk.tt==11:
                    if not (
                        tag_trk.IsWithinUS5Bar(data.Mufi, event.Digi_MuFilterHits) and
                        tag_trk.IsWithinDS3()
                    ): continue

                elif (tag_trk.tt==3 or tag_trk.tt==13):
                    if not tag_trk.IsWithinVetoBar(data.Mufi, event.Digi_MuFilterHits):
                        continue

                else: continue


                ref_tag = tag_trk.GetPointAtZ(z_ref[tt])
                x_tag = ref_tag.X()
                y_tag = ref_tag.Y()

                if (
                    ref_tag.X() >= xy_eff_range["min"]["x"] and
                    ref_tag.X() <= xy_eff_range["max"]["x"] and
                    ref_tag.Y() >= xy_eff_range["min"]["y"] and
                    ref_tag.Y() <= xy_eff_range["max"]["y"]
                ):
                    total[tt] += 1

                    for trk2 in event.Reco_MuonTracks:
                        trk2 = DdfTrack(Track=trk2, Event=event, IP1_Angle=1e9)

                        if not (
                            trk2.tt==tt and
                            trk2.IsIP1()
                        ): continue

                        xz_cand = 1e3*trk2.XZ
                        yz_cand = 1e3*trk2.YZ
                        ref2 = trk2.GetPointAtZ(z_ref[tt])

                        if not (
                            abs(ref_tag.X() - ref2.X()) <= 3. and
                            abs(ref_tag.Y() - ref2.Y()) <= 3.
                        ):
                            continue

                        passed[tt] += 1


                # flags = fillHistsTC(h, tag_trk, "muGun", z_ref[tt], tt)
                # print(flags)
                # if flags["passed"]:
                #     passed[tt] += 1
                # if flags["total"]:
                #     total[tt] += 1

    print(passed, total)

    # teff = getTEffDict(h, statOption='kfcp', suffix="tc")
    # eq = getFitEq(teff, "muGun.tc", track_types)
    # saveToRoot(teff, fout=fout, nested=False, print_filename=True)


    # for tt in (1, 11, 3, 13):
    #     h[tt]["dxRef"].Write()
    #     h[tt]["dyRef"].Write()
    #     h[tt]["dxz"].Write()
    #     h[tt]["dyz"].Write()

    saveEffsMugun(passed, total, fout=fout, statOption="clopper pearson", suffix="muGun.tc")
    print(fout_name)

    fout.Close()



if __name__=="__main__":
    main()
