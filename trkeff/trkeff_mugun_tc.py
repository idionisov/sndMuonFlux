import glob, os
from time import time
from array import array
from SndlhcGeo import GeoInterface
from ROOT import TH1F, TFile, TChain, TTree

from sndUtils import SndMCData, DdfTrack, system, algorithm, att
from ddfUtils import printStatus, getEffWithError
from ddfRoot import getTEffDict, saveToRoot

from helpers.hists import getHists, createHists, fillHistsTC, prpts_data_all
from helpers.trkSelection import isSelected
from helpers.misc import getFitEq, saveEffsMugun
from helpers.mugun import xy_full_range, xy_eff_range, getTrees, getRoundedE


def main():
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--files', type=str, default='*')
    parser.add_argument('-e', '--energy', type=float, default=200)
    parser.add_argument('-z', '--z_ref', nargs="+", type=float, default=[430., 450., 430., 450.])
    parser.add_argument('-xz', '--xz', type=float, default=12.5)
    parser.add_argument('-yz', '--yz', type=float, default=12.5)
    parser.add_argument('-o', '--fout', type=str, default="")
    parser.add_argument('--remote-eos', type=bool, default=False)

    args = parser.parse_args()
    e = getRoundedE(args.energy)
    track_types = (1, 11, 3, 13)
    files = args.files
    z_ref = {1: args.z_ref[0], 11: args.z_ref[1], 3: args.z_ref[2], 13: args.z_ref[3]}
    xz_min = -abs(args.xz)
    xz_max =  abs(args.xz)
    yz_min = -abs(args.yz)
    yz_max =  abs(args.yz)

    if args.fout:
        fout_name = args.fout
    else:
        fout_name = f"trkeff_muGun.{e}GeV_tc.root"

    if args.remote_eos:
        eos = os.getenv("eos")
    else:
        eos = "/EOS/user/i/idioniso"

    mfout = f"{eos}/mfout"
    fout_name = f"{mfout}/trkeff/trkeff-E/{fout_name}"

    inputDir = f"{eos}/1_Data/Monte_Carlo/muGun/reco"
    geofile = f'{glob.glob(f"{eos}/1_Data/Monte_Carlo/muGun/sim/pGun_muons_{e}-{e}GeV_0-*")[-1]}/geofile_full.PG_13-TGeant4.root'

    data = SndMCData(InputDir=inputDir, Files=f"muon_reco_MC.gun_{e}-{e}GeV*.root", Geofile=geofile)
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


                tag_trk = DdfTrack(Track=tag_trk, Event=event, IP1_Angle=20.)
                if not tag_trk.IsIP1():
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

                flags = fillHistsTC(h, tag_trk, "muGun", z_ref[tt], tt)
                if flags["passed"]:
                    passed[tt] += 1
                if flags["total"]:
                    total[tt] += 1


    teff = getTEffDict(h, statOption='kfcp', suffix="tc")
    eq = getFitEq(teff, "muGun", track_types)
    saveToRoot(teff, fout=fout, nested=False, print_filename=True)


    for tt in (1, 11, 3, 13):
        h[tt]["dxRef"].Write()
        h[tt]["dyRef"].Write()
        h[tt]["dxz"].Write()
        h[tt]["dyz"].Write()

    saveEffsMugun(passed, total, fout=fout, statOption="clopper pearson", suffix="muGun.tc")


    fout.Close()



if __name__=="__main__":
    main()
