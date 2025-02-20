import glob
from time import time
from array import array
from SndlhcGeo import GeoInterface
from ROOT import TH1F, TFile, TChain, TTree

from sndUtils import SndMCData, DdfTrack, DdfMCTrack, sfTrackIsReconstructible, dsTrackIsReconstructible, thereIsAMuon, system, algorithm
from ddfUtils import printStatus, getEffWithError
from ddfRoot import getTEffDict, saveToRoot

from helpers.hists import getHists, createHists, fillHistsRT, prpts_sim_all
from helpers.trkSelection import isSelected
from helpers.misc import getFitEq
from helpers.mugun import xy_full_range, xy_eff_range, getTrees, getRoundedE


def main():
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-e', '--energy', type=float, default=200)
    parser.add_argument('-z', '--z_ref', nargs="+", type=float, default=[430., 450., 430., 450.])
    parser.add_argument('-xz', '--xz', type=float, default=20.)
    parser.add_argument('-yz', '--yz', type=float, default=20.)
    parser.add_argument('-o', '--fout', type=str, default="")
    parser.add_argument('--remote-eos', type=bool, default=True)

    args = parser.parse_args()
    e = getRoundedE(args.energy)
    track_types = (1, 11, 3, 13)
    z_ref = {1: args.z_ref[0], 11: args.z_ref[1], 3: args.z_ref[2], 13: args.z_ref[3]}
    xz_min = -abs(args.xz)
    xz_max =  abs(args.xz)
    yz_min = -abs(args.yz)
    yz_max =  abs(args.yz)

    if args.fout:
        fout_name = args.fout
    else:
        fout_name = f"trkeff_muGun.{e}GeV_rt.root"

    if args.remote_eos:
        eos = "/eos/user/i/idioniso"
    else:
        eos = "/EOS/user/i/idioniso"

    mfout = f"{eos}/mfout"
    fout_name = f"{mfout}/trkeff/trkeff-E/{fout_name}"

    inputDir = f"{eos}/1_Data/Monte_Carlo/pGun/muons/{e}.GeV"
    data = SndMCData(InputDir=inputDir, Files=f"reco*{e}GeV.root")
    data.Print()

    fout = TFile(fout_name, "recreate")

    h = createHists("muGun", prpts=prpts_sim_all)
    passed = {1:0, 11:0, 3:0, 13:0}
    total  = {1:0, 11:0, 3:0, 13:0}

    start_time = time()
    count = 0
    tree = data.Tree
    nEntries = tree.GetEntries()
    a = {tt: {"passed": 0, "total": 0} for tt in (1, 11, 3, 13)}

    for i_event, event in enumerate(tree):
        count = printStatus(i_event, nEntries, start_time, count)
        if not (event.EventHeader.isIP1() and thereIsAMuon(event)):
            continue

        _sf = sfTrackIsReconstructible(event)
        _ds = dsTrackIsReconstructible(event)

        if _sf==False and _ds==False:
            continue

        reco = {1: _sf, 11: _sf, 3: _ds, 13: _ds}

        for tt in track_types:
            if not reco[tt]:
                continue

            total[tt] += 1

            for trk in event.Reco_MuonTracks:
                if trk.getTrackType() != tt:
                    continue

                passed[tt] += 1

            for mcTrack in event.MCTrack:

                ddfMCTrack = DdfMCTrack(mcTrack, Event=event, IP1_Angle=20.)
                if not (
                    ddfMCTrack.XZ <= xz_max/1e3 and
                    ddfMCTrack.XZ >= xz_min/1e3 and
                    ddfMCTrack.YZ <= yz_max/1e3 and
                    ddfMCTrack.YZ >= yz_min/1e3
                ): continue

                if not (ddfMCTrack.IsWithinDS3() and ddfMCTrack.IsWithinSF1()):
                    continue

                flag = fillHistsRT(h, ddfMCTrack, "muGun.rt", z_ref[tt], tt)
                if flag["passed"]: a[tt]["passed"] += 1
                if flag["total"]:  a[tt]["total"] += 1


    teff = getTEffDict(h, statOption='kfcp', suffix="rt")
    eq = getFitEq(teff, "muGun.rt", track_types)
    saveToRoot(teff, fout=fout, nested=False, print_filename=True)

    eff = {}
    for tt in (1, 11, 3, 13):
        h[tt]["dxRef"].Write()
        h[tt]["dyRef"].Write()
        h[tt]["dxz"].Write()
        h[tt]["dyz"].Write()

        eff[tt] = {}

        eff[tt]['tree'] = TTree(f"eff_{tt}_muGun.rt",  f"{system(tt)} {algorithm(tt)} efficiency dependence on energy")

        eff[tt]['eff']       = array('f', [ 0. ])
        eff[tt]['effErrUp']  = array('f', [ 0. ])
        eff[tt]['effErrLow'] = array('f', [ 0. ])

        eff[tt]['eff'][0], eff[tt]['effErrUp'][0], eff[tt]['effErrLow'][0] = getEffWithError(passed[tt], total[tt])
        print(f" >> {tt}:\t{(eff[tt]['eff'][0], eff[tt]['effErrUp'][0], eff[tt]['effErrLow'][0])}")

        a0, a0h, a0l = getEffWithError(a[tt]["passed"], a[tt]["total"])
        print(f" ~~ {tt}:\t{(a0, a0h, a0l)}")

        eff[tt]['tree'].Branch("eff",         eff[tt]['eff'],       "eff/F")
        eff[tt]['tree'].Branch("effErrUp",    eff[tt]['effErrUp'],  "effErrUp/F")
        eff[tt]['tree'].Branch("effErrLow",   eff[tt]['effErrLow'], "effErrLow/F")

        eff[tt]['tree'].Fill()
        eff[tt]['tree'].Write()

    fout.Close()



if __name__=="__main__":
    main()
