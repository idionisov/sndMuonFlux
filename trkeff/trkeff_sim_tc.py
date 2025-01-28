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
from helpers.mugun import xy_full_range, xy_eff_range, getTrees


def main():
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-z', '--z_ref', nargs="+", type=float, default=[430., 450., 430., 450.])
    parser.add_argument('-xz', '--xz', type=float, default=20.)
    parser.add_argument('-yz', '--yz', type=float, default=20.)
    parser.add_argument('-o', '--fout', type=str, default="")

    args = parser.parse_args()
    track_types = (1, 11, 3, 13)
    z_ref = {1: args.z_ref[0], 11: args.z_ref[1], 3: args.z_ref[2], 13: args.z_ref[3]}
    xz_min = -abs(args.xz)
    xz_max =  abs(args.xz)
    yz_min = -abs(args.yz)
    yz_max =  abs(args.yz)

    if args.fout:
        fout_name = args.fout
    else:
        fout_name = f"trkeff_sim_rt.root"
    mfout = "/eos/user/i/idioniso/mfout"
    fout_name = f"{mfout}/{fout_name}"

    inputDir = "/eos/user/i/idioniso/1_Data/Monte_Carlo"
    data = {
        'EMD': SndMCData(InputDir=inputDir, Files="muonReco_MC-EMD_PbPb.root", Geofile=f"{inputDir}/geofile_MC.PbPb.root"),
        'NI':  SndMCData(InputDir=inputDir, Files="muonReco_MC-NI_PbPb.root",  Geofile=f"{inputDir}/geofile_MC.PbPb.root")
    }

    mcFactor = {
        'EMD': 0.1388888888888889,
        'NI':  2.003205128205128
    }

    for k in data:
        data[k].InitGeo()
        data[k].Print()

    fout = TFile(fout_name, "recreate")

    h = createHists("sim", prpts=prpts_sim_all)
    passed = {1:0, 11:0, 3:0, 13:0}
    total  = {1:0, 11:0, 3:0, 13:0}

    start_time = time()
    count = 0
    nEntries = data["EMD"].Tree.GetEntries() + data["NI"].Tree.GetEntries()

    for key in data:
        tree = data[key].Tree

        if key == "NI":
            events = data["EMD"].Tree.GetEntries()
        else:
            events = 0

        for i_event, event in enumerate(tree):
            count = printStatus(i_event+events, nEntries, start_time, count)

            if not (event.EventHeader.isIP1() and thereIsAMuon(event)):
                continue

            _sf = sfTrackIsReconstructible(event)
            _ds = dsTrackIsReconstructible(event)

            if _sf==False and _ds==False:
                continue

            reco = {1: _sf, 11: _sf, 3: _ds, 13: _ds}
            weight = mcFactor[key] * event.MCTrack[0].GetWeight()

            for tt in track_types:
                if not reco[tt]:
                    continue

                total[tt] += weight
                for trk in event.Reco_MuonTracks:
                    if trk.getTrackType() != tt:
                        continue

                    passed[tt] += weight

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

                    flag = fillHistsRT(h, ddfMCTrack, "muGun.rt", z_ref[tt], tt, weight=weight)



    teff = getTEffDict(h, statOption='bayesian', suffix="rt")
    eq = getFitEq(teff, "sim.rt", track_types)
    saveToRoot(teff, fout=fout, nested=False, print_filename=True)

    eff = {}
    for tt in (1, 11, 3, 13):
        h[tt]["dxRef"].Write()
        h[tt]["dyRef"].Write()
        h[tt]["dxz"].Write()
        h[tt]["dyz"].Write()

        eff[tt] = {}

        eff[tt]['tree'] = TTree(f"eff_{tt}_sim.rt",  f"{system(tt)} {algorithm(tt)} efficiency dependence on energy")

        eff[tt]['eff']       = array('f', [ 0. ])
        eff[tt]['effErrUp']  = array('f', [ 0. ])
        eff[tt]['effErrLow'] = array('f', [ 0. ])

        eff[tt]['eff'][0], eff[tt]['effErrUp'][0], eff[tt]['effErrLow'][0] = getEffWithError(passed[tt], total[tt], statOption="bayesian")
        print(f" >> {tt}:\t{(eff[tt]['eff'][0], eff[tt]['effErrUp'][0], eff[tt]['effErrLow'][0])}")

        eff[tt]['tree'].Branch(f"eff",        eff[tt]['eff'],       "eff/F")
        eff[tt]['tree'].Branch("effErrUp",    eff[tt]['effErrUp'],  "effErrUp/F")
        eff[tt]['tree'].Branch("effErrLow",   eff[tt]['effErrLow'], "effErrLow/F")

        eff[tt]['tree'].Fill()
        eff[tt]['tree'].Write()

    fout.Close()



if __name__=="__main__":
    main()
