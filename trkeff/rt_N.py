import glob
from time import time
from array import array
from SndlhcGeo import GeoInterface
from ROOT import TH1F, TFile, TChain, TTree

from sndUtils import SndMCData, DdfTrack, DdfMCTrack, sfTrackIsReconstructible, dsTrackIsReconstructible, thereIsAMuon, system, algorithm
from ddfUtils import printStatus, getEffWithError
from ddfRoot import getTEffDict, saveToRoot

from helpers.hists import getHists, createHists, fillHistsRT, getEffNRT, prpts_sim_all
from helpers.trkSelection import isSelected
from helpers.misc import getFitEq
from helpers.mugun import xy_full_range, xy_eff_range, getTrees, getRoundedE


def main():
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-z', '--z_ref', nargs="+", type=float, default=[430., 450., 430., 450.])
    parser.add_argument('-xz', '--xz', type=float, default=20.)
    parser.add_argument('-yz', '--yz', type=float, default=20.)
    parser.add_argument('-o', '--fout', type=str, default="")
    parser.add_argument('--remote-eos', type=bool, default=True)

    energies = [10, 12.5, 15, 17.5, 20, 25, 30, 35, 55, 100, 200, 300, 450, 600, 800, 1010]

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
        fout_name = f"trkeff_muGun_rt.root"

    if args.remote_eos:
        eos = "/eos/user/i/idioniso"
    else:
        eos = "/EOS/user/i/idioniso"

    mfout = f"{eos}/mfout"
    fout_name = f"{mfout}/trkeff/trkeff-E/{fout_name}"
    fout = TFile(fout_name, "recreate")

    h = createHists("muGun", prpts=("n"))
    passed = {1:0, 11:0, 3:0, 13:0}
    total  = {1:0, 11:0, 3:0, 13:0}

    start_time = time()
    count = 0
    for i_e, e in enumerate(energies):
        inputDir = f"{eos}/1_Data/Monte_Carlo/pGun/muons/{e}.GeV"
        print(2*"\n" + 25*"-" + f"[{i_e*100/len(energies):.02f}%] Starting dataset with energy {e} GeV")
        data = SndMCData(InputDir=inputDir, Files=f"reco*{e}GeV.root")
        data.Print()

        tree = data.Tree
        nEntries = tree.GetEntries()

        for i_event, event in enumerate(tree):
            count = printStatus(i_event, nEntries, start_time, count)
            if not (event.EventHeader.isIP1() and thereIsAMuon(event)):
                continue

            flag = getEffNRT(
                event = event,
                h = h,
                mcSet = "muGun.rt",
                z_ref = z_ref,
                weight = 1,
                track_types = (1, 11, 3, 13),
                xz_min = xz_min,
                yz_min = yz_min,
                xz_max = xz_max,
                yz_max = yz_max
            )
            for tt in (1, 11, 3, 13):
                if flag["passed"][tt]:
                    passed[tt] += 1
                if flag["total"][tt]:
                    total[tt] += 1

    teff = getTEffDict(h, statOption='kfcp', suffix="rt")
    saveToRoot(teff, fout=fout, nested=False, print_filename=True)

    eff = {}
    for tt in (1, 11, 3, 13):
        h[tt]["dxRef"].Write()
        h[tt]["dyRef"].Write()
        h[tt]["dxz"].Write()
        h[tt]["dyz"].Write()

        eff[tt] = {}

        eff[tt]["v"], eff[tt]["e+"], eff[tt]["e-"] = getEffWithError(passed[tt], total[tt])
        print(f" >> {tt}:\t{(eff[tt]['v'], eff[tt]['e-'], eff[tt]['e+'])}")

    fout.Close()



if __name__=="__main__":
    main()
