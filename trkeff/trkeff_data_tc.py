from time import time
from SndlhcGeo import GeoInterface
from ROOT import TH1F, TFile, TChain

from sndUtils import SndData, DdfTrack, att
from ddfUtils import printStatus
from ddfRoot import getTEffDict, saveToRoot

from helpers.hists import getHists, createHists, fillHistsTC, prpts_data_all
from helpers.trkSelection import isSelected
from helpers.misc import getFitEq, saveEffsData



def main():
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-r', '--run', type=int, default=7080)
    parser.add_argument('-f', '--files', type=str, default='*f10*.root')
    parser.add_argument('-z', '--z_ref', nargs="+", type=float, default=[430., 450., 430., 450.])
    parser.add_argument('-xz', '--xz', type=float, default=20.)
    parser.add_argument('-yz', '--yz', type=float, default=20.)
    parser.add_argument('-o', '--fout', type=str, default="")
    parser.add_argument('--chi2ndf', nargs="+", type=float, default=[1e6, 1e6, 1e6, 1e6])

    args = parser.parse_args()
    run = args.run
    track_types = (1, 11, 3, 13)
    files = args.files
    z_ref = {1: args.z_ref[0], 11: args.z_ref[1], 3: args.z_ref[2], 13: args.z_ref[3]}
    chi2ndf = {1: args.chi2ndf[0], 11: args.chi2ndf[1], 3: args.chi2ndf[2], 13: args.chi2ndf[3]}
    xz_min = -abs(args.xz)
    xz_max =  abs(args.xz)
    yz_min = -abs(args.yz)
    yz_max =  abs(args.yz)

    mfout = "/eos/user/i/idioniso/mfout"
    if args.fout:
        fout_name = args.fout
    else:
        fout_name = f"{mfout}/trkeff_Run{run}_tc.root"


    #geofile = "/eos/experiment/sndlhc/convertedData/physics/2023/geofile_sndlhc_TI18_V4_2023.root"
    geofile = "/eos/experiment/sndlhc/convertedData/physics/2024/run_2412/geofile_sndlhc_TI18_V12_2024.root"
    data = SndData(Run=run, InputDir="/eos/user/i/idioniso/1_Data/Tracks", Files=files, Geofile=geofile)
    data.InitGeo()
    data.Print()

    fout = TFile(fout_name, "recreate")

    h = createHists(run, prpts=prpts_data_all)

    start_time = time()
    count = 0
    tree = data.Tree
    nEntries = tree.GetEntries()

    for i_event, event in enumerate(tree):
        count = printStatus(i_event, nEntries, start_time, count)
        if not event.EventHeader.isIP1():
            continue

        eventNum = event.EventHeader.GetEventNumber()
        for tt in track_types:
            for tag_trk in event.Reco_MuonTracks:
                tag_trk = DdfTrack(Track=tag_trk, Event=event, IP1_Angle=0.02)

                if not (tag_trk.tt != att(tt) and tag_trk.IsGood(xz_min=-0.02, xz_max=0.02, yz_min=-0.02, yz_max=0.02)):
                    continue

                ref_tag = tag_trk.GetPointAtZ(z_ref[tt])
                x_tag = ref_tag.X()
                y_tag = ref_tag.Y()

                if not (
                    x_tag < -10. and x_tag > -42. and
                    y_tag <  48. and y_tag > 19.
                ):
                    continue



                if tag_trk.tt==1 or tag_trk.tt==11:
                    if not (
                        tag_trk.IsWithinUS5Bar(data.Mufi, event.Digi_MuFilterHits) and
                        tag_trk.IsWithinDS3() and
                        tag_trk.Chi2Ndf <= chi2ndf[tag_trk.tt]
                    ): continue

                elif (tag_trk.tt==3 or tag_trk.tt==13):
                    if not (
                        tag_trk.IsWithinVetoBar(data.Mufi, event.Digi_MuFilterHits) and
                        tag_trk.Chi2Ndf <= chi2ndf[tag_trk.tt]
                    ): continue

                else: continue


                flags = fillHistsTC(h, tag_trk, run, z_ref[tt], tt)



    teff = getTEffDict(h, statOption='kfcp', suffix="tc")
    eq = getFitEq(teff, run, track_types)
    saveToRoot(teff, fout=fout, nested=False, print_filename=True)

    for tt in (1, 11, 3, 13):
        h[tt]["dxRef"].Write()
        h[tt]["dyRef"].Write()
        h[tt]["dxz"].Write()
        h[tt]["dyz"].Write()


    saveEffsData(teff,
        fout=fout,
        statOption = "clopper pearson",
        suffix = "data.tc"
    )

    fout.Close()



if __name__=="__main__":
    main()
