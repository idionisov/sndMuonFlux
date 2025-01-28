import argparse
from typing import Union
from datetime import datetime
import ROOT
from time import time
from ddfRoot import saveToRoot
from sndUtils import getN, DdfTrack, SndMCData
from ddfUtils import printStatus
from helpers.hists import get_h_xz, get_h_yz, get_h_xzyz, get_h_x, get_h_y, get_h_xy, get_h_n, get_h_trkP, get_h_chi2, get_h_chi2ndf, get_h, prpts_all



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-tt', '--track-types', nargs='+', type=int, default=[1, 11, 3, 13])
    parser.add_argument('-i', '--input-dir', type=str, default="/eos/user/i/idioniso/1_Data/Monte_Carlo")
    parser.add_argument('-o', '--output-dir', type=str, default="/eos/user/i/idioniso/mfout")
    parser.add_argument('--fout', type=str, default="prptsSim.root")

    args = parser.parse_args()

    track_types = args.track_types
    z_ref = {1: 430., 11: 450., 3: 430., 13: 450.}
    inputDir = args.input_dir
    outputDir = args.output_dir
    fout = args.fout
    mcFactor = {
        'EMD': 0.1388888888888889,
        'NI':  2.003205128205128
    }

    data = {a: SndMCData(InputDir=inputDir, Files=f"muonReco_MC-{a}_PbPb.root") for a in ("EMD", "NI")}
    for a in data: data[a].Print()

    nEvents = data["EMD"].Tree.GetEntries() + data["NI"].Tree.GetEntries()

    fout = ROOT.TFile(f"{outputDir}/{fout}", "recreate")

    h = {
        tt: {
            prpt: get_h[prpt]("sim", tt) for prpt in prpts_all
        } for tt in track_types
    }


    start_time = time()
    count = 0

    for key in data:
        if key=="NI":
            events = data["EMD"].Tree.GetEntries()
        else:
            events = 0

        for i_event, event in enumerate(data[key].Tree):
            count = printStatus(i_event+events, nEvents, start_time, count)

            for i_trk, trk in enumerate(event.Reco_MuonTracks):
                if trk.getTrackMom().Z()==0:
                    continue

                weight = mcFactor[key] * event.MCTrack[0].GetWeight()

                trk = DdfTrack(Track=trk, Event=event)
                tt = trk.tt

                x  = trk.GetPointAtZ(Z = z_ref[tt]).X()
                y  = trk.GetPointAtZ(Z = z_ref[tt]).Y()
                xz = 1e3*trk.XZ
                yz = 1e3*trk.YZ
                chi2 = trk.Chi2
                chi2ndf = trk.Chi2Ndf
                trkP = trk.GetPoints().size()
                n = getN(tt, event=event)

                h[tt]["x"].Fill(x, weight)
                h[tt]["y"].Fill(y, weight)
                h[tt]["x.y"].Fill(x, y, weight)
                h[tt]["xz"].Fill(xz, weight)
                h[tt]["yz"].Fill(yz, weight)
                h[tt]["xz.yz"].Fill(xz, yz, weight)
                h[tt]["n"].Fill(n, weight)
                h[tt]["trkP"].Fill(trkP, weight)
                h[tt]["chi2"].Fill(chi2, weight)
                h[tt]["chi2ndf"].Fill(chi2ndf, weight)


    saveToRoot(h, fout=fout, nested=False, print_filename=True)
    fout.Close()

if __name__=="__main__":
    main()
