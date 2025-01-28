import argparse
from typing import Union
from datetime import datetime
import ROOT
from time import time
from ddfRoot import saveToRoot
from sndUtils import getN, DdfTrack, SndData
from ddfUtils import printStatus
from helpers.hists import get_h_xz, get_h_yz, get_h_xzyz, get_h_x, get_h_y, get_h_xy, get_h_n, get_h_trkP, get_h_chi2, get_h_chi2ndf, get_h, prpts_all



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run', type=int, default=7080)
    parser.add_argument('-f', '--files', type=str, default="*f10*.root")
    parser.add_argument('-tt', '--track_types', nargs='+', type=int, default=[1, 11, 3, 13])
    parser.add_argument('-i', '--input-dir', type=str, default="")
    parser.add_argument('-o', '--output_dir', type=str, default="")
    parser.add_argument('--fout', type=str, default="")

    args = parser.parse_args()

    run = args.run
    files = args.files
    track_types = args.track_types
    z_ref = {1: 430., 11: 450., 3: 430., 13: 450.}

    if not args.input_dir:
        input_dir = f"/eos/user/i/idioniso/1_Data/Tracks"
    else:
        input_dir = args.input_dir

    if not args.output_dir:
        output_dir = "/eos/user/i/idioniso/mfout"
    else:
        output_dir = args.input_dir

    if not args.fout:
        fout_name = f"prptsDataRun{run}.root"
    else:
        fout_name = args.fout

    data = SndData(Run=run, InputDir=input_dir, Files=files)
    data.Print()

    cbmsim = data.Tree
    nevents = cbmsim.GetEntries()

    fout = ROOT.TFile(f"{output_dir}/{fout_name}", "recreate")

    h = {
        tt: {
            prpt: get_h[prpt](run, tt) for prpt in prpts_all
        } for tt in track_types
    }


    start_time = time()
    count = 0

    for i_event, event in enumerate(cbmsim):
        count = printStatus(i_event, nevents, start_time, count)

        for i_trk, trk in enumerate(event.Reco_MuonTracks):
            if trk.getTrackMom().Z()==0:
                continue

            trk = DdfTrack(Track=trk, Event=event)
            tt = trk.tt

            x  = trk.GetPointAtZ(Z = z_ref[tt]).X()
            y  = trk.GetPointAtZ(Z = z_ref[tt]).Y()
            xz = 1e3*trk.XZ
            yz = 1e3*trk.YZ
            chi2 = trk.Chi2
            chi2ndf = trk.Chi2Ndf
            trkP = trk.getPoints().size()
            n = getN(tt, event=event)

            h[tt]["x"].Fill(x)
            h[tt]["y"].Fill(y)
            h[tt]["xy"].Fill(x, y)
            h[tt]["xz"].Fill(xz)
            h[tt]["yz"].Fill(yz)
            h[tt]["xzyz"].Fill(xz, yz)
            h[tt]["n"].Fill(n)
            h[tt]["trkP"].Fill(trkP)
            h[tt]["chi2"].Fill(chi2)
            h[tt]["chi2ndf"].Fill(chi2ndf)


    saveToRoot(h, fout=fout, nested=False, print_filename=True)
    fout.Close()

if __name__=="__main__":
    main()
