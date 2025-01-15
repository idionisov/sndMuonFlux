import argparse
from typing import Union
from ROOT import TFile, TH1F, TChain, TFile
from time import time
from ddf.root.misc import save_to_root
from ddf.snd.snd import get_N
from ddf.snd.trk import get_intersection
from ddf.pyfuncs import get_current_dmy, print_status
from hists import get_h_xz, get_h_yz, get_h_xzyz, get_h_x, get_h_y, get_h_xy, get_h_n, get_h_trkP, get_h_chi2, get_h_chi2ndf, get_h, prpts_all

parser = argparse.ArgumentParser(description="Script for plotting the track properties of an SND@LHC run")
parser.add_argument('--run', type=int, default=7080)
parser.add_argument('--selection', type=str, default="*f10*")
parser.add_argument('--track_types', nargs='+', type=int, default=[1, 11, 3, 13])
parser.add_argument('--input_dir', type=str, default="")
parser.add_argument('--output_dir', type=str, default="")
parser.add_argument('--fout', type=str, default="")

args = parser.parse_args()

run = args.run
selection = args.selection
track_types = args.track_types

if not args.input_dir:
    input_dir = f"/eos/user/i/idioniso/1_Data/Tracks/run_{run:06d}"
else:
    input_dir = args.input_dir

if not args.output_dir:
    output_dir = "/eos/user/i/idioniso/mfout"
else:
    output_dir = args.input_dir

if not args.fout:
    fout_name = f"prpts_data_run{run}.root"
else:
    fout_name = args.fout

input_files = f"{input_dir}/{selection}.root"

cbmsim = TChain("cbmsim")
cbmsim.Add(input_files)
nevents = cbmsim.GetEntries()

fout = TFile(f"{output_dir}/{fout_name}", "recreate")

h = {
    tt: {
        prpt: get_h[prpt](run, tt) for prpt in prpts_all
    } for tt in track_types
}


start_time = time()
count = 0

for i_event, event in enumerate(cbmsim):
    count = print_status(i_event, nevents, start_time, count)

    for i_trk, trk in enumerate(event.Reco_MuonTracks):
        if trk.getTrackMom().Z()==0: continue
        tt = trk.getTrackType()

        x  = get_intersection(trk, 490.).X()
        y  = get_intersection(trk, 490.).Y()
        xz = 1e3*trk.getAngleXZ()
        yz = 1e3*trk.getAngleYZ()
        chi2 = trk.getChi2()
        chi2ndf = trk.getChi2Ndf()
        trkP = trk.getTrackPoints().size()
        n = get_N(tt, event=event)

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


save_to_root(h, tfile=fout)
fout.Close()
