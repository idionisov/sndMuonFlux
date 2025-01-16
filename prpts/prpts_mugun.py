import ROOT
import numpy as np
import os, glob
from time import time

from helpers.mugun import getRoundedE
from helpers.hists import prpts_all, prpts_all_sim, get_h
from sndUtils import SndMCData, DdfTrack
from ddfUtils import printStatus
from ddfRoot import getN, saveToRoot



import argparse
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--files', type=str, default='*')
parser.add_argument('-e', '--energy', type=float, default=200)
parser.add_argument('-o', '--fout', type=str, default="")
parser.add_argument('--remote-eos', type=bool, default=True)

args = parser.parse_args()
e = getRoundedE(args.energy)
track_types = (1, 11, 3, 13)
files = args.files
prpts = prpts_all

if args.remote_eos:
    eos = os.getenv("eos")
else:
    eos = "/EOS/user/i/idioniso"

if args.fout:
    fout_name = args.fout
else:
    fout_name = f"prpts_muGun.{e}GeV.root"
mfout = f"{eos}/mfout"
fout_name = f"{mfout}/prpts_varE/{fout_name}"

inputDir = f"{eos}/1_Data/Monte_Carlo/muGun/reco"

data = SndMCData(InputDir=inputDir, Files=f"muon_reco_MC.gun_{e}-{e}GeV*.root")
data.Print()

fout = ROOT.TFile(fout_name, "recreate")

h = {
    tt: {
        prpt: get_h[prpt]("muGun", tt) for prpt in prpts
    } for tt in track_types
}


start_time = time()
count = 0
tree = data.Tree
nEntries = tree.GetEntries()


for i_event, event in enumerate(tree):
    count = printStatus(i_event, nEntries, start_time, count)

    for i_trk, trk in enumerate(event.Reco_MuonTracks):
        if trk.getTrackMom().Z()==0: continue

        trk = DdfTrack(Track=trk, Event=event)
        tt=trk.tt

        x  = trk.GetPointAtZ(450.).X()
        y  = trk.GetPointAtZ(450.).Y()

        xz = 1e3*trk.XZ
        yz = 1e3*trk.YZ

        chi2 = trk.Chi2
        chi2ndf = trk.Chi2Ndf

        trkP = trk.sndRecoTrack.getTrackPoints().size()

        n = getN(tt, event=event)

        if "x"       in prpts: h[tt]["x"].Fill(x)
        if "y"       in prpts: h[tt]["y"].Fill(y)
        if "x.y"     in prpts: h[tt]["x.y"].Fill(x, y)
        if "xz"      in prpts: h[tt]["xz"].Fill(xz)
        if "yz"      in prpts: h[tt]["yz"].Fill(yz)
        if "xz.yz"   in prpts: h[tt]["xz.yz"].Fill(xz, yz)
        if "n"       in prpts: h[tt]["n"].Fill(n)
        if "trkP"    in prpts: h[tt]["trkP"].Fill(trkP)
        if "chi2"    in prpts: h[tt]["chi2"].Fill(chi2)
        if "chi2ndf" in prpts: h[tt]["chi2ndf"].Fill(chi2ndf)


saveToRoot(h, fout=fout, nested=False, print_filename=True)
fout.Close()
