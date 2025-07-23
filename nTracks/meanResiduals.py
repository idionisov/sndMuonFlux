import ROOT, uproot, argparse, roostyling
import numpy as np
import pandas as pd

from sndUtils import *
from ddfRoot import *
from ddfUtils import *

xy_eff_range = {
    'min': {'x': -42., 'y': 19.},
    'max': {'x': -10., 'y': 48.}
}

roostyling.setDdfpubStyle()

parser = argparse.ArgumentParser()

parser.add_argument('-r', '--runs', nargs="+", type=int, default=[7080])
parser.add_argument('-f', '--files', type=str, default='*f10*.root')
parser.add_argument('-i', '--input-dir', type=str, default='/eos/user/i/idioniso/1_Data/Tracks')
parser.add_argument('-o', '--fout', type=str, default='')
parser.add_argument('-z', '--z-ref', nargs="+", type=float, default=[430., 430., 450., 450.])
parser.add_argument('-xz', '--xz', type=float, default=80.)
parser.add_argument('-yz', '--yz', type=float, default=80.)
parser.add_argument('--xz-min', type=float, default=0.)
parser.add_argument('--xz-max', type=float, default=0.)
parser.add_argument('--yz-min', type=float, default=0.)
parser.add_argument('--yz-max', type=float, default=0.)
parser.add_argument('--scale', type=int, default=1)

args = parser.parse_args()
runs = args.runs
scale = args.scale
input = args.input_dir
fout = args.fout

TTs = (1, 11, 3, 13)
z_ref = {1: args.z_ref[0], 11: args.z_ref[1], 3: args.z_ref[2], 13: args.z_ref[3]}

if args.xz_min:
    xz_min = -abs(args.xz_min)
else:
    xz_min = -abs(args.xz)
if args.xz_max:
    xz_max = abs(args.xz_max)
else:
    xz_max = abs(args.xz)
if args.yz_min:
    yz_min = -abs(args.yz_min)
else:
    yz_min = -abs(args.yz)
if args.yz_max:
    yz_max =  abs(args.yz_max)
else:
    yz_max =  abs(args.yz)

data={}
tree={}
nEntries={}
for run in runs:
    if run!=8329:
        data[run] = SndData(Run=run, InputDir=input, Files=args.files)
    else:
        data[run] = SndData(Run=run, InputDir="/eos/experiment/sndlhc/users/sii/2024", TopDir=str(run), Files="*1*.root")
    tree[run] = data[run].Tree
    nEntries[run] = tree[run].GetEntries()
    data[run].Print()

start_time = time()
count = 0



def getMeanRes(trk):
    meanResidual = 0
    points = trk.GetPoints()
    for point in points:
        extrapolatedPoint = trk.GetPointAtZ(point.Z())
        meanResidual += abs(extrapolatedPoint.X() - point.X()) + abs(extrapolatedPoint.Y() - point.Y())
    meanResidual /= len(points)
    return meanResidual

h = {}
for run in runs:
    h[run] = {}
    for tt in TTs:
        h[run][tt] = ROOT.TH1F(f"h_meanRes_{tt}_{run}", ";Mean Residual [um];Counts", 240, 0, 12)

for i_run, run in enumerate(runs):
    print(f"Starting run {run}! [{i_run+1}/{len(runs)}]")

    for i_event, event in enumerate(tree[run]):
        count = printStatus(i_event, nEntries[run], start_time, count)

    

        for trk in event.Reco_MuonTracks:
            if trk.getTrackMom().Z()==0 or trk.getTrackFlag()==False:
                continue

            trk = DdfTrack(Track=trk, Event=event, IP1_Angle=1e12)

            if not trk.IsWithinAref(
                Zref=z_ref[trk.tt],
                xmin = xy_eff_range["min"]["x"],
                xmax = xy_eff_range["max"]["x"],
                ymin = xy_eff_range["min"]["y"],
                ymax = xy_eff_range["max"]["y"]
            ): continue

            aa=1e4*getMeanRes(trk)
            h[run][trk.tt].Fill(aa)


c = {}
for i_run, run in enumerate(runs):
    c[run] = ROOT.TCanvas(f"c_meanResiduals_{run}", "", 1000, 600)
    
    h[run][1].SetLineColor(60)
    h[run][1].Draw("hist")
    h[run][11].SetLineColor(206)
    h[run][11].Draw("hist same")
    h[run][3].SetLineColor(8)
    h[run][3].Draw("hist same")
    h[run][13].SetLineColor(13)
    h[run][13].Draw("hist same")



saveToRoot(c, h, fout="/eos/user/i/idioniso/mfout/trkMeanRes.root")
