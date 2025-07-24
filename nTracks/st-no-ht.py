import ROOT
import numpy as np
import pandas as pd
import argparse


from time import time
from sndUtils import *
from ddfUtils import *
from ddfRoot import *

parser = argparse.ArgumentParser()

parser.add_argument('-r', '--run', type=int, default=7080)
parser.add_argument('-f', '--files', type=str, default='*f10*.root')
parser.add_argument('-i', '--input-dir', type=str, default='/eos/user/i/idioniso/1_Data/Tracks')
parser.add_argument('-z', '--z_ref', nargs="+", type=float, default=[430., 430., 450., 450.])
parser.add_argument('-xz', '--xz', type=float, default=80.)
parser.add_argument('-yz', '--yz', type=float, default=80.)
parser.add_argument('-c', '--cut', type=int, default=1000)
parser.add_argument('--chi2ndf', nargs="+", type=float, default=[1e12, 1e12, 1e12, 1e12])

args = parser.parse_args()
run = args.run
input_dir = args.input_dir
track_types = (1, 11, 3, 13)
files = args.files
z_ref = {1: args.z_ref[0], 11: args.z_ref[1], 3: args.z_ref[2], 13: args.z_ref[3]}
chi2ndf = {1: args.chi2ndf[0], 11: args.chi2ndf[1], 3: args.chi2ndf[2], 13: args.chi2ndf[3]}
xz_min = -abs(args.xz)
xz_max =  abs(args.xz)
yz_min = -abs(args.yz)
yz_max =  abs(args.yz)
cutoff = args.cut

if run!=8329:
    data = SndData(Run=run, InputDir=input_dir, Files=files)
else:
    data = SndData(Run=run, InputDir="/eos/experiment/sndlhc/users/sii/2024", TopDir=str(run), Files=files)

data.Print()


EVENTS = {"sf": [], "ds": []}
start_time = time()
count = 0
tree = data.Tree
nEntries = tree.GetEntries()
for i_event, event in enumerate(tree):
    if i_event==cutoff:
        break


    count = printStatus(i_event, nEntries, start_time, count)
    if not event.EventHeader.isIP1():
        continue



    trkExists = {tt: False for tt in (1, 11, 3, 13)}

    for trk in event.Reco_MuonTracks:
        trk = DdfTrack(Track=trk, Event=event)

        if not trk.IsIP1():
            trkExists[trk.tt] = False
            continue

        ref = trk.GetPointAtZ(z_ref[trk.tt])
        x = ref.X()
        y = ref.Y()

        #if not (
        #    x < -10. and x > -42. and
        #    y <  48. and y > 19. and
        #    trk.XZ <= xz_max/1e3 and
        #    trk.XZ >= xz_min/1e3 and
        #    trk.YZ <= yz_max/1e3 and
        #    trk.YZ >= yz_min/1e3
        #):
        #    trkExists[trk.tt] = False
        #    continue

        trkExists[trk.tt] = True

    #print(trkExists)
    eventNum = event.EventHeader.GetEventNumber()
    if trkExists[1] and not trkExists[11]:   
        EVENTS["sf"].append(eventNum)
        print(eventNum, trkExists)
    
    if trkExists[3] and not trkExists[13]:
        EVENTS["ds"].append(eventNum)


