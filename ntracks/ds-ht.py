import argparse
from time import time
from ROOT import TChain, sndRecoTrack
from ddf.snd.snd import get_snd_years
from ddf.snd.trk import get_intersection, sys, alg
from ddf.pyfuncs import print_status


parser = argparse.ArgumentParser(description="Script for getting the number of tracks from IP1.")

parser.add_argument('--run', type=int, default=7080)
parser.add_argument('--selection', type=str, default='*')


args = parser.parse_args()
run = args.run
selection = args.selection

#xy_eff_range: dict = {
#    'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
#    'ds': {'min': {'x': -46., 'y': 14.}, 'max': {'x': -6.,  'y': 54.}}
#}

xy_eff_range: dict = {
    'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
    'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
}

cbmsim = TChain("cbmsim")
cbmsim.Add(f"/eos/user/i/idioniso/1_Data/Tracks/run_{run:06d}/{selection}.root")
nentries = cbmsim.GetEntries()

def is_within_xy_eff_range(
    track: sndRecoTrack,
    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.},
    xy_eff_range: dict = {
        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
        'ds': {'min': {'x': -46., 'y': 14.}, 'max': {'x': -6.,  'y': 54.}}
    }
) -> bool:
    tt = track.getTrackType()
    x = get_intersection(track, z_ref[tt]).X()
    y = get_intersection(track, z_ref[tt]).Y()

    if (
        x >= xy_eff_range[sys(tt)]['min']['x'] and
        x <= xy_eff_range[sys(tt)]['max']['x'] and
        y >= xy_eff_range[sys(tt)]['min']['y'] and
        y <= xy_eff_range[sys(tt)]['max']['y']
    ): return True
    else: return False

def is_ip1(track: sndRecoTrack, event: TChain) -> bool:
    xz = track.getAngleXZ()
    yz = track.getAngleYZ()

    if (
        event.EventHeader.isIP1() and
        xz <= 0.02 and xz >= -0.02 and
        yz <= 0.02 and yz >=-0.02
    ): return True
    else: return False

eventNr = 0

start_time = time()
for i_event, event in enumerate(cbmsim):
    if eventNr >= 100: break


    tt1=False
    tt11=False
    tt3=False
    tt13=False
    for i_trk, trk in enumerate(event.Reco_MuonTracks):
        if not trk.getTrackMom().Z(): continue

        tt = trk.getTrackType()
        xz = trk.getAngleXZ()
        yz = trk.getAngleYZ()

        if not (
            is_within_xy_eff_range(trk, xy_eff_range=xy_eff_range) and
            xz <= 0.02 and xz >= -0.02 and yz <= 0.02 and yz >= -0.02 
        ): continue
        
        if not is_ip1(trk, event): continue

        if tt==1:  tt1=True
        if tt==11: tt11=True
        if tt==3:  tt3=True
        if tt==13: tt13=True
    
    if (tt1==True and tt3==False) and (tt11==True and tt13==False):
        print(eventNr, i_event, "\tSciFi tracks of both types exist, but a DS one does not!")
        eventNr+=1
    elif tt1==True and tt3==False:
        print(eventNr, i_event, "\tSF ST track exists, but a DS ST does not!")
        eventNr+=1
    elif tt11==True and tt13==False:
        print(eventNr, i_event, "\tSF HT track exists, but a DS HT does not!")
        eventNr+=1

