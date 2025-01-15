import argparse
from time import time
from ROOT import TChain, sndRecoTrack
from ddf.snd.snd import get_snd_years
from ddf.snd.trk import get_intersection, sys, alg, get_anti_tt
from ddf.pyfuncs import print_status


parser = argparse.ArgumentParser(description="Script for getting the number of tracks from IP1.")

parser.add_argument('--run', type=int, default=7080)
parser.add_argument('--selection', type=str, default='*')
parser.add_argument('--track_types', nargs='+', type=int, default=[1, 11, 3, 13])


args = parser.parse_args()
run = args.run
selection = args.selection
track_types = args.track_types

#xy_eff_range: dict = {
#    'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
#    'ds': {'min': {'x': -46., 'y': 14.}, 'max': {'x': -6.,  'y': 54.}}
#}

xy_eff_range: dict = {
    'sf': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}},
    'ds': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}}
}
z_ref = {1: 450., 11: 430., 3: 430., 13: 450.}

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
        xz <= 0.02 and xz >=-0.02 and
        yz <= 0.02 and yz >=-0.02
    ): return True
    else: return False


#def should_be_same_muon(
#    sf_trk: sndRecoTrack,
#    ds_trk: sndRecoTrack,
#    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.},
#    Aref: dict = {
#        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
#        'ds': {'min': {'x': -46., 'y': 14.}, 'max': {'x': -6.,  'y': 54.}}
#    }
#) -> bool:

#def is_deflected_out_of_Aref(
#    event: TChain,
#    ds_trk: sndRecoTrack,
#    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.},
#    Aref: dict = {
#        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
#        'ds': {'min': {'x': -46., 'y': 14.}, 'max': {'x': -6.,  'y': 54.}}
#    }
#) -> bool:
#    if (
#        is_ip1(ds_trk, event) and
#        is_within_xy_eff_range(ds_trk, z_ref, Aref)
#    ) or ds_trk.getTrackMom()==0: return False
#
#    ds_tt = ds_trk.getTrackType()
#
#    sf_trk = None
#    for trk in event.Reco_MuonTracks:
#        tt = trk.getTrackType()
#
#        if tt!=get_anti_tt(ds_tt)  or trk.getTrackMom().Z()==0: continue
#        sf_trk = trk
#        sf_tt = tt
#        
#    if sf_trk is None: return False
#    if not (
#        is_ip1(sf_trk, event) and
#        is_within_xy_eff_range(sf_trk, z_ref, Aref)
#    ): return False
#
#    x_ds = get_intersection(ds_trk, z_ref[tt]).X()
#    y_ds = get_intersection(ds_trk, z_ref[tt]).Y()
#    x_sf = get_intersection(sf_trk, z_ref[tt]).X()
#    y_sf = get_intersection(sf_trk, z_ref[tt]).Y()
#    
#    tt_ds = ds_trk.getTrackType()
#    tt_sf = sf_trk.getTrackType()
#
#    if (x_ds < Aref[sys(tt_ds)]['min']['x'] and x_ds >= Aref[sys(tt_ds)]['min']['x']-1):
#        if (x_sf > Aref[sys(tt_ds)]['min']['x'] and x_sf < Aref[sys])

ntracks = {
    'IP1': {
        tt: 0 for tt in track_types
    },
    'all': {
        tt: 0 for tt in track_types
    }
}
nevents = {'IP1': 0, 'all': nentries}

start_time = time()
count = 0
for i_event, event in enumerate(cbmsim):
    count = print_status(i_event, nentries, start_time, count)

    if event.EventHeader.isIP1(): nevents['IP1'] += 1

    for i_trk, trk in enumerate(event.Reco_MuonTracks):
        if not trk.getTrackMom().Z(): continue
        
        tt = trk.getTrackType()
        xz = trk.getAngleXZ()
        yz = trk.getAngleYZ()

        if not (
            is_within_xy_eff_range(trk, xy_eff_range=xy_eff_range) and
            xz <= 0.08 and xz >= -0.08 and yz <= 0.08 and yz >= -0.08 
        ): continue
        
        ntracks['all'][tt] += 1
        
        if is_ip1(trk, event):
            ntracks['IP1'][tt] += 1


print("\n\n" + 40*"-" + f"\nRUN {run}")
print("  TRACKS")
for i0, j0 in ntracks.items():
    print(f"    {i0}")
    for i1, j1 in j0.items():
        print(f"      {i1}:\t{j1}")

print("\n  EVENTS")
for i0, j0 in nevents.items():
    print(f"    {i0}:\t{j0}")
print(25*"-")
