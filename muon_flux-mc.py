import argparse
import numpy as np
from time import time
from ROOT import TChain, sndRecoTrack
from ddf.snd.trk import get_intersection, sys, alg

parser = argparse.ArgumentParser(description="Script for evaluating simulation scaling factors and muon flux.")

parser.add_argument('--track_types', nargs='+', type=int, default=[1, 11, 3, 13])
parser.add_argument('--input_dir', type=str, default="/eos/user/i/idioniso/1_Data/Monte_Carlo")

args = parser.parse_args()

mfout = "/eos/user/i/idioniso/mfout"
track_types = args.track_types
input_dir = args.input_dir
xy = {
    'sf': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}},
    'ds': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}}
}
z_ref = {1: 430., 11: 450., 3: 430., 13: 450.}


input_files = {
    "EMD": f"{input_dir}/*EMD*PbPb*.root",
    "NI":  f"{input_dir}/*NI*PbPb*.root"
}
cbmsim = {
    "EMD": TChain("cbmsim"),
    "NI":  TChain("cbmsim")
}
area = {
    sys: abs(xy[sys]['max']['x']-xy[sys]['min']['x']) * abs(xy[sys]['max']['y']-xy[sys]['min']['y']) for sys in ("sf", "ds")
}


L_LHC = 6.4e-6
sigma = {"NI": 7.8e9, "EMD": 4.5e11}
N_rate = {"NI": 1e5, "EMD": 4e5}
L_MC = {i: N_rate[i]/sigma[i] for i in ("EMD", "NI")}

factor = {i: L_MC[i]/L_LHC for i in ("EMD", "NI")}

for i in ("EMD", "NI"):
    cbmsim[i].Add(input_files[i])
    print(f"  >> {i}:\t{cbmsim[i].GetEntries():,}")

eps = {
    1:  0.9576654730093384,
    11: 0.9543752351191904,
    3:  0.8495184882706716,
    13: 0.8480419728914257
}
deps_up = {
    1:  0.013751442782017786,
    11: 0.014977287008961948,
    3:  0.032926004277221255,
    13: 0.033125564174966815
}
deps_low = {
    1:  0.074738303590004,
    11: 0.07551826360610492,
    3:  0.05207095669326034,
    13: 0.052189537705818134
}


eventMCmu = {i: {} for i in ("EMD", "NI")}
Nrate = {
    i: {tt: 0 for tt in track_types} for i in ("EMD", "NI")
}


def is_ip1(track: sndRecoTrack, event: TChain) -> bool:
    tt = track.getTrackType()
    xz = track.getAngleXZ()
    yz = track.getAngleYZ()

    if (
        event.EventHeader.isIP1() and
        xz <= 0.02 and xz >= -0.02 and
        yz <= 0.02 and yz >= -0.02
    ): return True
    else: return False

def is_within_area(
    track: sndRecoTrack,
    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.},
    area: dict = {
        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
        'ds': {'min': {'x': -46., 'y': 14.}, 'max': {'x': -6.,  'y': 54.}}
    }
) -> bool:
    tt = track.getTrackType()
    ref = get_intersection(track, Z=z_ref[tt])

    if (
        ref.X() >= xy[sys(tt)]["min"]["x"] and
        ref.Y() >= xy[sys(tt)]["min"]["y"] and
        ref.X() <= xy[sys(tt)]["max"]["x"] and
        ref.Y() <= xy[sys(tt)]["max"]["y"]
    ): return True
    else: return False 

total_events = cbmsim["EMD"].GetEntries() + cbmsim["NI"].GetEntries()
for i in ("EMD", "NI"):
    for i_event, event in enumerate(cbmsim[i]):
        if i=="EMD": current_event = i_event
        else:        current_event = i_event + cbmsim["EMD"].GetEntries()

        if current_event%10000==0: print(f"{current_event*100/total_events:.02f}%")

        for mctrack in event.MCTrack:
            if mctrack.GetMotherId()==-1:
                eventMCmu[i][i_event] = mctrack.GetWeight()/factor[i]

        for trk in event.Reco_MuonTracks:   
            if (not trk.getTrackFlag() or trk.getTrackMom().Z()==0): continue
            
            tt = trk.getTrackType()
            ref = get_intersection(trk, Z=z_ref[tt])

            if (
                is_within_area(trk, z_ref=z_ref, area=xy) and
                is_ip1(trk, event)
            ): Nrate[i][tt] += eventMCmu[i][i_event]


mf = {}
for i in ("EMD", "NI"):
    mf[i] = {}

    for tt in track_types:
        A = area[sys(tt)]
        Nr = Nrate[i][tt]

        mf[i][tt] = (1e6*Nr)/(6.4*A*eps[tt])


dNrate = {
    i: {
        tt: np.sqrt(Nrate[i][tt]) for tt in track_types
    } for i in ("EMD", "NI")
}

dmf_up = {}
dmf_low = {}
for i in ("EMD", "NI"):
    dmf_up[i] = {}
    dmf_low[i] = {}

    for tt in track_types:
        A = area[sys(tt)]
        Nr = Nrate[i][tt]
        dNr = dNrate[i][tt]

        dmf_up[i][tt] = (1e6/(A*eps[tt]*6.4)) * np.sqrt(
            (Nr**2 * deps_up[tt]**2)/eps[tt]**2 + Nr
        )
        dmf_low[i][tt] = (1e6/(A*eps[tt]*6.4)) * np.sqrt(
            (Nr**2 * deps_low[tt]**2)/eps[tt]**2 + Nr
        )


for tt in track_types:
    flux = mf["EMD"][tt] + mf["NI"][tt]
    dflux_up = dmf_up["EMD"][tt] + dmf_up["NI"][tt]
    dflux_low = dmf_low["EMD"][tt] + dmf_low["NI"][tt]

    print(f"{tt}:\t{flux/1e4}\t±\t(+{dflux_up/1e4}, -{dflux_low/1e4})\t[10⁴ nb/cm²]")
