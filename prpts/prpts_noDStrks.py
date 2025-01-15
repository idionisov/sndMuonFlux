import argparse
import numpy as np
from typing import Union
from ROOT import TFile, TH1F, TChain, TFile, sndRecoTrack
from time import time
from ddf.root.misc import save_to_root
from ddf.root.teff import calculate_eff
from ddf.snd.snd import get_N
from ddf.snd.trk import get_intersection, sys, alg, sys_name, alg_name
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
z_ref = {1: 440., 11: 435., 3: 460., 13: 475.}

if not args.input_dir:
    input_dir = f"/eos/user/i/idioniso/1_Data/Tracks/run_{run:06d}"
else:
    input_dir = args.input_dir

if not args.output_dir:
    output_dir = "/eos/user/i/idioniso/mfout"
else:
    output_dir = args.input_dir

if not args.fout:
    fout_name = f"prpts-no.DS.tracks_data_run{run}.root"
else:
    if args.fout.endswith(".root"): fout_name = args.fout
    else: fout_name = f"{args.fout}.root"

input_files = f"{input_dir}/{selection}.root"

cbmsim = TChain("cbmsim")
cbmsim.Add(input_files)
nevents = cbmsim.GetEntries()

fout = TFile(f"{output_dir}/{fout_name}", "recreate")

h = {
    tt: {
        prpt: get_h[prpt](run, tt) for prpt in prpts_all
    } for tt in (3, 13)
}


A: dict = {
    'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
    'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
}
z_ref = {1: 440., 11: 435., 3: 460., 13: 475.}


def is_within_A(
    track: sndRecoTrack,
    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.},
    A: dict = {
        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
        'ds': {'min': {'x': -46., 'y': 14.}, 'max': {'x': -6.,  'y': 54.}}
    }
) -> bool:
    tt = track.getTrackType()
    x = get_intersection(track, z_ref[tt]).X()
    y = get_intersection(track, z_ref[tt]).Y()

    if (
        x >= A[sys(tt)]['min']['x'] and
        x <= A[sys(tt)]['max']['x'] and
        y >= A[sys(tt)]['min']['y'] and
        y <= A[sys(tt)]['max']['y']
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


def with_missing_DS_track(
    event: TChain,
    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.},
    A: dict = {
        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
        'ds': {'min': {'x': -46., 'y': 14.}, 'max': {'x': -6.,  'y': 54.}}
    }
) -> dict:
    if not event.EventHeader.isIP1(): return {3: False, 13: False}

    trks     = {tt: None for tt in (1, 11, 3, 13)}
    lt20mrad = {tt: False for tt in (1, 11, 3, 13)}
    withinA  = {tt: False for tt in (1, 11, 3, 13)}

    has_missing_ds_track = {3: False, 13: False}

    for trk in event.Reco_MuonTracks:
        if not trk.getTrackMom().Z(): continue
        tt = trk.getTrackType()
        
        if   tt==1:  trks[1]=trk;  lt20mrad[1]
        elif tt==11: trks[11]=trk
        elif tt==3:  trks[3]=trk
        elif tt==13: trks[13]=trk

        if (trks[1] is not None and trks[3] is not None):
            xz = trks[1].getAngleXZ()
            yz = trks[1].getAngleYZ()

            if (
                event.EventHeader.isIP1() and
                xz <= 0.02 and xz >=-0.02 and
                yz <= 0.02 and yz >=-0.02
            ): lt20mrad[1]
        if (
            is_ip1(trks[1], event) and
            is_within_A(trks[1], z_ref, xy_eff_range)
        ) and not (
            is_ip1(trks[3], event) and
            is_within_A(trks[3], z_ref, xy_eff_range)
        ):
            has_missing_ds_track[3] = True
    
    if (trks[11] is not None and trks[13] is not None):
        if (
            is_ip1(trks[11], event) and
            is_within_A(trks[11], z_ref, xy_eff_range)
        ) and not (
            is_ip1(trk13, event) and
            is_within_A(trks[13], z_ref, xy_eff_range)
        ):
            has_missing_ds_track[13] = True

    return has_missing_ds_track


class PotentialTrack:
    def __init__(
        self, event: TChain, reco_trk: sndRecoTrack
    ):
        self.event = event
        self.reco_trk = reco_trk

    def isIP1(self):
        return self.event.EventHeader.isIP1()

    def isLt20mrad(self):
        xz = self.reco_trk.getAngleXZ()
        yz = self.reco_trk.getAngleYZ()
        
        if (abs(xz) <= 0.02 and abs(yz) <= 0.02): return True
        else: return False

    def isWithinA(self):
        tt = self.reco_trk.getTrackType()
        x = get_intersection(self.reco_trk, z_ref[tt]).X()
        y = get_intersection(self.reco_trk, z_ref[tt]).Y()

        if (
            x >= A[sys(tt)]['min']['x'] and
            x <= A[sys(tt)]['max']['x'] and
            y >= A[sys(tt)]['min']['y'] and
            y <= A[sys(tt)]['max']['y']
        ): return True
        else: return False



start_time = time()
count = 0

total    = {3: 0, 13: 0}
gt20mrad = {3: 0, 13: 0}
outsideA = {3: 0, 13: 0}
for i_event, event in enumerate(cbmsim):
    count = print_status(i_event, nevents, start_time, count)

    if not event.EventHeader.isIP1(): continue
    trks = {1: None, 11: None, 3: None, 13: None}

    for i_trk, trk in enumerate(event.Reco_MuonTracks):
        if trk.getTrackMom().Z()==0: continue
        tt = trk.getTrackType()

        trks[tt] = PotentialTrack(event, trk)

    if (trks[1] is not None and trks[3] is not None):
        if (
            trks[1].isIP1() and
            trks[1].isLt20mrad() and
            trks[1].isWithinA()
        ):
            total[3] += 1
            if   not trks[3].isLt20mrad(): gt20mrad[3] += 1
            elif not trks[3].isWithinA():  outsideA[3] += 1
    
    if (trks[11] is not None and trks[13] is not None):
        if (
            trks[11].isIP1() and
            trks[11].isLt20mrad() and
            trks[11].isWithinA()
        ):
            total[13] += 1
            if   not trks[13].isLt20mrad(): gt20mrad[13] += 1
            elif not trks[13].isWithinA():  outsideA[13] += 1

    for tt in (3, 13):
        if trks[tt] is None: continue

        if (
            trks[tt].isIP1() and trks[tt].isLt20mrad() and trks[tt].isWithinA()
        ): continue

        trk = trks[tt].reco_trk

        x  = get_intersection(trk, z_ref[tt]).X()
        y  = get_intersection(trk, z_ref[tt]).Y()
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

gt20={}; d_gt20_up={}; d_gt20_low={}
outA={}; d_outA_up={}; d_outA_low={}

for tt in (3, 13):
    gt20[tt], d_gt20_up[tt], d_gt20_low[tt] = calculate_eff(
        gt20mrad[tt], total[tt], stat_option="clopper_pearson"
    )
    outA[tt], d_outA_up[tt], d_outA_low[tt] = calculate_eff(
        outsideA[tt], total[tt], stat_option="clopper_pearson"
    )

noDs = {tt: gt20[tt] + outA[tt] for tt in (3, 13)}
dNoDs_up = {tt: d_gt20_up[tt] + d_outA_up[tt] for tt in (3, 13)}
dNoDs_low = {tt: d_gt20_low[tt] + d_outA_low[tt] for tt in (3, 13)}

f={}; df={}
for tt in (3, 13):
    f[tt] = outA[tt]/gt20[tt]

    df[tt] = (outsideA[tt]/gt20mrad[tt]) * np.sqrt( 1/outsideA[tt] + 1/gt20mrad[tt] )


print("\n\n" + 40*"-")
print(f"f         = Nᵒᵘᵗᴬ / Nᵍᵗ²⁰ᵐʳᵃᵈ")
print(f"Nˢᶜᵃᵗ     = Nᵍᵗ²⁰ᵐʳᵃᵈ + Nᵒᵘᵗᴬ")
print(f"ϵ         = Nˢᶜᵃᵗ / Nₜₒₜ")
print(f"ϵᵍᵗ²⁰ᵐʳᵃᵈ = Nᵍᵗ²⁰ᵐʳᵃᵈ / Nₜₒₜ")
print(f"ϵᵒᵘᵗᴬ     = Nᵒᵘᵗᴬ / Nₜₒₜ")

print(40*"-")
print("DS Simple tracking:\n")
print(f"  Nₜₒₜ:      \t{total[3]}")
print(f"  Nᵍᵗ²⁰ᵐʳᵃᵈ: \t{gt20mrad[3]}")
print(f"  Nᵒᵘᵗᴬ:     \t{outsideA[3]}")
print(f"  Nˢᶜᵃᵗ:     \t{gt20mrad[3] + outsideA[3]}")
print(f"  ϵ:         \t{noDs[3]} (+{dNoDs_up[3]}, -{dNoDs_low[3]})")
print(f"  ϵᵍᵗ²⁰ᵐʳᵃᵈ: \t{gt20[3]} (+{d_gt20_up[3]}, -{d_gt20_low[3]})")
print(f"  ϵᵒᵘᵗᴬ:     \t{outA[3]} (+{d_outA_up[3]}, -{d_outA_low[3]})")
print(f"  f:         \t{f[3]} ± {df[3]}")

print(40*"-")
print("DS Hough transform:\n")
print(f"  Nₜₒₜ:      \t{total[13]}")
print(f"  Nᵍᵗ²⁰ᵐʳᵃᵈ: \t{gt20mrad[13]}")
print(f"  Nᵒᵘᵗᴬ:     \t{outsideA[13]}")
print(f"  Nˢᶜᵃᵗ:     \t{gt20mrad[13] + outsideA[13]}")
print(f"  ϵ:         \t{noDs[13]} (+{dNoDs_up[13]}, -{dNoDs_low[13]})")
print(f"  ϵᵍᵗ²⁰ᵐʳᵃᵈ: \t{gt20[13]} (+{d_gt20_up[13]}, -{d_gt20_low[13]})")
print(f"  ϵᵒᵘᵗᴬ:     \t{outA[13]} (+{d_outA_up[13]}, -{d_outA_low[13]})")
print(f"  f:         \t{f[13]} ± {df[13]}")
print(40*"-" + "\n\n")
