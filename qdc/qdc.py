import argparse
from time import time
from ROOT import TFile, TChain, sndRecoTrack, TH1F, TH2F, TProfile
from ddf.snd.snd import get_snd_years
from ddf.snd.trk import get_intersection, sys, alg
from ddf.pyfuncs import print_status


parser = argparse.ArgumentParser(description="Script for getting the number of tracks from IP1.")

parser.add_argument('--run', type=int, default=7080)
parser.add_argument('--selection', type=str, default='*f10*')
parser.add_argument('--track_types', nargs='+', type=int, default=[1, 11, 3, 13])


args = parser.parse_args()
run = args.run
selection = args.selection
track_types = args.track_types

#A: dict = {
#    'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
#    'ds': {'min': {'x': -46., 'y': 14.}, 'max': {'x': -6.,  'y': 54.}}
#}

A: dict = {
    'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
    'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
}
z_ref = {1: 440., 11: 435., 3: 460., 13: 475.}


input_dir = f"/eos/user/i/idioniso/1_Data/Tracks/run_{run:06d}"
input_files = f"{input_dir}/{selection}.root"
print(f"INPUT:\t{input_files}")

output_dir = f"/eos/user/i/idioniso/mfout"
output_file = f"{output_dir}/qdc_Run{run}.root"
print(f"OUTPUT:\t{output_file}")


cbmsim = TChain("cbmsim")
cbmsim.Add(input_files)
nentries = cbmsim.GetEntries()

fout = TFile(output_file, "recreate")



h_qdc = TH1F(f"h_qdc_all_{run}", f"All IP1 events (run {run});QDC;", 400, 100, 4100);

h_qdc_neither = TH1F(f"h_qdc_neither_{run}", f"Events without missing DS tracks (run {run});QDC;Fraction of all entries", 400, 10, 4010)

h_qdc_both = TH1F(f"h_qdc_both_{run}", f"Events with missing DS tracks of both types (run {run});QDC;Fraction of all entries", 400, 10, 4010)
h_qdc_3    = TH1F(f"h_qdc_3_{run}", f"Events with missing DS ST tracks (run {run});QDC;Fraction of all Entries", 400, 10, 4010)
h_qdc_13   = TH1F(f"h_qdc_13_{run}", f"Events with missing DS HT tracks (run {run});QDC;Fraction of all entries", 400, 10, 4010)

h_dxz3_qdc  = TH2F(f"h_qdc.dxz_3_{run}",  f"Events with missing DS ST tracks (run {run});Difference in XZ angle d#theta_{{XZ}} [mrad];QDC;Fraction of all entries", 50, -50, 50, 100, 10, 4010)
h_dyz3_qdc  = TH2F(f"h_qdc.dyz_3_{run}",  f"Events with missing DS ST (run {run});Difference in YZ angle d#theta_{{YZ}} [mrad];QDC;Fraction of all entries", 50, -50, 50, 100, 10, 4010)
h_dxz13_qdc = TH2F(f"h_qdc.dxz_13_{run}", f"Events with missing DS HT (run {run});Difference in XZ angle d#theta_{{XZ}} [mrad];QDC;Fraction of all entries", 50, -50, 50, 100, 10, 4010)
h_dyz13_qdc = TH2F(f"h_qdc.dyz_13_{run}", f"Events with missing DS HT (run {run});Difference in YZ angle d#theta_{{YZ}} [mrad];QDC;Fraction of all entries", 50, -50, 50, 100, 10, 4010)

pr_dxz3_qdc  = TProfile(f"pr_dxz.qdc_3_{run}",  f"Events with missing DS ST tracks (run {run});Difference in XZ angle d#theta_{{XZ}} [mrad];Average QDC", 50, -50, 50)
pr_dyz3_qdc  = TProfile(f"pr_dyz.qdc_3_{run}",  f"Events with missing DS ST tracks (run {run});Difference in YZ angle d#theta_{{YZ}} [mrad];Average QDC", 50, -50, 50)
pr_dxz13_qdc = TProfile(f"pr_dxz.qdc_13_{run}", f"Events with missing DS HT tracks (run {run});Difference in XZ angle d#theta_{{XZ}} [mrad];Average QDC", 50, -50, 50)
pr_dyz13_qdc = TProfile(f"pr_dyz.qdc_13_{run}", f"Events with missing DS HT tracks (run {run});Difference in YZ angle d#theta_{{YZ}} [mrad];Average QDC", 50, -50, 50)

pr_qdc_dxz3  = TProfile(f"pr_qdc.dxz_3_{run}",  f"Events with missing DS ST tracks (run {run});QDC;Average |d#theta_{{XZ}}| [mrad]", 100, 10, 4010)
pr_qdc_dyz3  = TProfile(f"pr_qdc.dyz_3_{run}",  f"Events with missing DS ST tracks (run {run});QDC;Average |d#theta_{{YZ}}| [mrad]", 100, 10, 4010)
pr_qdc_dxz13 = TProfile(f"pr_qdc.dxz_13_{run}", f"Events with missing DS HT tracks (run {run});QDC;Average |d#theta_{{XZ}}| [mrad]", 100, 10, 4010)
pr_qdc_dyz13 = TProfile(f"pr_qdc.dyz_13_{run}", f"Events with missing DS HT tracks (run {run});QDC;Average |d#theta_{{YZ}}| [mrad]", 100, 10, 4010)

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

    trk1  = None
    trk11 = None
    trk3  = None
    trk13 = None
    has_missing_ds_track = {3: False, 13: False}

    for trk in event.Reco_MuonTracks:
        tt = trk.getTrackType()
        
        if   tt==1:  trk1=trk
        elif tt==11: trk11=trk
        elif tt==3:  trk3=trk
        elif tt==13: trk13=trk

    if (trk1 is not None and trk3 is not None):
        if (
            is_ip1(trk1, event) and
            is_within_A(trk1, z_ref, xy_eff_range)
        ) and not (
            is_ip1(trk3, event) and
            is_within_A(trk3, z_ref, xy_eff_range)
        ):
            has_missing_ds_track[3] = True
    
    if (trk11 is not None and trk13 is not None):
        if (
            is_ip1(trk11, event) and
            is_within_A(trk11, z_ref, xy_eff_range)
        ) and not (
            is_ip1(trk13, event) and
            is_within_A(trk13, z_ref, xy_eff_range)
        ):
            has_missing_ds_track[13] = True

    return has_missing_ds_track


def get_us_qdc(event: TChain) -> float:
    qdc = 0
    for mf_hit in event.Digi_MuFilterHits:
        if mf_hit.GetSystem() != 2: continue
        qdc += mf_hit.SumOfSignals()["Sum"]
    return qdc



start_time = time()
count = 0
for i_event, event in enumerate(cbmsim):
    count = print_status(i_event, nentries, start_time, count)

    if not event.EventHeader.isIP1(): continue

    qdc = get_us_qdc(event)
    h_qdc.Fill(qdc)




    trk1  = None
    trk11 = None
    trk3  = None
    trk13 = None
    no_ds_track = {3: False, 13: False}

    for trk in event.Reco_MuonTracks:
        tt = trk.getTrackType()
        
        if   tt==1:  trk1=trk
        elif tt==11: trk11=trk
        elif tt==3:  trk3=trk
        elif tt==13: trk13=trk

    if (trk1 is not None and trk3 is not None):
        if (
            is_ip1(trk1, event) and
            is_within_A(trk1, z_ref, A)
        ):
            if not (
                is_ip1(trk3, event) and
                is_within_A(trk3, z_ref, A)
            ): no_ds_track[3] = True

            elif (
                is_ip1(trk3, event) and
                is_within_A(trk3, z_ref, A)
            ): h_qdc_neither.Fill(qdc)
    
    if (trk11 is not None and trk13 is not None):
        if (
            is_ip1(trk11, event) and
            is_within_A(trk11, z_ref, A)
        ):
            if not (
                is_ip1(trk13, event) and
                is_within_A(trk13, z_ref, A)
            ): no_ds_track[13] = True
        
            elif(
                is_ip1(trk13, event) and
                is_within_A(trk13, z_ref, A)
            ): h_qdc_neither.Fill(qdc)

    if no_ds_track[3] and no_ds_track[13]:
        h_qdc_both.Fill(qdc)
        
    
    if no_ds_track[3]:
        h_qdc_3.Fill(qdc)
        
        dxz3 = 1e3 * (trk3.getAngleXZ() - trk1.getAngleXZ())
        dyz3 = 1e3 * (trk3.getAngleYZ() - trk1.getAngleYZ())

        h_dxz3_qdc.Fill(dxz3, qdc)
        pr_dxz3_qdc.Fill(dxz3, qdc)
        pr_qdc_dxz3.Fill(qdc, abs(dxz3))

        h_dyz3_qdc.Fill(dyz3, qdc)
        pr_dyz3_qdc.Fill(dyz3, qdc)
        pr_qdc_dyz3.Fill(qdc, abs(dyz3))

    if no_ds_track[13]:
        h_qdc_13.Fill(qdc)

        dxz13 = 1e3 * (trk13.getAngleXZ() - trk11.getAngleXZ())
        dyz13 = 1e3 * (trk13.getAngleYZ() - trk11.getAngleYZ())

        h_dxz13_qdc.Fill(dxz13, qdc)
        pr_dxz13_qdc.Fill(dxz13, qdc)
        pr_qdc_dxz13.Fill(qdc, abs(dxz13))

        h_dyz13_qdc.Fill(dyz13, qdc)
        pr_dyz13_qdc.Fill(dyz13, qdc)
        pr_qdc_dyz13.Fill(qdc, abs(dyz13))
    
h_qdc.Scale(1./h_qdc.Integral())
h_qdc_both.Scale(1./h_qdc_both.Integral())
h_qdc_neither.Scale(1./h_qdc_neither.Integral())
h_qdc_3.Scale(1./h_qdc_3.Integral())
h_qdc_13.Scale(1./h_qdc_13.Integral())

h_dxz3_qdc.Scale(1./h_dxz3_qdc.Integral())
h_dxz13_qdc.Scale(1./h_dxz13_qdc.Integral())
h_dyz3_qdc.Scale(1./h_dyz3_qdc.Integral())
h_dyz13_qdc.Scale(1./h_dyz13_qdc.Integral())


fout.Write()
fout.Close()
