import argparse
from ROOT import TChain, TFile, TH1F, sndRecoTrack, gPad, TCanvas, kRed, kBlue, TLegend
from ddf.snd.trk import get_intersection, sys, alg, sys_name, alg_name
from ddf.pyfuncs import print_status
from time import time

parser = argparse.ArgumentParser(description="Script for getting angle differences between SciFi and DS")
parser.add_argument('--run', type=int, default=7080)
parser.add_argument('--selection', type=str, default="*f10*")
parser.add_argument('--fout', type=str, default="")

args = parser.parse_args()
run = args.run
selection = args.selection
if not args.fout:
    fout_name = f"AngleDifference_Run{run}.root"
else:
    if args.fout.endswith(".root"): fout_name = args.fout
    else: fout_name = f"{args.fout}.root"

xy_eff_range: dict = {
    'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
    'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
}

input_dir = f"/eos/user/i/idioniso/1_Data/Tracks/run_{run:06d}"
mfout = "/eos/user/i/idioniso/mfout"
input_data = f"{input_dir}/{selection}.root"
print(f"INPUT: {input_data}")


cbmsim = TChain("cbmsim")
cbmsim.Add(f"{input_dir}/{selection}.root") 
nevents = cbmsim.GetEntries()

print(f"Events: {nevents:,}")

fout = TFile(f"{mfout}/{fout_name}", "recreate")

h_st_xz = TH1F(f"h_st_AngleDifferenceXZ", "Angle Difference between SciFi and DS (simple tracking);d#theta_{XZ} [mrad];", 1000, -100, 100)
h_st_yz = TH1F(f"h_st_AngleDifferenceYZ", "Angle Difference between SciFi and DS (simple tracking);d#theta_{YZ} [mrad];", 1000, -100, 100)
h_ht_xz = TH1F(f"h_ht_AngleDifferenceXZ", "Angle Difference between SciFi and DS (Hough transform);d#theta_{XZ} [mrad];", 1000, -100, 100)
h_ht_yz = TH1F(f"h_ht_AngleDifferenceYZ", "Angle Difference between SciFi and DS (Hough transform);d#theta_{YZ} [mrad];", 1000, -100, 100)



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



start_time = time()
count = 0
for i_event, event in enumerate(cbmsim):
    count = print_status(i_event, nevents, start_time, count)

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
            is_ip1(trk, event)
        ): continue
 
        if   tt==1:  tt1=True;  sf_st = trk
        elif tt==11: tt11=True; sf_ht = trk
        elif tt==3:  tt3=True;  ds_st = trk
        elif tt==13: tt13=True; ds_ht = trk
        else: continue

    if tt1 and tt3:
        h_st_xz.Fill(1e3*(ds_st.getAngleXZ()-sf_st.getAngleXZ()))
        h_st_yz.Fill(1e3*(ds_st.getAngleYZ()-sf_st.getAngleYZ()))
    if tt1 and tt3:
        h_ht_xz.Fill(1e3*(ds_ht.getAngleXZ()-sf_ht.getAngleXZ()))
        h_ht_yz.Fill(1e3*(ds_ht.getAngleYZ()-sf_ht.getAngleYZ()))

c = TCanvas(f"c_AngleDifference_Run{run}", "", 1100, 500)
c.Divide(2, 1)

c.cd(1)
gPad.SetLogy()
gPad.SetGridx()
h_ht_xz.SetLineColor(kRed)
h_ht_xz.SetLineWidth(2)
h_ht_xz.SetFillColorAlpha(kRed, 0.1)
h_ht_xz.Draw("hist")

h_st_xz.SetLineColor(kBlue)
h_st_xz.SetLineWidth(2)
h_st_xz.SetFillColorAlpha(kBlue, 0.1)
h_st_xz.Draw("hist same")

leg1 = TLegend(0.7, 0.7, 0.9, 0.9)
leg1.AddEntry(h_st_xz, "Simple tracking", "lf")
leg1.AddEntry(h_ht_xz, "Hough transform", "lf")
leg1.Draw("same")


c.cd(2)
gPad.SetLogy()
gPad.SetGridx()
h_ht_yz.SetLineColor(kRed)
h_ht_yz.SetLineWidth(2)
h_ht_yz.SetFillColorAlpha(kRed, 0.1)
h_ht_yz.Draw("hist")

h_st_yz.SetLineColor(kBlue)
h_st_yz.SetLineWidth(2)
h_st_yz.SetFillColorAlpha(kBlue, 0.1)
h_st_yz.Draw("hist same")
leg1.Draw("same")

c.Write()


fout.Write()
fout.Close()

