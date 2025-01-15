from typing import Union
from time import time
from ROOT import TChain, sndRecoTrack, MuFilter, MuFilterHit, TClonesArray, Scifi, SNDLHCEventHeader, TF1, ShipMCTrack, TFile, TH1I, TH1F, TH2F, TMath, TVector3
import numpy as np
from ddf.pyfuncs import print_status, get_eff_with_error
from ddf.snd.mc import *
from ddf.snd.trk import xy_eff_range, get_intersection, is_good, is_within_ds3, is_within_us5_bar, is_within_veto_bar, sys, alg, sys_name, alg_name, get_anti_tt
from ddf.snd.trkeff import ref1_and_ref2_are_within_allowed_distance, ref1_is_within_eff_area
from IneffF import *

input_dir = "/eos/user/i/idioniso/1_Data/Monte_Carlo"
A = {
    'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
    'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
}
z_ref = {1: 440., 11: 435., 3: 460., 13: 475.}

cbmsim = {"emd": TChain("cbmsim"), "ni": TChain("cbmsim")}
cbmsim["emd"].Add(f"{input_dir}/muonReco_MC-EMD_PbPb.root")
cbmsim["ni"].Add(f"{input_dir}/muonReco_MC-NI_PbPb.root")

nentries = cbmsim["emd"].GetEntries() + cbmsim["ni"].GetEntries()
print(f"Entries:\t{nentries:,}")

f = {
    'emd': 0.1388888888888889,
    'ni':  2.003205128205128
}

mfout = "/eos/user/i/idioniso/mfout"
fout = TFile(f"{mfout}/ineff.ds-sim.root", "recreate")

h = {
    tt: {
        "x":     TH1F(f"h_x_ineff_{tt}_sim",   f"Scattered {sys(tt).upper()} {alg(tt).upper()} tracks X distribution;X [cm]",  40, -70, 10),
        "y":     TH1F(f"h_y_ineff_{tt}_sim",   f"Scattered {sys(tt).upper()} {alg(tt).upper()} tracks Y distribution;Y [cm]",  40, -10, 70),
        "x.y":   TH2F(f"h_x.y_ineff_{tt}_sim", f"Scattered {sys(tt).upper()} {alg(tt).upper()} tracks XY distribution;X [cm]; Y [cm]", 40, -70, 10, 40, -10, 70),
        "xz": TH1F(f"h_xz_ineff_{tt}_sim",  f"Scattered {sys(tt).upper()} {alg(tt).upper()} XZ angle distribution;#theta_{{XZ}} [mrad]", 100, -100, 100),
        "yz": TH1F(f"h_yz_ineff_{tt}_sim",  f"Scattered {sys(tt).upper()} {alg(tt).upper()} YZ angle distribution;#theta_{{YZ}} [mrad]", 100, -100, 100),
        "xz.yz": TH2F(f"h_xz.yz_ineff_{tt}_sim",  f"Scattered {sys(tt).upper()} {alg(tt).upper()} angle distribution;#theta_{{XZ}} [mrad];#theta_{{YZ}} [mrad]", 100, -100, 100, 100, -100, 100)
    } for tt in (3, 13)
}


for mc in ("emd", "ni"):
    for i_event, event in enumerate(cbmsim[mc]):
        if not event.EventHeader.isIP1(): continue

        aTrkShouldExist = should_be_a_track(event)
        if not (aTrkShouldExist["sf"] and aTrkShouldExist["ds"]): continue

        dsP = get_ds_points(event)[0]
        x  = {tt: get_intersection_z(dsP, z_ref[tt]).X() for tt in (3, 13)}
        y  = {tt: get_intersection_z(dsP, z_ref[tt]).Y() for tt in (3, 13)}
        xz = 1e3*get_xz(dsP)
        yz = 1e3*get_yz(dsP)

        w = event.MCTrack[0].GetWeight() * f[mc]
    
        us_scat = us_scatter(event, A, z_ref)
        for tt in (3, 13):
            if not us_scat[tt]: continue
            h[tt]["x"].Fill(x[tt], w)
            h[tt]["y"].Fill(y[tt], w)
            h[tt]["x.y"].Fill(x[tt], y[tt], w)
        
            h[tt]["xz"].Fill(xz, w)
            h[tt]["yz"].Fill(xz, w)
            h[tt]["xz.yz"].Fill(xz, yz, w)
             
fout.Write()
fout.Close()
