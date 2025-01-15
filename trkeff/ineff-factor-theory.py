import scipy.stats as stats
import matplotlib.pyplot as plt
from typing import Union
from time import time
from ROOT import TChain, sndRecoTrack, MuFilter, MuFilterHit, TClonesArray, Scifi, SNDLHCEventHeader, TF1, ShipMCTrack, TFile, TH1I, TMath, TVector3, TProfile, TProfile2D, TH2F, TH1F
import numpy as np
from ddf.pyfuncs import print_status, get_eff_with_error
from ddf.snd.mc import *
from ddf.snd.trk import xy_eff_range, get_intersection, is_good, is_within_ds3, is_within_us5_bar, is_within_veto_bar, sys, alg, sys_name, alg_name, get_anti_tt
from ddf.snd.trkeff import ref1_and_ref2_are_within_allowed_distance, ref1_is_within_eff_area
from selection import is_selected
from IneffF import *

input_dir = "/eos/user/i/idioniso/1_Data/Monte_Carlo"
A = {
    'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
    'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
}
z_ref = {1: 410., 11: 430., 3: 480., 13: 480.}

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
fout = TFile(f"{mfout}/ThetaRMS.root", "recreate")

h_ThetaRms = TH1F("h_ThetaRMS_sim", ";Entries;Outgoing angle root mean square $theta_{RMS} [mrad]", 30, 0, 10)

pr_xz_thetaRms = TProfile("pr_ThetaXZ.ThetaRMS-XZ_sim", ";#theta_{{XZ}} [mrad];#theta_{{XZ}} _{{RMS}} [mrad]", 20, -20, 20)
pr_yz_thetaRms = TProfile("pr_ThetaYZ.ThetaRMS-YZ_sim", ";#theta_{{YZ}} [mrad];#theta_{{YZ}} _{{RMS}} [mrad]", 20, -20, 20)
h_xz_thetaRms = TH2F("h_ThetaXZ.ThetaRMS-XZ_sim", ";#theta_{{XZ}} [mrad];#theta_{{XZ}} _{{RMS}} [mrad]", 20, -20, 20, 20, 1e-6, 1e-5)
h_yz_thetaRms = TH2F("h_ThetaYZ.ThetaRMS-YZ_sim", ";#theta_{{YZ}} [mrad];#theta_{{YZ}} _{{RMS}} [mrad]", 20, -20, 20, 20, 1e-6, 1e-5)


pr_xz_scat = {
    tt: TProfile(
        f"pr_xz.scat_{tt}_sim",
        f"{sys_name(tt)} {alg_name(tt)};#theta_{{XZ}} [mrad];Probability of misdirection",
        20, -20, 20
    ) for tt in (3, 13)
}
pr_yz_scat = {
    tt: TProfile(
        f"pr_yz.scat_{tt}_sim",
        f"{sys_name(tt)} {alg_name(tt)};#theta_{{YZ}} [mrad];Probability of misdirection",
        20, -20, 20
    ) for tt in (3, 13)
}
pr_xz_yz_scat = {
    tt: TProfile2D(
        f"pr_xz.yz.scat_{tt}_sim",
        f"{sys_name(tt)} {alg_name(tt)};#theta_{{XZ}} [mrad];#theta_{{YZ}} [mrad];Probability of misdirection",
        20, -20, 20, 20, -20, 20
    ) for tt in (3, 13)
}
h_xz_scat  = {
    tt: TH2F(
        f"h_xz.scat_{tt}_sim",
        f"{sys_name(tt)} {alg_name(tt)};#theta_{{XZ}} [mrad];Probability of misdirection",
        20, -20, 20, 20, 0, 1
    ) for tt in (3, 13)
}
h_yz_scat = {
    tt: TH2F(
        f"h_yz.scat_{tt}_sim", 
        f"{sys_name(tt)} {alg_name(tt)};#theta_{{YZ}} [mrad];Probability of misdirection",
        20, -20, 20, 20, 0, 1
    ) for tt in (3, 13)
}

pr_x_scat = {
    tt: TProfile(
        f"pr_x.scat_{tt}_sim",
        f"{sys_name(tt)} {alg_name(tt)};X_{{ref}} [cm];Probability of misdirection",
        20, -60, 0
    ) for tt in (3, 13)
}
pr_y_scat = {
    tt: TProfile(
        f"pr_y.scat_{tt}_sim",
        f"{sys_name(tt)} {alg_name(tt)};Y_{{ref}} [cm];Probability of misdirection",
        20, 0, 60
    ) for tt in (3, 13)
}
pr_x_y_scat = {
    tt: TProfile2D(
        f"pr_x.y.scat_{tt}_sim",
        f"{sys_name(tt)} {alg_name(tt)};X_{{ref}} [cm];Y_{{ref}} [cm];Probability of misdirection",
        20, -60, 0, 20, 0, 60
    ) for tt in (3, 13)
}
h_x_scat  = {
    tt: TH2F(
        f"h_x.scat_{tt}_sim",
        f"{sys_name(tt)} {alg_name(tt)};X_{{ref}} [cm];Probability of misdirection",
        20, -60, 0, 20, 0, 1
    ) for tt in (3, 13)
}
h_y_scat = {
    tt: TH2F(
        f"h_y.scat_{tt}_sim",
        f"{sys_name(tt)} {alg_name(tt)};Y_{{ref}} [cm];Probability of misdirection",
        20, 0, 60, 20, 0, 1
    ) for tt in (3, 13)
}

N_scat = {
    "v": {tt: 0 for tt in (3, 13)},
    "e": {tt: 0 for tt in (3, 13)}
}
N_total = {"v": 0, "e": 0}
saved = False
for mc in ("emd", "ni"):
    for i_event, event in enumerate(cbmsim[mc]):
        if not event.EventHeader.isIP1(): continue

        aTrkShouldExist = should_be_a_track(event)
        if not (aTrkShouldExist["sf"] and aTrkShouldExist["ds"]): continue

        w = event.MCTrack[0].GetWeight() * f[mc]
    
        sfp = get_sf_points(event)
        mc_point = sfp[-1]

        xz = get_xz(mc_point)
        yz = get_yz(mc_point)

        p = get_mom(mc_point)                         # [GeV/c]
        dz = abs(160/(TMath.Cos(xz) * TMath.Cos(yz))) # [cm]
        beta = get_beta(p)                     # []

        theta_rms = get_theta_rms(p, dz, X0, beta=beta)
        h_ThetaRms.Fill(1e3*theta_rms, w)
        #print(f"\n\nθᶦⁿ: {xz*1e3:.03f} [mrad]\nθʳᵐˢ: {theta_rms*1e3:.03f} [mrad]\np: {p:.03f} GeV\nbeta: {beta}\n\n")
        if p>24 and p<25 and not saved:
            Mom = p
            thetaXZ = -0.011652

            thetaRMS = theta_rms
            dThetaRMS = 0.11*thetaRMS            

            Xstart, Zstart = -40.25, 350.
            L = 130.
            dX_max = -11-Xstart                        # [cm]
            dX_min = -42-Xstart                        # [cm]

            thetaXZ_max = TMath.Tan(dX_max/L)          # [rad]
            thetaXZ_min = TMath.Tan(dX_min/L)          # [rad]

            if thetaXZ_max >  0.02: thetaXZ_max =  0.02
            if thetaXZ_min < -0.02: thetaXZ_min = -0.02

            saved = True

        pr_xz_thetaRms.Fill(1e3*xz, 1e3*theta_rms, w)
        pr_yz_thetaRms.Fill(1e3*yz, 1e3*theta_rms, w)
        h_xz_thetaRms.Fill(1e3*xz, 1e3*theta_rms, w)
        h_yz_thetaRms.Fill(1e3*yz, 1e3*theta_rms, w)


        #p_scat = {
        #    "v": {tt: get_p_in_Ax(mc_point, tt, A, z_ref, 0.11, theta_rms)[0] for tt in (3, 13)},
        #    "e": {tt: get_p_in_Ax(mc_point, tt, A, z_ref, 0.11, theta_rms)[1] for tt in (3, 13)}
        #}
        p_scat = {
            "v": {tt: get_p_scat(mc_point, tt, A, z_ref, 0.11)[0] for tt in (3, 13)},
            "e": {tt: get_p_scat(mc_point, tt, A, z_ref, 0.11)[1] for tt in (3, 13)}
        }
        for tt in (3, 13):
            N_scat["v"][tt] += w*(1-p_scat["v"][tt])
            N_scat["e"][tt] += w*p_scat["e"][tt] * w*p_scat["e"][tt]
        N_total["v"] += w
        N_total["e"] += w*w
        
        x = {}; y = {}
        for tt in (3, 13):
            x[tt] = get_intersection_z(mc_point, z_ref[tt]).X()
            y[tt] = get_intersection_z(mc_point, z_ref[tt]).Y()

            pr_xz_scat[tt].Fill(1e3*xz, 1-p_scat["v"][tt], w)
            pr_yz_scat[tt].Fill(1e3*yz, 1-p_scat["v"][tt], w)
            pr_xz_yz_scat[tt].Fill(1e3*xz, 1e3*yz, 1-p_scat["v"][tt], w)
            pr_x_scat[tt].Fill(x[tt], 1-p_scat["v"][tt], w)
            pr_y_scat[tt].Fill(y[tt], 1-p_scat["v"][tt], w)
            pr_x_y_scat[tt].Fill(x[tt], y[tt], 1-p_scat["v"][tt], w)
        
            h_xz_scat[tt].Fill(1e3*xz, 1-p_scat["v"][tt], w)
            h_yz_scat[tt].Fill(1e3*yz, 1-p_scat["v"][tt], w)
            h_x_scat[tt].Fill(x[tt], 1-p_scat["v"][tt], w)
            h_y_scat[tt].Fill(y[tt], 1-p_scat["v"][tt], w)

mu = 1e3*thetaXZ
sigma = 1e3*thetaRMS
sigma_eh = 1e3*(sigma + dThetaRMS/2)
sigma_el = 1e3*(sigma - dThetaRMS/2)

x = np.linspace(-21., 21., 500)
y = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
y_eh = (1 / (sigma_eh * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma_eh) ** 2)
y_el = (1 / (sigma_el * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma_el) ** 2)

plt.style.use("root")

plt.figure(figsize=(7, 3.5))
plt.plot(x, y, linestyle="-")
plt.fill_between(x, y_el, y_eh, alpha=0.2, label="11% uncertainty")
plt.vlines(1e3*thetaXZ_min, 0, 1.1*np.max(y), color='r', linestyle=":", linewidth=2)
plt.vlines(1e3*thetaXZ_max, 0, 1.1*np.max(y), color='r', linestyle=":", linewidth=2)
plt.text(0.4, 0.780, f"$\\theta_{{XZ\;in}}$\t$\\approx {mu:.03f}$ mrad", transform=plt.gcf().transFigure)
plt.text(0.4, 0.705, f"$\\theta_{{XZ\;RMS}}$\t$\\approx {sigma:.03f}$ mrad $\\pm$ 11%", transform=plt.gcf().transFigure)
plt.text(0.4, 0.630, f"$|\\vec{{p_{{\\mu}}}}|$\t$\\approx {Mom:.03f}$ GeV/c", transform=plt.gcf().transFigure)
plt.text(0.4, 0.560, f"US entrance point $X\\approx {Xstart:.01f}$ cm", transform=plt.gcf().transFigure)
plt.text(0.4, 0.490, f"$X_{{ref}} \\in$ [$-42$ cm, $-11$ cm]", transform=plt.gcf().transFigure)
plt.xlabel("Angle of outgoing muon $\\theta_{XZ\;out}$ [mrad]")
plt.ylabel('$pdf(\\theta_{XZ\;out})$')
plt.xlim(-21, 21)
plt.ylim(0, np.max(y)*1.1)

plt.tight_layout()
plt.savefig(f'{mfout}/outXZprob.png')
plt.close()

fout.Write()
fout.Close()

# N_ang - Number of DS tracks that acquire angles > 20 mrad
# N_A   - Number of DS tracks that acquire angles < 20 mrad, but are scattered out of fiducial area
# k = N_A / N_ang
#   Calculated by from data in `/afs/cern.ch/work/i/idioniso/0_Workdir/muon_flux/prpts/prpts_noDStrks.py`
k  = {3: 0.21080847463861405,   13: 0.22846914623007492 }
dk = {3: 0.0010833184054139875, 13: 0.001173995668478104}

N_total["e"] = np.sqrt(N_total["e"])
for tt in (3, 13):
    N_scat["e"][tt] = np.sqrt(N_scat["e"][tt])


alpha = {3: {}, 13: {}}
for tt in (3, 13):
    alpha[tt]['v'], alpha[tt]['eh'], alpha[tt]['el'] = get_eff_with_error(N_scat["v"][tt], N_total["v"], stat_option="bayesian")

eps_ds = {
    tt: {
        'v':  alpha[tt]['v'],
        'el': alpha[tt]['el'],
        'eh': alpha[tt]['eh']
    } for tt in (3, 13)
}

print("\n\n"+40*"-")
print(f"Nᵖᵃˢˢᵉᵈ:\t{N_scat['v'][3]} ± {N_scat['e'][3]}")
print(f"Nᵗᵒᵗᵃˡ: \t{N_total['v']} ± {N_total['e']}")
print(f"Simple tracking:\tϵᵈˢ = {eps_ds[3]['v']} ± (+{eps_ds[3]['eh']}, -{eps_ds[3]['el']})")
print(f"Hough transform:\tϵᵈˢ = {eps_ds[13]['v']} ± (+{eps_ds[13]['eh']}, -{eps_ds[13]['el']})")
print(40*"-"+"\n")
