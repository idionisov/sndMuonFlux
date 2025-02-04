import scipy.stats as stats
import matplotlib.pyplot as plt
import ROOT
from typing import Union
from time import time
import numpy as np
from ddfUtils import printStatus, getEffWithError
from sndUtils import SndMCData, DdfTrack, DdfMCTrack, system, algorithm, sys, alg, \
    sfTrackIsReconstructible, dsTrackIsReconstructible, thereIsAMuon
from helpers import *
from hists import *

inputDir = "/eos/user/i/idioniso/1_Data/Monte_Carlo"
mfout = "/eos/user/i/idioniso/mfout"
A = {
    'sf': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}},
    'ds': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}}
}
zRef = {1: 430., 11: 450., 3: 430., 13: 450.}
f = {
    'EMD': 0.1388888888888889,
    'NI':  2.003205128205128
}


def main():
    data = {
        a: SndMCData(
            InputDir=inputDir,
            Files=f"muonReco_MC-{a}_PbPb.root",
            Geofile=f"{inputDir}/geofile_MC.PbPb.root"
        ) for a in ("EMD", "NI")
    }
    nEntries = data["EMD"].Tree.GetEntries() + data["NI"].Tree.GetEntries()
    for a in ("EMD", "NI"):
        data[a].Print()


    fout = ROOT.TFile(f"{mfout}/ThetaRMS.root", "recreate")

    start_time = time()
    count = 0
    N_scat = {
        "v": {tt: 0 for tt in (3, 13)},
        "e": {tt: 0 for tt in (3, 13)}
    }
    N_total = {"v": 0, "e": 0}

    events = 0
    saved = False
    thetaXZ_min=0.02
    thetaXZ_max=0.02
    for mc in data:

        if mc == list(data.keys())[0]:
            events = 0
        else:
            events += data[mc].Tree.GetEntries()

        for i_event, event in enumerate(data[mc].Tree):
            count = printStatus(i_event+events, nEntries, start_time, count)

            if not (event.EventHeader.isIP1() and thereIsAMuon(event)):
                continue

            _sf = sfTrackIsReconstructible(event)
            _ds = dsTrackIsReconstructible(event)

            if _sf==False and _ds==False:
                continue

            reco = {1: _sf, 11: _sf, 3: _ds, 13: _ds}
            w = f[mc] * event.MCTrack[0].GetWeight()

            sfp = getPointsSF(event)
            if not sfp:
                continue
            mcPoint = sfp[-1]

            xz = getXZ(mcPoint)
            yz = getYZ(mcPoint)

            p = getMom(mcPoint)                         # [GeV/c]
            dz = abs(160/(TMath.Cos(xz) * TMath.Cos(yz))) # [cm]
            beta = getBeta(p)                     # []

            theta_rms = getThetaRMS(p, dz, X0, beta=beta)
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
            #    "v": {tt: get_p_in_Ax(mcPoint, tt, A, zRef, 0.11, theta_rms)[0] for tt in (3, 13)},
            #    "e": {tt: get_p_in_Ax(mcPoint, tt, A, zRef, 0.11, theta_rms)[1] for tt in (3, 13)}
            #}
            p_scat = {
                "v": {tt: getProbScatter(mcPoint, tt, A, zRef, 0.11)[0] for tt in (3, 13)},
                "e": {tt: getProbScatter(mcPoint, tt, A, zRef, 0.11)[1] for tt in (3, 13)}
            }
            for tt in (3, 13):
                N_scat["v"][tt] += w*(1-p_scat["v"][tt])
                N_scat["e"][tt] += w*p_scat["e"][tt] * w*p_scat["e"][tt]
            N_total["v"] += w
            N_total["e"] += w*w

            x = {}; y = {}
            for tt in (3, 13):
                x[tt] = getPointAtZ(mcPoint, zRef[tt]).X()
                y[tt] = getPointAtZ(mcPoint, zRef[tt]).Y()

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
    plt.text(0.4, 0.490, f"$X_{{ref}} \\in$ [$-42$ cm, $-10$ cm]", transform=plt.gcf().transFigure)
    plt.xlabel("Angle of outgoing muon $\\theta_{XZ\;out}$ [mrad]")
    plt.ylabel('$pdf(\\theta_{XZ\;out})$')
    plt.xlim(-21, 21)
    plt.ylim(0, np.max(y)*1.1)

    plt.tight_layout()
    plt.savefig(f'{mfout}/outXZprob.png')
    plt.close()


    h_ThetaRms.Write()
    pr_xz_thetaRms.Write()
    pr_yz_thetaRms.Write()
    h_xz_thetaRms.Write()
    h_yz_thetaRms.Write()
    for tt in (3, 13):
        pr_xz_scat[tt].Write()
        pr_xz_scat[tt].Write()
        pr_yz_scat[tt].Write()
        pr_xz_yz_scat[tt].Write()
        h_xz_scat[tt].Write()
        h_yz_scat[tt].Write()
        pr_x_scat[tt].Write()
        pr_y_scat[tt].Write()
        pr_x_y_scat[tt].Write()
        h_x_scat[tt].Write()
        h_y_scat[tt].Write()
    fout.Close()

    # N_ang - Number of DS tracks that acquire angles > 20 mrad
    # N_A   - Number of DS tracks that acquire angles < 20 mrad, but are scattered out of fiducial area
    # k = N_A / N_ang
    k  = {3: 0.21080847463861405,   13: 0.22846914623007492 }
    dk = {3: 0.0010833184054139875, 13: 0.001173995668478104}

    N_total["e"] = np.sqrt(N_total["e"])
    for tt in (3, 13):
        N_scat["e"][tt] = np.sqrt(N_scat["e"][tt])

    alpha = {3: {}, 13: {}}
    for tt in (3, 13):
        alpha[tt]['v'], alpha[tt]['eh'], alpha[tt]['el'] = getEffWithError(N_scat["v"][tt], N_total["v"], statOption="bayesian")

    eps_ds = {
        tt: {
            'v':  alpha[tt]['v'],
            'el': alpha[tt]['el'],
            'eh': alpha[tt]['eh']
        } for tt in (3, 13)
    }

    print("\n\n"+40*"-")
    print(f"Nᵖᵃˢˢᵉᵈ:         {N_scat['v'][3]} ± {N_scat['e'][3]}")
    print(f"Nᵗᵒᵗᵃˡ:          {N_total['v']} ± {N_total['e']}")
    print(f"Simple tracking: fᵈˢ [%] = {100*(1-eps_ds[3]['v']):.02f} (+{100*eps_ds[3]['el']:.02f}, -{100*eps_ds[3]['eh']:.02f})")
    print(f"Hough transform: fᵈˢ [%] = {100*(1-eps_ds[13]['v']):.02f} (+{100*eps_ds[13]['el']:.02f}, -{100*eps_ds[13]['eh']:.02f})")
    print(40*"-"+"\n")



if __name__=="__main__":
    main()
