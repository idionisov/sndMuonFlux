from typing import Union
import numpy as np
from ROOT import TH1F, TH2F
import roostyling
from sndUtils import sys, system, alg, algorithm, nName, att

ang_bin_edges_sf = np.array([-90, -70, -50, -30, -22.5, -15, -10, -7.5, -5, -2.5, 2.5, 5, 7.5, 10, 15, 22.5, 30, 50, 70, 90], dtype=np.float64)
ang_bin_edges_ds = np.array([-90, -50, -30, -12.5, -5, 5, 12.5, 30, 50, 90], dtype=np.float64)


def get_h_yz_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if sys(tt)=="sf":
        bin_edges = ang_bin_edges_sf
    elif sys(tt)=="ds":
        bin_edges = ang_bin_edges_sf
    else: raise ValueError(f"{tt} is an invalid track type!")

    h_yz_eff = tuple(
        TH1F(
            f"h{i}_yz_{tt}_{run_or_mcSet}",
            f"{i} {tt} {run_or_mcSet};#theta_{{YZ}} [mrad];Tracking Efficiency",
            len(bin_edges)-1, bin_edges
        ) for i in range(1, 3)
    )
    roostyling.axes(h_yz_eff[0], h_yz_eff[1])
    return h_yz_eff


def get_h_yz_x_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if sys(tt)=="sf":
        bin_edges = ang_bin_edges_sf
    elif sys(tt)=="ds":
        bin_edges = ang_bin_edges_sf
    else: raise ValueError(f"{tt} is an invalid track type!")

    h_yz_x_eff = tuple(
        TH2F(
            f"h{i}_yz.x_{tt}_{run_or_mcSet}",
            f"{i} {tt} {run_or_mcSet};#theta_{{YZ}} [mrad];X [cm];Tracking Efficiency",
            len(bin_edges)-1, bin_edges, 40, -75, 5
        ) for i in range(1, 3)
    )
    roostyling.axes(h_yz_x_eff[0], h_yz_x_eff[1])
    return h_yz_x_eff


def get_h_yz_y_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if sys(tt)=="sf":
        bin_edges = ang_bin_edges_sf
    elif sys(tt)=="ds":
        bin_edges = ang_bin_edges_sf
    else: raise ValueError(f"{tt} is an invalid track type!")

    h_yz_y_eff = tuple(
        TH2F(
            f"h{i}_yz.y_{tt}_{run_or_mcSet}",
            f"{i} {tt} {run_or_mcSet};#theta_{{YZ}} [mrad];Y [cm];Tracking Efficiency",
            len(bin_edges)-1, bin_edges, 40, 5, 75
        ) for i in range(1, 3)
    )
    roostyling.axes(h_yz_y_eff[0], h_yz_y_eff[1])
    return h_yz_y_eff


def get_h_yz_chi2_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if tt==1:
        chi2_max = 60.; nbins_chi2 = 30
        bin_edges = ang_bin_edges_sf
    elif tt==11:
        chi2_max = 90.; nbins_chi2 = 45
        bin_edges = ang_bin_edges_sf
    elif tt==3:
        chi2_max = 1000.; nbins_chi2 = 200
        bin_edges = ang_bin_edges_sf
    elif tt==13:
        chi2_max = 1200.; nbins_chi2 = 300
        bin_edges = ang_bin_edges_sf
    else: raise ValueError(f"{tt} is an invalid track type!")


    h_yz_chi2_eff = tuple(
        TH2F(
            f"h{i}_yz.chi2_{tt}_{run_or_mcSet}",
            f"{system(tt)} {algorithm(tt)};#theta_{{YZ}} [mrad];#chi2;Tracking Efficiency",
            len(bin_edges)-1, bin_edges, nbins_chi2, 0, chi2_max
        ) for i in range(1, 3)
    )

    roostyling.axes(h_yz_chi2_eff[0], h_yz_chi2_eff[1])
    return h_yz_chi2_eff


def get_h_yz_chi2ndf_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if tt==1:
        chi2ndf_max = 60.; chi2ndf_nbins = 30
        bin_edges   = ang_bin_edges_sf
    elif tt==11:
        chi2ndf_max = 35.; chi2ndf_nbins = 35
        bin_edges   = ang_bin_edges_sf
    elif tt==3:
        chi2ndf_max = 200.; chi2ndf_nbins = 100
        bin_edges   = ang_bin_edges_sf
    elif tt==13:
        chi2ndf_max = 150.; chi2ndf_nbins = 50
        bin_edges   = ang_bin_edges_sf
    else:
        raise ValueError(f"{tt} is an invalid track type!")



    h_yz_chi2ndf_eff = tuple(
        TH2F(
            f"h{i}_yz.chi2ndf_{tt}_{run_or_mcSet}",
            f"{system(tt)} {algorithm(tt)};#theta_{{YZ}} [mrad];#chi2/ndf;Tracking Efficiency",
            len(bin_edges)-1, bin_edges, chi2ndf_nbins, 0, chi2ndf_max
        ) for i in range(1, 3)
    )

    roostyling.axes(h_yz_chi2ndf_eff[0], h_yz_chi2ndf_eff[1])
    return h_yz_chi2ndf_eff



def get_h_yz_n_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if tt==1:
        n_max=55; nbins_n=50
        bin_edges = ang_bin_edges_sf
    elif tt==11:
        n_max=125; nbins_n=120
        bin_edges = ang_bin_edges_sf
    elif tt==3:
        n_max=20; nbins_n=15
        bin_edges = ang_bin_edges_sf
    elif tt==13:
        n_max=60; nbins_n=55
        bin_edges = ang_bin_edges_sf
    else: raise ValueError(f"{tt} is an invalid track type!")


    h_yz_n_eff = tuple(
        TH2F(
            f"h{i}_yz.n_{tt}_{run_or_mcSet}",
            f"{system(tt)} {algorithm(tt)};#theta_{{YZ}} [mrad];Number of {nName(tt)};Tracking Efficiency",
            len(bin_edges)-1, bin_edges, nbins_n, 5, n_max
        ) for i in range(1, 3)
    )

    roostyling.axes(h_yz_n_eff[0], h_yz_n_eff[1])
    return h_yz_n_eff


def get_h_yz_trkP_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if tt==1:
        trkP_min=6; trkP_max=10; nbins_trkP=4
        bin_edges = ang_bin_edges_sf
    elif tt==11:
        trkP_min=6; trkP_max=21; nbins_trkP=15
        bin_edges = ang_bin_edges_sf
    elif tt==3:
        trkP_min=6; trkP_max=26; nbins_trkP=20
        bin_edges = ang_bin_edges_sf
    elif tt==13:
        trkP_min=6; trkP_max=41; nbins_trkP=35
        bin_edges = ang_bin_edges_sf
    else: raise ValueError(f"{tt} is an invalid track type!")

    _att = att(tt)

    h_yz_trkP_eff = tuple(
        TH2F(
            f"h{i}_yz.trkP_{tt}_{run_or_mcSet}",
            f"{system(tt)} {algorithm(tt)};#theta_{{YZ}} [mrad];Number of tag track's points (utilized {nName(_att)});Tracking Efficiency",
            len(bin_edges)-1, bin_edges, nbins_trkP, trkP_min, trkP_max
        ) for i in range(1, 3)
    )

    roostyling.axes(h_yz_trkP_eff[0], h_yz_trkP_eff[1])
    return h_yz_trkP_eff







def get_h_dyz(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if tt not in (1, 11, 3, 13): raise ValueError(f"{tt} is an invalid track type!")

    h_dyz = TH1F(
        f"h_dYZ_{tt}_{run_or_mcSet}",
        f"{system(tt)} {algorithm(tt)} ({run_or_mcSet});YZ angle difference [mrad];",
        50, -50, 50
    )
    roostyling.axes(h_dyz)
    return h_dyz
