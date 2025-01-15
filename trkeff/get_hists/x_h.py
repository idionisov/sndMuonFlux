from typing import Union
import numpy as np
from ROOT import TH1F, TH2F
import roostyling
from sndUtils import sys, system, alg, algorithm, nName, att

xy_eff_range = {
    'min': {'x': -42., 'y': 19.},
    'max': {'x': -10., 'y': 48.}
}
xy_full_range = {
    'min': {'x': -70, 'y': -5},
    'max': {'x':  10, 'y':  75}
}
ang_bin_edges_sf = np.array([-90, -70, -50, -30, -22.5, -15, -10, -7.5, -5, -2.5, 2.5, 5, 7.5, 10, 15, 22.5, 30, 50, 70, 90], dtype=np.float64)
ang_bin_edges_ds = np.array([-90, -50, -30, -12.5, -5, 5, 12.5, 30, 50, 90], dtype=np.float64)



def get_h_x_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if tt not in (1, 11, 3, 13): raise ValueError(f"{tt} is an invalid track type!")

    h_x_eff = tuple(
        TH1F(
            f"h{i}_x_{tt}_{run_or_mcSet}",
            f"{i} {tt} {run_or_mcSet};X [cm];Tracking Efficiency",
            10, xy_eff_range["min"]["x"], xy_eff_range["max"]["x"]
        ) for i in range(1, 3)
    )
    roostyling.axes(h_x_eff[0], h_x_eff[1])
    return h_x_eff


def get_h_x_y_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if tt not in (1, 11, 3, 13): raise ValueError(f"{tt} is an invalid track type!")

    h_x_y_eff = tuple(
        TH2F(
            f"h{i}_x.y_{tt}_{run_or_mcSet}",
            f"{i} {tt} {run_or_mcSet};X [cm];Y [cm]",
            40, xy_full_range["min"]["x"], xy_full_range["max"]["x"],
            40, xy_full_range["min"]["y"], xy_full_range["max"]["y"]
        ) for i in range(1, 3)
    )
    roostyling.axes(h_x_y_eff[0], h_x_y_eff[1])
    return h_x_y_eff

def get_h_x_xz_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if sys(tt)=="sf":
        bin_edges = ang_bin_edges_sf
    elif sys(tt)=="ds":
        bin_edges = ang_bin_edges_ds
    else: raise ValueError(f"{tt} is an invalid track type!")

    h_x_xz_eff = tuple(
        TH2F(
            f"h{i}_x.xz_{tt}_{run_or_mcSet}",
            f"{system(tt)} {algorithm(tt)};X [cm];#theta_{{XZ}} [mrad];Tracking Efficiency",
            10, xy_eff_range["min"]["x"], xy_eff_range["max"]["x"],
            len(bin_edges)-1, bin_edges
        ) for i in range(1, 3)
    )

    roostyling.axes(h_x_xz_eff[0], h_x_xz_eff[1])

    return h_x_xz_eff

def get_h_x_yz_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if sys(tt)=="sf":
        bin_edges = ang_bin_edges_sf
    elif sys(tt)=="ds":
        bin_edges = ang_bin_edges_ds
    else: raise ValueError(f"{tt} is an invalid track type!")

    h_x_yz_eff = tuple(
        TH2F(
            f"h{i}_x.yz_{tt}_{run_or_mcSet}",
            f"{system(tt)} {algorithm(tt)};X [cm];#theta_{{YZ}} [mrad];Tracking Efficiency",
            10, xy_eff_range["min"]["x"], xy_eff_range["max"]["x"],
            len(bin_edges)-1, bin_edges
        ) for i in range(1, 3)
    )

    roostyling.axes(h_x_yz_eff[0], h_x_yz_eff[1])
    return h_x_yz_eff

def get_h_x_n_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if   tt==1:  n_max=55;  nbins=50
    elif tt==11: n_max=125; nbins=120
    elif tt==3:  n_max=20;  nbins=15
    elif tt==13: n_max=60;  nbins=55
    else: raise ValueError(f"{tt} is an invalid track type!")


    h_x_n_eff = tuple(
        TH2F(
            f"h{i}_x.n_{tt}_{run_or_mcSet}",
            f"{system(tt)} {algorithm(tt)};X [cm];Number of {nName(tt)};Tracking Efficiency",
            10, xy_eff_range["min"]["x"], xy_eff_range["max"]["x"], nbins, 5, n_max
        ) for i in range(1, 3)
    )

    roostyling.axes(h_x_n_eff[0], h_x_n_eff[1])
    return h_x_n_eff



def get_h_x_chi2_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if   tt==1:  chi2_max = 60.;   nbins_chi2 = 30
    elif tt==11: chi2_max = 90.;   nbins_chi2 = 45
    elif tt==3:  chi2_max = 1000.; nbins_chi2 = 200
    elif tt==13: chi2_max = 1200.; nbins_chi2 = 300
    else: raise ValueError(f"{tt} is an invalid track type!")



    h_x_chi2_eff = tuple(
        TH2F(
            f"h{i}_x.chi2_{tt}_{run_or_mcSet}",
            f"{system(tt)} {algorithm(tt)};X [cm];#chi2;Tracking Efficiency",
            10, xy_eff_range["min"]["x"], xy_eff_range["max"]["x"],
            nbins_chi2, 0, chi2_max
        ) for i in range(1, 3)
    )

    roostyling.axes(h_x_chi2_eff[0], h_x_chi2_eff[1])
    return h_x_chi2_eff



def get_h_x_chi2ndf_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if   tt==1:  chi2ndf_max = 60.;  nbins_chi2ndf = 30
    elif tt==11: chi2ndf_max = 35.;  nbins_chi2ndf = 35
    elif tt==3:  chi2ndf_max = 200.; nbins_chi2ndf = 100
    elif tt==13: chi2ndf_max = 150.; nbins_chi2ndf = 50
    else: raise ValueError(f"{tt} is an invalid track type!")


    h_x_chi2ndf_eff = tuple(
        TH2F(
            f"h{i}_x.chi2ndf_{tt}_{run_or_mcSet}",
            f"{system(tt)} {algorithm(tt)};X [cm];#chi2/ndf;Tracking Efficiency",
            10, xy_eff_range["min"]["x"], xy_eff_range["max"]["x"],
            nbins_chi2ndf, 0, chi2ndf_max
        ) for i in range(1, 3)
    )

    roostyling.axes(h_x_chi2ndf_eff[0], h_x_chi2ndf_eff[1])
    return h_x_chi2ndf_eff




def get_h_x_trkP_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if   tt==1:  trkP_min = 6; trkP_max = 10; trkP_nbins = 4
    elif tt==11: trkP_min = 6; trkP_max = 21; trkP_nbins = 15
    elif tt==3:  trkP_min = 6; trkP_max = 26; trkP_nbins = 20
    elif tt==13: trkP_min = 6; trkP_max = 41; trkP_nbins = 35
    else: raise valueerror(f"{tt} is an invalid track type!")

    _att = att(tt)
    h_x_trkP_eff = tuple(
        TH2F(
            f"h{i}_x.trkP_{tt}_{run_or_mcSet}",
            f"{system(tt)} {algorithm(tt)};X [cm];Number of tag rack's points (utilized {nName(_att)});Tracking Efficiency",
            10, xy_eff_range["min"]["x"], xy_eff_range["max"]["x"],
            trkP_nbins, trkP_min, trkP_max
        ) for i in range(1, 3)
    )

    roostyling.axes(h_x_trkP_eff[0], h_x_trkP_eff[1])
    return h_x_trkP_eff



def get_h_dxRef(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if tt not in (1, 11, 3, 13): raise ValueError(f"{tt} is an invalid track type!")

    h_dxRef = TH1F(
        f"h_dxRef_{tt}_{run_or_mcSet}",
        f"{system(tt)} {algorithm(tt)};X-Xtag [cm];",
        100, -10, 10
    )
    roostyling.axes(h_dxRef)
    return h_dxRef
