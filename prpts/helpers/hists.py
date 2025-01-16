import numpy as np
from typing import Union
from ROOT import TH1F, TH2F, TH1I
from sndUtils import system, algorithm, nName
import roostyling

def get_h_x(
    run_or_mcSet: Union[int, str] = 7080,
    tt:    int   =  1,
    xmin:  float = -70.,
    xmax:  float =  10.,
    xBins: int   =  40
):
    h_x = TH1F(
        f"h_x_{tt}_{run_or_mcSet}",
        f"{system(tt)} {algorithm(tt)} ({run_or_mcSet});X [cm];",
        xBins, xmin, xmax
    )
    roostyling.axes(h_x)
    return h_x

def get_h_y(
    run_or_mcSet: Union[int, str] = 7080,
    tt:    int   =  1,
    ymin:  float = -5.,
    ymax:  float =  75.,
    yBins: int   =  80
):
    h_y = TH1F(
        f"h_y_{tt}_{run_or_mcSet}",
        f"{system(tt)} {algorithm(tt)} ({run_or_mcSet});Y [cm];",
        yBins, ymin, ymax
    )
    roostyling.axes(h_y)
    return h_y


def get_h_xy(
    run_or_mcSet: Union[int, str] = 7080,
    tt:    int   =  1,
    xmin:  float = -70.,
    xmax:  float =  10.,
    ymin:  float = -5.,
    ymax:  float =  75.,
    xBins: int   =  40,
    yBins: int   =  40
):
    h_xy = TH2F(
        f"h_xy_{tt}_{run_or_mcSet}",
        f"{system(tt)} {algorithm(tt)} ({run_or_mcSet});X [cm];y [cm];",
        xBins, xmin, xmax, yBins, ymin, ymax
    )
    roostyling.axes(h_xy)
    return h_xy

def get_h_xz(
    run_or_mcSet: Union[int, str] = 7080,
    tt: int = 1
):
    center = 4.2
    sfBinWidth = 2.5
    dsBinWidth = 20

    xz_min = -300. + center
    xz_max =  300. + center

    if tt==1 or tt==11:
        xz_bin_edges = np.arange(xz_min, xz_max, sfBinWidth)
    elif tt==3 or tt==13:
        xz_bin_edges = np.arange(xz_min, xz_max, dsBinWidth)
    else: raise ValueError(f"Invalid track type: {tt}")

    h_xz = TH1F(
        f"h_xz_{tt}_{run_or_mcSet}",
        f"{system(tt)} {algorithm(tt)} ({run_or_mcSet});#theta_{{XZ}} [mrad];",
        len(xz_bin_edges)-1, xz_bin_edges
    )
    roostyling.axes(h_xz)
    return h_xz


def get_h_yz(
    run_or_mcSet: Union[int, str] = 7080,
    tt: int = 1
):
    center = 3.6
    sfBinWidth = 2.5
    dsBinWidth = 20

    yz_min = -300. + center
    yz_max =  300. + center

    if tt==1 or tt==11:
        yz_bin_edges = np.arange(yz_min, yz_max, sfBinWidth)
    elif tt==3 or tt==13:
        yz_bin_edges = np.arange(yz_min, yz_max, dsBinWidth)
    else: raise ValueError(f"Invalid track type: {tt}")

    h_yz = TH1F(
        f"h_yz_{tt}_{run_or_mcSet}",
        f"{system(tt)} {algorithm(tt)} ({run_or_mcSet});#theta_{{YZ}} [mrad];",
        len(yz_bin_edges)-1, yz_bin_edges
    )
    roostyling.axes(h_yz)
    return h_yz


def get_h_xzyz(
    run_or_mcSet: Union[int, str] = 7080,
    tt: int = 1
):
    xz_min = -120.
    xz_max =  120.
    yz_min = -120.
    yz_max =  120.

    if tt==1 or tt==11:
        xz_bins=60; yz_bins=60
    elif tt==3 or tt==13:
        xz_bins=30; yz_bins=30
    else: raise ValueError(f"Invalid track type: {tt}!")

    h_xzyz = TH2F(
        f"h_xz.yz_{tt}_{run_or_mcSet}",
        f"{system(tt)} {algorithm(tt)} ({run_or_mcSet});#theta_{{XZ}} [mrad];#theta_{{YZ}} [mrad];",
        xz_bins, xz_min, xz_max, yz_bins, yz_min, yz_max
    )
    roostyling.axes(h_xzyz)
    return h_xzyz


def get_h_n(
    run_or_mcSet: Union[int, str] = 7080,
    tt: int = 1
):
    if   tt==1:  nMax=55
    elif tt==11: nMax=125
    elif tt==3:  nMax=20
    elif tt==13: nMax=60
    else: raise ValueError(f"Invalid track type: {tt}")

    nMin = 5
    nBins = int(nMax-nMin)

    h_n = TH1I(
        f"h_n_{tt}_{run_or_mcSet}",
        f"{system(tt)} {algorithm(tt)} ({run_or_mcSet});Number of {nName(tt)};",
        nBins, nMin, nMax
    )
    roostyling.axes(h_n)
    return h_n


def get_h_trkP(
    run_or_mcSet: Union[int, str] = 7080,
    tt: int = 1
):
    if   tt==1:  trkP_max=18
    elif tt==11: trkP_max=40
    elif tt==3:  trkP_max=8
    elif tt==13: trkP_max=16
    else: raise ValueError(f"Invalid track type: {tt}")

    trkP_min=6
    trkP_bins=int(trkP_max-trkP_min)

    h_trkP = TH1I(
        f"h_trkP_{tt}_{run_or_mcSet}",
        f"{system(tt)} {algorithm(tt)} ({run_or_mcSet});Number of utilized {nName(tt)};",
        trkP_bins, trkP_min, trkP_max
    )
    roostyling.axes(h_trkP)
    return h_trkP


def get_h_chi2(
    run_or_mcSet: Union[int, str] = 7080,
    tt: int = 1
):
    if   tt==1:  chi2_min=0.; chi2_max=700.;  chi2_bins=175
    elif tt==11: chi2_min=0.; chi2_max=1000.; chi2_bins=250
    elif tt==3:  chi2_min=0.; chi2_max=35.;   chi2_bins=140
    elif tt==13: chi2_min=0.; chi2_max=70.;   chi2_bins=140
    else: raise ValueError(f"Invalid track type: {tt}")

    h_chi2 = TH1F(
        f"h_chi2_{tt}_{run_or_mcSet}",
        f"{system(tt)} {algorithm(tt)} ({run_or_mcSet});#chi^{{2}};",
        chi2_bins, chi2_min, chi2_max
    )
    roostyling.axes(h_chi2)
    return h_chi2

def get_h_chi2ndf(
    run_or_mcSet: Union[int, str] = 7080,
    tt: int = 1
):
    if   tt==1:  chi2ndf_min=0.; chi2ndf_max=125.; chi2ndf_bins=125
    elif tt==11: chi2ndf_min=0.; chi2ndf_max=250.; chi2ndf_bins=125
    elif tt==3:  chi2ndf_min=0.; chi2ndf_max=30.;  chi2ndf_bins=120
    elif tt==13: chi2ndf_min=0.; chi2ndf_max=20.;  chi2ndf_bins=120
    else: raise ValueError(f"Invalid track type: {tt}")

    h_chi2ndf = TH1F(
        f"h_chi2ndf_{tt}_{run_or_mcSet}",
        f"{system(tt)} {algorithm(tt)} ({run_or_mcSet});#chi^{{2}}/ndf;",
        chi2ndf_bins, chi2ndf_min, chi2ndf_max
    )
    roostyling.axes(h_chi2ndf)
    return h_chi2ndf


get_h = {
    'x':       get_h_x,
    'y':       get_h_y,
    'x.y':     get_h_xy,
    'xz':      get_h_xz,
    'yz':      get_h_yz,
    'xz.yz':   get_h_xzyz,
    'n':       get_h_n,
    'trkP':    get_h_trkP,
    'chi2':    get_h_chi2,
    'chi2ndf': get_h_chi2ndf
}



prpts_all = tuple(get_h.keys())
prpts_all_sim = ("x", "y", "x.y", "xz", "yz", "xz.yz", "n")
