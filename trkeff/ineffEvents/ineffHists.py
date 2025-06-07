import ROOT
from sndUtils import *




def getIneffNHists(runOrMC = 7080, TTs = (1, 11, 3, 13)):
    h = {}
    for tt in TTs:
        if   tt==1:  n_min = 5; n_max = 40;  nbins_n = 35
        elif tt==11: n_min = 5; n_max = 80;  nbins_n = 75
        elif tt==3:  n_min = 5; n_max = 55;  nbins_n = 50
        else:        n_min = 5; n_max = 125; nbins_n = 120

        h[tt] = {}
        for rsn in ("all", "eff", "ineff", "notBuilt", "notPaired", "cand>20mrad"):
            h[tt][rsn] = ROOT.TH1F(
                f"h_ineff_n_{tt}_{runOrMC}_{rsn}",
                f";# {nName(tt)};",
                nbins_n, n_min, n_max
            )
    return h

def getIneffXYHists(runOrMC = 7080, TTs = (1, 11, 3, 13)):

    h = {}
    for tt in TTs:
        h[tt] = {}

        for rsn in ("all", "eff", "ineff", "notBuilt", "notPaired", "cand>20mrad"):
            h[tt][rsn] = ROOT.TH2F(
                f"h_ineff_x.y_{tt}_{runOrMC}_{rsn}",
                f";X_{{ref}} (cm);Y_{{ref}} (cm);",
                80, -75, 5, 80, -5, 75
            )

    return h


def getIneffMaxPlaneHitsHists(runOrMC = 7080, TTs = (1, 11, 3, 13)):
    h = {}
    for tt in TTs:
        if   tt==1:  n_min = 0; n_max = 40; nbins_n = 40
        elif tt==11: n_min = 0; n_max = 40; nbins_n = 40
        elif tt==3:  n_min = 0; n_max = 20; nbins_n = 40
        else:        n_min = 0; n_max = 20; nbins_n = 40

        h[tt] = {}
        for rsn in ("all", "eff", "ineff", "notBuilt", "notPaired", "cand>20mrad"):
            h[tt][rsn] = ROOT.TH1F(
                f"h_ineff_maxHits_{tt}_{runOrMC}_{rsn}",
                f";# {nName(tt)};",
                nbins_n, n_min, n_max
            )
    return h
