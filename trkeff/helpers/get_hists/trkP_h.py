from typing import Union
import numpy as np
from ROOT import TH1F, TH2F
import roostyling
from ddf.snd.trk import sys, sys_name, alg, alg_name, n_name, get_anti_tt

from ddf.snd.trk import xy_eff_range, xy_full_range



def get_h_trkP_eff(run_or_mcSet: Union[int, str], tt: int):
    if   tt==1:  trkP_min = 6; trkP_max = 8;  nbins_trkP = 2
    elif tt==11: trkP_min = 6; trkP_max = 16; nbins_trkP = 10
    elif tt==3:  trkP_min = 6; trkP_max = 16; nbins_trkP = 10
    elif tt==13: trkP_min = 6; trkP_max = 41; nbins_trkP = 35
    else: raise ValueError(f"{tt} is an invalid track type!")
    
    att = get_anti_tt(tt)
    h_trkPoints_eff = tuple(
        TH1F(
            f"h{i}_trkP_{tt}_{run_or_mcSet}",
            f"{sys_name(tt)} {alg_name(tt)} {run_or_mcSet};Number of tag track's points (utilized {n_name(att)});Tracking Efficiency",
            nbins_trkP, trkP_min, trkP_max
        ) for i in range(1, 3)
    )
    roostyling.axes(h_trkPoints_eff[0], h_trkPoints_eff[1])
    return h_trkPoints_eff
