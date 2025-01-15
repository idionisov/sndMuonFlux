from typing import Union
import numpy as np
from ROOT import TH1F, TH2F
import roostyling
from ddf.snd.trk import sys, sys_name, alg, alg_name, n_name, get_anti_tt

from ddf.snd.trk import xy_eff_range, xy_full_range


def get_h_n_eff(run_or_mcSet: Union[int, str], tt: int):
    if   tt==1:  n_min = 5; n_max = 20;  nbins_n = 15
    elif tt==11: n_min = 5; n_max = 60;  nbins_n = 55
    elif tt==3:  n_min = 5; n_max = 55;  nbins_n = 50
    elif tt==13: n_min = 5; n_max = 125; nbins_n = 120
    else: raise ValueError(f"{tt} is an invalid track type!")
    
    h_n_eff = tuple(
        TH1F(
            f"h{i}_n_{tt}_{run_or_mcSet}",
            f"{n_name(tt)} {run_or_mcSet};Number of {n_name(tt)};Tracking Efficiency",
            nbins_n, n_min, n_max
        ) for i in range(1, 3)
    )

    roostyling.axes(h_n_eff[0], h_n_eff[1])
    return h_n_eff
