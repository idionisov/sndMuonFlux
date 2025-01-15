from typing import Union
import numpy as np
from ROOT import TH1F, TH2F
import roostyling
from ddf.snd.trk import sys, sys_name, alg, alg_name, n_name, get_anti_tt

from ddf.snd.trk import xy_eff_range, xy_full_range



def get_h_chi2_eff(run_or_mcSet: Union[int, str], tt: int):
    if   tt==1:  chi2_max = 60.;   nbins_chi2 = 30
    elif tt==11: chi2_max = 100.;  nbins_chi2 = 50
    elif tt==3:  chi2_max = 600.;  nbins_chi2 = 150
    elif tt==13: chi2_max = 800.;  nbins_chi2 = 200
    else: raise ValueError(f"{tt} is an invalid track type!")

    h_chi2_eff = tuple(
        TH1F(
            f"h{i}_chi2_{tt}_{run_or_mcSet}",
            f"{sys_name(tt)} {alg_name(tt)} {run_or_mcSet};#chi2;Tracking Efficiency",
            nbins_chi2, 0, chi2_max
        ) for i in range(1, 3)
    )

    roostyling.axes(h_chi2_eff[0], h_chi2_eff[1])
    return h_chi2_eff
