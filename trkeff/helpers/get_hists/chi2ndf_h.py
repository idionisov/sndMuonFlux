from typing import Union
import numpy as np
from ROOT import TH1F, TH2F
import roostyling
from ddf.snd.trk import sys, sys_name, alg, alg_name, n_name, get_anti_tt

from ddf.snd.trk import xy_eff_range, xy_full_range




def get_h_chi2ndf_eff(run_or_mcSet: Union[int, str], tt: int):
    if   tt==1:  chi2ndf_max = 30.;  nbins_chi2ndf = 15
    elif tt==11: chi2ndf_max = 20.;  nbins_chi2ndf = 20
    elif tt==3:  chi2ndf_max = 125.; nbins_chi2ndf = 50
    elif tt==13: chi2ndf_max = 250.; nbins_chi2ndf = 50
    else: raise ValueError(f"{tt} is an invalid track type!")
    
    h_chi2ndf_eff = tuple(
        TH1F(
            f"h{i}_chi2ndf_{tt}_{run_or_mcSet}",
            f"{sys_name(tt)} {alg_name(tt)} {run_or_mcSet};#chi2/ndf;Tracking Efficiency",
            nbins_chi2ndf, 0, chi2ndf_max
        ) for i in range(1, 3)
    )

    roostyling.axes(h_chi2ndf_eff[0], h_chi2ndf_eff[1])
    return h_chi2ndf_eff
