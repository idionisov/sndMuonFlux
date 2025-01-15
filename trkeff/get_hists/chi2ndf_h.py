from typing import Union
import numpy as np
from ROOT import TH1F, TH2F
import roostyling
from sndUtils import sys, system, alg, algorithm, nName, att



def get_h_chi2ndf_eff(run_or_mcSet: Union[int, str], tt: int):
    if   tt==1:  chi2ndf_max = 60.;  nbins_chi2ndf = 30
    elif tt==11: chi2ndf_max = 35.;  nbins_chi2ndf = 35
    elif tt==3:  chi2ndf_max = 200.; nbins_chi2ndf = 100
    elif tt==13: chi2ndf_max = 150.; nbins_chi2ndf = 50
    else: raise ValueError(f"{tt} is an invalid track type!")

    h_chi2ndf_eff = tuple(
        TH1F(
            f"h{i}_chi2ndf_{tt}_{run_or_mcSet}",
            f"{system(tt)} {algorithm(tt)} {run_or_mcSet};Tag track #chi2/ndf ({sys(att(tt)).upper()} {sys(att(tt)).upper()});Tracking Efficiency",
            nbins_chi2ndf, 0, chi2ndf_max
        ) for i in range(1, 3)
    )

    roostyling.axes(h_chi2ndf_eff[0], h_chi2ndf_eff[1])
    return h_chi2ndf_eff
