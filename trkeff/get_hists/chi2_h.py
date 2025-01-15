from typing import Union
import numpy as np
from ROOT import TH1F, TH2F
import roostyling
from sndUtils import sys, system, alg, algorithm, nName, att


def get_h_chi2_eff(run_or_mcSet: Union[int, str], tt: int):
    if   tt==1:  chi2_max = 60.;   nbins_chi2 = 30
    elif tt==11: chi2_max = 90.;   nbins_chi2 = 45
    elif tt==3:  chi2_max = 1000.; nbins_chi2 = 200
    elif tt==13: chi2_max = 1200.; nbins_chi2 = 300
    else: raise ValueError(f"{tt} is an invalid track type!")

    h_chi2_eff = tuple(
        TH1F(
            f"h{i}_chi2_{tt}_{run_or_mcSet}",
            f"{system(tt)} {algorithm(tt)} {run_or_mcSet};Tag track #chi2 ({sys(att(tt)).upper()} {alg(att(tt)).upper()});Tracking Efficiency",
            nbins_chi2, 0, chi2_max
        ) for i in range(1, 3)
    )

    roostyling.axes(h_chi2_eff[0], h_chi2_eff[1])
    return h_chi2_eff
