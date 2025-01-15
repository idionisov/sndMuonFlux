from typing import Union
import numpy as np
from ROOT import TH1F, TH2F
import roostyling
from sndUtils import sys, system, alg, algorithm, nName, att


def get_h_trkP_eff(run_or_mcSet: Union[int, str], tt: int):
    if   tt==1:  trkP_min = 6; trkP_max = 10; nbins_trkP = 4
    elif tt==11: trkP_min = 6; trkP_max = 21; nbins_trkP = 15
    elif tt==3:  trkP_min = 6; trkP_max = 26; nbins_trkP = 20
    elif tt==13: trkP_min = 6; trkP_max = 41; nbins_trkP = 35
    else: raise valueerror(f"{tt} is an invalid track type!")

    _att = att(tt)
    h_trkPoints_eff = tuple(
        TH1F(
            f"h{i}_trkP_{tt}_{run_or_mcSet}",
            f"{system(tt)} {algorithm(tt)} {run_or_mcSet};Number of tag track's points (utilized {nName(_att)});Tracking Efficiency",
            nbins_trkP, trkP_min, trkP_max
        ) for i in range(1, 3)
    )
    roostyling.axes(h_trkPoints_eff[0], h_trkPoints_eff[1])
    return h_trkPoints_eff
