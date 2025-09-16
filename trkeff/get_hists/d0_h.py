from typing import Union
import numpy as np
from ROOT import TH1F, TH2F
import roostyling
from sndUtils import sys, system, alg, algorithm, nName, att


d0_min:   float = 12.
d0_max:   float = 76.
nbins_d0: int   = 32


def get_h_d0_eff(
    run_or_mcSet: Union[int, str],
    tt: int
):
    if tt not in (1, 11, 3, 13): raise ValueError(f"{tt} is an invalid track type!")

    h_d0_eff = tuple(
        TH1F(
            f"h{i}_d0_{tt}_{run_or_mcSet}",
            f"{i} {tt} {run_or_mcSet};Distance to beam axis [cm];Tracking Efficiency",
            nbins_d0, d0_min, d0_max
        ) for i in range(1, 3)
    )
    roostyling.axes(h_d0_eff[0], h_d0_eff[1])
    return h_d0_eff



# def get_h_dRef(
#     run_or_mcSet: Union[int, str],
#     tt: int
# ):
#     if tt not in (1, 11, 3, 13):
#         raise ValueError(f"{tt} is an invalid track type!")

#     h_dRef = TH1F(
#         f"h{i}_dRef_{tt}_{run_or_mcSet}",
#         f"{i} {tt} {run_or_mcSet};Tag and candidate track distance at #Z_{{ref}} [cm];",
#         100, 0, 10
#     )

#     roostyling.axes(h_dRef)
#     return h_dRef
