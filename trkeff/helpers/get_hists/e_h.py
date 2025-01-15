from typing import Union
import numpy as np
from ROOT import TH1F, TH2F
import roostyling
from ddf.snd.trk import sys, sys_name, alg, alg_name, n_name, get_anti_tt

from ddf.snd.trk import xy_eff_range, xy_full_range


energies = np.array([10, 15, 20, 35, 55, 100, 200, 300, 600, 1010], dtype=np.float64)


def get_h_e_eff(
    run_or_mcSet: Union[int, str],
    tt: int    
):
    if tt not in (1, 11, 3, 13): raise ValueError(f"{tt} is an invalid track type!")
    
    h_e_eff = tuple(
        TH1F(
            f"h{i}_e_{tt}_{run_or_mcSet}",
            f"{sys_name(tt)} {alg_name(tt)} {run_or_mcSet};Muon Energy [GeV];Tracking Efficiency",
            len(energies)-1, energies
        ) for i in range(1, 3)
    )

    roostyling.axes(h_e_eff[0], h_e_eff[1])
    return h_e_eff
