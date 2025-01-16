import ROOT
from typing import Union

xy_full_range = {
    'min': {'x': -75., 'y': -5. },
    'max': {'x':   5., 'y':  75.}
}
xy_eff_range = {
    'min': {'x': -42., 'y': 19.},
    'max': {'x': -10., 'y': 48.}
}



def getTrees(
    track_types = (1, 11, 3, 13)
) -> dict:
    eff = {}

    efficiency = {}
    deff_up    = {}
    deff_low   = {}
    E          = {}
    dE_up      = {}
    dE_low     = {}

    for tt in track_types:
        eff[tt] = ROOT.TTree(f"eff_e_{tt}",  f"{sys_name(tt)} {alg_name(tt)} efficiency dependence on energy")

        efficiency[tt] = array('f', [ 0. ])
        deff_up[tt]    = array('f', [ 0. ])
        deff_low[tt]   = array('f', [ 0. ])
        E[tt]          = array('f', [ 0. ])
        dE_up[tt]      = array('f', [ 0. ])
        dE_low[tt]     = array('f', [ 0. ])

        eff[tt].Branch("efficiency", efficiency[tt], "efficiency/F")
        eff[tt].Branch("deff_up",    deff_up[tt],    "deff_up/F")
        eff[tt].Branch("deff_low",   deff_low[tt],   "deff_low/F")
        eff[tt].Branch("E",          E[tt],          "E/F")
        eff[tt].Branch("dE_up",      dE_up[tt],      "dE_up/F")
        eff[tt].Branch("dE_low_",    dE_low[tt],     "dE_low/F")

    return eff


def getRoundedE(e: float) -> Union[int, float]:
    e = round(e*2)/2
    e = int(e) if e.is_integer() else e
    return e
