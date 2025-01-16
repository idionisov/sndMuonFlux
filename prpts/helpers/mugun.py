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



def getRoundedE(e: float) -> Union[int, float]:
    e = round(e*2)/2
    e = int(e) if e.is_integer() else e
    return e
