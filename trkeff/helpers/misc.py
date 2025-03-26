from typing import Union
from array import array
from ddfRoot import DdfEff, getGraphFromTEff, getAsNumpy
from ddfUtils import getEffWithError
from hists import xy_eff_range
from sndUtils import system, algorithm
import ROOT
import numpy as np

def getFitEq(
    teff: dict, runOrMC: int=7080, track_types=(1, 11, 3, 13)
):
    eq = {
        tt: {
            xy: ROOT.TF1(f"fit_{xy}_{tt}_{runOrMC}", "[0]", xy_eff_range["min"][xy], xy_eff_range["max"][xy])
            for xy in ["x", "y"]
        } for tt in track_types
    }

    for tt in track_types:
        for xy in ["x", "y"]:
            eq[tt][xy].SetParameters(0.9)
            if isinstance(teff[tt][xy], DdfEff):
                teff[tt][xy].GetGraph().Fit(eq[tt][xy], "SRQ0+")
            elif isinstance(teff[tt][xy], ROOT.TEfficiency):
                getGraphFromTEff(teff[tt][xy]).Fit(eq[tt][xy], "SRQ0+")
            else:
                continue
    return eq


# def saveEffsData(
#     passed: dict,
#     total: dict,
#     fout: ROOT.TFile,
#     statOption: str = "clopper pearson",
#     suffix: str = "tc"
# ):

#     if passed.keys()!=total.keys():
#         raise ValueError("Incompatible passed/total dictionaries!")

#     track_types = tuple(passed.keys())

#     eff = {}
#     for tt in track_types:
#         eff[tt] = {}

#         if not suffix:
#             treeName = f"eff_{tt}"
#         else:
#             treeName = f"eff_{tt}_{suffix}"

#         eff[tt]['tree'] = ROOT.TTree(treeName,  f"{system(tt)} {algorithm(tt)} efficiency")

#         eff[tt]['eff']      = array('f', [ 0. ])
#         eff[tt]['deff_up']  = array('f', [ 0. ])
#         eff[tt]['deff_low'] = array('f', [ 0. ])

#         eff[tt]['eff'][0], eff[tt]['deff_up'][0], eff[tt]['deff_low'][0] = getEffWithError(passed[tt], total[tt], statOption=statOption)
#         print(f" >> {tt}:\t{(eff[tt]['eff'][0], eff[tt]['deff_up'][0], eff[tt]['deff_low'][0])}")

#         eff[tt]['tree'].Branch(f"efficiency", eff[tt]['eff'],      "efficiency/F")
#         eff[tt]['tree'].Branch("deff_up",     eff[tt]['deff_up'],  "deff_up/F")
#         eff[tt]['tree'].Branch("deff_low",    eff[tt]['deff_low'], "deff_low/F")

#         eff[tt]['tree'].Fill()
#         fout.cd()
#         eff[tt]['tree'].Write()






def saveEffsData(
    hists: dict,
    fout: ROOT.TFile,
    statOption: str = "clopper pearson",
    suffix: str = "tc"
):
    track_types = (1, 11, 3, 13)

    eff = {}
    for tt in track_types:
        eff[tt] = {}

        hist = hists[tt]['x.y']
        if not suffix:
            treeName = f"eff_{tt}"
        else:
            treeName = f"eff_{tt}_{suffix}"

        eff[tt]['tree'] = ROOT.TTree(treeName,  f"{system(tt)} {algorithm(tt)} efficiency")

        eff[tt]['eff']    = array('f', [ 0. ])
        eff[tt]['effErr'] = array('f', [ 0. ])

        _eff, _deff, _eff_xEdges, _eff_yEdges = getAsNumpy(hist,
            xmin=xy_eff_range["min"]["x"],
            xmax=xy_eff_range["max"]["x"],
            ymin=xy_eff_range["min"]["y"],
            ymax=xy_eff_range["max"]["y"]
        )

        eff[tt]['eff'][0] = np.mean(_eff)
        eff[tt]['effErr'][0] = np.std(_eff)
        print(f" >> {tt}:\t{eff[tt]['eff'][0]:.03f} ± {eff[tt]['effErr'][0]:.03f}")

        eff[tt]['tree'].Branch("eff",    eff[tt]['eff'],     "eff/F")
        eff[tt]['tree'].Branch("effErr", eff[tt]['effErr'],  "effErr/F")

        eff[tt]['tree'].Fill()
        fout.cd()
        eff[tt]['tree'].Write()



def saveEffsMugun(
    passed:     dict,
    total:      dict,
    fout:       ROOT.TFile,
    statOption: str = "clopper pearson",
    suffix:     str = "muGun.tc"
):
    if (
        list(passed.keys()) != [1, 11, 3, 13] or
        list(total.keys()) != [1, 11, 3, 13]
    ):
        raise ValueError("Invalid passed/total dictionaries!")

    eff = {}
    for tt in passed:
        eff[tt] = {}

        if not suffix:
            treeName = f"eff_{tt}"
        else:
            treeName = f"eff_{tt}_{suffix}"

        eff[tt]['tree'] = ROOT.TTree(treeName,  f"{system(tt)} {algorithm(tt)} efficiency")

        eff[tt]['eff']       = array('f', [ 0. ])
        eff[tt]['effErrUp']  = array('f', [ 0. ])
        eff[tt]['effErrLow'] = array('f', [ 0. ])

        _eff, _effErrLow, _effErrUp = getEffWithError(passed[tt], total[tt], statOption)

        eff[tt]['eff'][0] = _eff
        eff[tt]['effErrLow'][0] = _effErrLow
        eff[tt]['effErrUp'][0] = _effErrUp
        print(f" >> {tt}:\t{eff[tt]['eff'][0]:.03f} (+{eff[tt]['effErrUp'][0]:.03f}, -{eff[tt]['effErrLow'][0]:.03f})")

        eff[tt]['tree'].Branch("eff",       eff[tt]['eff'],       "eff/F")
        eff[tt]['tree'].Branch("effErrLow", eff[tt]['effErrLow'], "effErrLow/F")
        eff[tt]['tree'].Branch("effErrUp",  eff[tt]['effErrUp'],  "effErrUp/F")

        eff[tt]['tree'].Fill()
        fout.cd()
        eff[tt]['tree'].Write()
