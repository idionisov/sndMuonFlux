import ROOT
from typing import Union
import numpy as np
from sndUtils import DdfTrack, DdfMCTrack, sfTrackIsReconstructible, dsTrackIsReconstructible
from ddfRoot import getN
from get_hists.x_h import *
from get_hists.y_h import *
from get_hists.xz_h import *
from get_hists.yz_h import *
from get_hists.chi2_h import *
from get_hists.chi2ndf_h import *
from get_hists.trkP_h import *
from get_hists.n_h import *
from get_hists.d0_h import *
from get_hists.e_h import *

xy_range = {
    'min': {'x': -70., 'y': 0. },
    'max': {'x':  10., 'y': 80.}
}
xy_eff_range = {
    'min': {'x': -42., 'y': 19.},
    'max': {'x': -10., 'y': 48.}
}

getHists: dict = {
    'x':         get_h_x_eff,
    'y':         get_h_y_eff,
    'xz':        get_h_xz_eff,
    'yz':        get_h_yz_eff,
    'chi2':      get_h_chi2_eff,
    'chi2ndf':   get_h_chi2ndf_eff,
    'trkP':      get_h_trkP_eff,
    'n':         get_h_n_eff,
    'd0':        get_h_d0_eff,
    'e':         get_h_e_eff,

    'x.y':       get_h_x_y_eff,
    'x.xz':      get_h_x_xz_eff,
    'x.yz':      get_h_x_yz_eff,
    'x.chi2':    get_h_x_chi2_eff,
    'x.chi2ndf': get_h_x_chi2ndf_eff,
    'x.trkP':    get_h_x_trkP_eff,
    'x.n':       get_h_x_n_eff,

    'y.xz':      get_h_y_xz_eff,
    'y.yz':      get_h_y_yz_eff,
    'y.chi2':    get_h_y_chi2_eff,
    'y.chi2ndf': get_h_y_chi2ndf_eff,
    'y.trkP':    get_h_y_trkP_eff,
    'y.n':       get_h_y_n_eff,

    'xz.x':       get_h_xz_x_eff,
    'xz.y':       get_h_xz_y_eff,
    'xz.yz':      get_h_xz_yz_eff,
    'xz.chi2':    get_h_xz_chi2_eff,
    'xz.chi2ndf': get_h_xz_chi2ndf_eff,
    'xz.trkP':    get_h_xz_trkP_eff,
    'xz.n':       get_h_xz_n_eff,

    'yz.x':       get_h_yz_x_eff,
    'yz.y':       get_h_yz_y_eff,
    'yz.chi2':    get_h_yz_chi2_eff,
    'yz.chi2ndf': get_h_yz_chi2ndf_eff,
    'yz.trkP':    get_h_yz_trkP_eff,
    'yz.n':       get_h_yz_n_eff,

    'dxRef':      get_h_dxRef,
    'dyRef':      get_h_dyRef,

    'dxz':        get_h_dxz,
    'dyz':        get_h_dyz
}

prpts_data_all = (
    'x', 'y', 'xz', 'yz', 'chi2', 'chi2ndf', 'trkP', 'n', 'd0',
    'x.y', 'x.xz', 'x.yz', 'x.chi2', 'x.chi2ndf', 'x.trkP', 'x.n',
    'y.xz', 'y.yz', 'y.chi2', 'y.chi2ndf', 'y.trkP', 'y.n',
    'xz.x', 'xz.y', 'xz.yz', 'xz.chi2', 'xz.chi2ndf', 'xz.trkP', 'xz.n',
    'yz.x', 'yz.y', 'yz.chi2', 'yz.chi2ndf', 'yz.trkP', 'yz.n',
    'dxRef', 'dyRef', 'dxz', 'dyz'
)

prpts_sim_all = (
    'x', 'y', 'xz', 'yz', 'n', 'd0',
    'x.y', 'x.xz', 'x.yz', 'x.n',
    'y.xz', 'y.yz', 'y.n',
    'xz.x', 'xz.y', 'xz.yz', 'xz.n',
    'yz.x', 'yz.y', 'yz.n',
    'dxRef', 'dyRef', 'dxz', 'dyz'
)


def createHists(
    run_or_mcSet = 7080,
    track_types = (1, 11, 3, 13),
    prpts = prpts_data_all
):
    h = {}
    for tt in track_types:
        h[tt] = {}

        for prpt in prpts:
            h[tt][prpt] = getHists[prpt](run_or_mcSet, tt)
    return h

def isWithinFiducialArea(
    point: ROOT.TVector3,
    xmin: float = -42.,
    xmax: float = -10.,
    ymin: float =  19.,
    ymax: float =  48.
) -> bool:
    if (
        point.X() >= xmin and
        point.X() <= xmax and
        point.Y() >= ymin and
        point.Y() <= ymax
    ): return True
    else: return False



def areWithinAllowedDistance(
    point1: ROOT.TVector3,
    point2: ROOT.TVector3,
    distance: float = 3.
) -> bool:
    if (
        abs(point2.X()-point1.X()) <= distance and
        abs(point2.Y()-point1.Y()) <= distance
    ): return True
    else: return False





def fillHistsTC(
    h:            dict,
    tag_trk:      DdfTrack,
    run_or_mcSet: Union[int, str] = 7080,
    z_ref:        float = 450.,
    tt:           int = 1,
    ip1_angle:    float = 0.02,
    weight:       Union[int, float] = 1,
) -> dict:
    result = {"total": False, "passed": False}
    ##if tag_trk.tt != att(tt):
    #    print(f"      NO! {tt}/{tag_trk.tt} \t {tag_trk.att()}|{att(tag_trk.tt)}")
    #    return result
    #else:
    event = tag_trk.Event

    ref_tag = tag_trk.GetPointAtZ(z_ref)
    x_tag = ref_tag.X()
    y_tag = ref_tag.Y()

    d0 = np.sqrt(x_tag*x_tag + y_tag*y_tag)

    xz_tag = tag_trk.XZ*1e3
    yz_tag = tag_trk.YZ*1e3

    chi2_tag = tag_trk.Chi2
    chi2ndf_tag = tag_trk.Chi2Ndf

    trkP_tag = tag_trk.GetPoints().size()
    n = getN(tt, event)

    prpts = h[tt].keys()

    if "x.y" in prpts: h[tt]["x.y"][1].Fill(x_tag, y_tag, weight)
    if isWithinFiducialArea(ref_tag,
        xmin=xy_eff_range["min"]["x"],
        xmax=xy_eff_range["max"]["x"],
        ymin=xy_eff_range["min"]["y"],
        ymax=xy_eff_range["max"]["y"]
    ):
        if "x"       in prpts: h[tt]["x"][1].Fill(x_tag, weight)
        if "y"       in prpts: h[tt]["y"][1].Fill(y_tag, weight)
        if "xz"      in prpts: h[tt]['xz'][1].Fill(xz_tag, weight)
        if "yz"      in prpts: h[tt]['yz'][1].Fill(yz_tag, weight)
        if "chi2"    in prpts: h[tt]['chi2'][1].Fill(chi2_tag, weight)
        if "chi2ndf" in prpts: h[tt]['chi2ndf'][1].Fill(chi2ndf_tag, weight)
        if "n"       in prpts: h[tt]['n'][1].Fill(n, weight)
        if "trkP"    in prpts: h[tt]['trkP'][1].Fill(trkP_tag, weight)
        if "d0"      in prpts: h[tt]['d0'][1].Fill(d0, weight)

        if "x.xz"       in prpts: h[tt]['x.xz'][1].Fill(x_tag, xz_tag, weight)
        if "x.yz"       in prpts: h[tt]['x.yz'][1].Fill(x_tag, yz_tag, weight)
        if "x.chi2"     in prpts: h[tt]['x.chi2'][1].Fill(x_tag, chi2_tag, weight)
        if "x.chi2ndf"  in prpts: h[tt]['x.chi2ndf'][1].Fill(x_tag, chi2ndf_tag, weight)
        if "x.trkP"     in prpts: h[tt]['x.trkP'][1].Fill(x_tag, trkP_tag, weight)
        if "x.n"        in prpts: h[tt]['x.n'][1].Fill(x_tag, n, weight)

        if "y.xz"       in prpts: h[tt]['y.xz'][1].Fill(y_tag, xz_tag, weight)
        if "y.yz"       in prpts: h[tt]['y.yz'][1].Fill(y_tag, yz_tag, weight)
        if "y.chi2"     in prpts: h[tt]['y.chi2'][1].Fill(y_tag, chi2_tag, weight)
        if "y.chi2ndf"  in prpts: h[tt]['y.chi2ndf'][1].Fill(y_tag, chi2ndf_tag, weight)
        if "y.trkP"     in prpts: h[tt]['y.trkP'][1].Fill(y_tag, trkP_tag, weight)
        if "y.n"        in prpts: h[tt]['y.n'][1].Fill(y_tag, n, weight)

        if "xz.x"       in prpts: h[tt]['xz.x'][1].Fill(xz_tag, x_tag, weight)
        if "xz.y"       in prpts: h[tt]['xz.y'][1].Fill(xz_tag, y_tag, weight)
        if "xz.yz"      in prpts: h[tt]["xz.yz"][1].Fill(xz_tag, yz_tag, weight)
        if "xz.chi2"    in prpts: h[tt]['xz.chi2'][1].Fill(xz_tag, chi2_tag, weight)
        if "xz.chi2ndf" in prpts: h[tt]['xz.chi2ndf'][1].Fill(xz_tag, chi2ndf_tag, weight)
        if "xz.trkP"    in prpts: h[tt]['xz.trkP'][1].Fill(xz_tag, trkP_tag, weight)
        if "xz.n"       in prpts: h[tt]['xz.n'][1].Fill(xz_tag, n, weight)

        if "yz.x"       in prpts: h[tt]['yz.x'][1].Fill(yz_tag, x_tag, weight)
        if "yz.y"       in prpts: h[tt]['yz.y'][1].Fill(yz_tag, y_tag, weight)
        if "yz.chi2"    in prpts: h[tt]['yz.chi2'][1].Fill(yz_tag, chi2_tag, weight)
        if "yz.chi2ndf" in prpts: h[tt]['yz.chi2ndf'][1].Fill(yz_tag, chi2ndf_tag, weight)
        if "yz.trkP"    in prpts: h[tt]['yz.trkP'][1].Fill(yz_tag, trkP_tag, weight)
        if "yz.n"       in prpts: h[tt]['yz.n'][1].Fill(yz_tag, n, weight)

        result["total"] = True

    for trk2 in event.Reco_MuonTracks:
        trk2 = DdfTrack(Track=trk2, Event=event, IP1_Angle=0.08)

        if not (
            trk2.tt==tt and
            trk2.IsIP1()
            #abs(trk2.XZ) <= 0.02 and
            #abs(trk2.YZ) <= 0.02
            #abs(tag_trk.XZ-trk2.XZ) <= 0.02 and
            #abs(tag_trk.YZ-trk2.YZ) <= 0.02
        ): continue

        xz_cand = 1e3*trk2.XZ
        yz_cand = 1e3*trk2.YZ
        ref2 = trk2.GetPointAtZ(z_ref)


        if "dxRef" in prpts: h[tt]['dxRef'].Fill(ref2.X()-x_tag, weight)
        if "dyRef" in prpts: h[tt]['dyRef'].Fill(ref2.Y()-y_tag, weight)
        if "dxz"   in prpts: h[tt]['dxz'].Fill(xz_cand - xz_tag, weight)
        if "dyz"   in prpts: h[tt]['dyz'].Fill(yz_cand - yz_tag, weight)

        if not areWithinAllowedDistance(ref_tag, ref2, 3.):
            continue

        if "x.y" in prpts: h[tt]["x.y"][0].Fill(x_tag, y_tag, weight)

        if not isWithinFiducialArea(ref_tag,
            xmin=xy_eff_range["min"]["x"],
            xmax=xy_eff_range["max"]["x"],
            ymin=xy_eff_range["min"]["y"],
            ymax=xy_eff_range["max"]["y"]
        ): continue

        if "x"       in prpts: h[tt]["x"][0].Fill(x_tag, weight)
        if "y"       in prpts: h[tt]["y"][0].Fill(y_tag, weight)
        if "xz"      in prpts: h[tt]['xz'][0].Fill(xz_tag, weight)
        if "yz"      in prpts: h[tt]['yz'][0].Fill(yz_tag, weight)
        if "chi2"    in prpts: h[tt]['chi2'][0].Fill(chi2_tag, weight)
        if "chi2ndf" in prpts: h[tt]['chi2ndf'][0].Fill(chi2ndf_tag, weight)
        if "n"       in prpts: h[tt]['n'][0].Fill(n, weight)
        if "trkP"    in prpts: h[tt]['trkP'][0].Fill(trkP_tag, weight)
        if "d0"      in prpts: h[tt]['d0'][0].Fill(d0, weight)

        if "x.xz"       in prpts: h[tt]['x.xz'][0].Fill(x_tag, xz_tag, weight)
        if "x.yz"       in prpts: h[tt]['x.yz'][0].Fill(x_tag, yz_tag, weight)
        if "x.chi2"     in prpts: h[tt]['x.chi2'][0].Fill(x_tag, chi2_tag, weight)
        if "x.chi2ndf"  in prpts: h[tt]['x.chi2ndf'][0].Fill(x_tag, chi2ndf_tag, weight)
        if "x.trkP"     in prpts: h[tt]['x.trkP'][0].Fill(x_tag, trkP_tag, weight)
        if "x.n"        in prpts: h[tt]['x.n'][0].Fill(x_tag, n, weight)

        if "y.xz"       in prpts: h[tt]['y.xz'][0].Fill(y_tag, xz_tag, weight)
        if "y.yz"       in prpts: h[tt]['y.yz'][0].Fill(y_tag, yz_tag, weight)
        if "y.chi2"     in prpts: h[tt]['y.chi2'][0].Fill(y_tag, chi2_tag, weight)
        if "y.chi2ndf"  in prpts: h[tt]['y.chi2ndf'][0].Fill(y_tag, chi2ndf_tag, weight)
        if "y.trkP"     in prpts: h[tt]['y.trkP'][0].Fill(y_tag, trkP_tag, weight)
        if "y.n"        in prpts: h[tt]['y.n'][0].Fill(y_tag, n, weight)

        if "xz.x"       in prpts: h[tt]['xz.x'][0].Fill(xz_tag, x_tag, weight)
        if "xz.y"       in prpts: h[tt]['xz.y'][0].Fill(xz_tag, y_tag, weight)
        if "xz.yz"      in prpts: h[tt]["xz.yz"][0].Fill(xz_tag, yz_tag, weight)
        if "xz.chi2"    in prpts: h[tt]['xz.chi2'][0].Fill(xz_tag, chi2_tag, weight)
        if "xz.chi2ndf" in prpts: h[tt]['xz.chi2ndf'][0].Fill(xz_tag, chi2ndf_tag, weight)
        if "xz.trkP"    in prpts: h[tt]['xz.trkP'][0].Fill(xz_tag, trkP_tag, weight)
        if "xz.n"       in prpts: h[tt]['xz.n'][0].Fill(xz_tag, n, weight)

        if "yz.x"       in prpts: h[tt]['yz.x'][0].Fill(yz_tag, x_tag, weight)
        if "yz.y"       in prpts: h[tt]['yz.y'][0].Fill(yz_tag, y_tag, weight)
        if "yz.chi2"    in prpts: h[tt]['yz.chi2'][0].Fill(yz_tag, chi2_tag, weight)
        if "yz.chi2ndf" in prpts: h[tt]['yz.chi2ndf'][0].Fill(yz_tag, chi2ndf_tag, weight)
        if "yz.trkP"    in prpts: h[tt]['yz.trkP'][0].Fill(yz_tag, trkP_tag, weight)
        if "yz.n"       in prpts: h[tt]['yz.n'][0].Fill(yz_tag, n, weight)

        result["passed"] = True
    return result




def fillHistsRT(
    h:            dict,
    flag:         dict,
    tag_trk:      DdfMCTrack,
    mcSet:        str,
    z_ref:        float = 450.,
    tt:           int = 1,
    ip1_angle:    float = 20.,
    weight:       Union[int, float] = 1
):
    if not flag["total"][tt]:
        return

    event = tag_trk.Event
    ref_tag = tag_trk.GetPointAtZ(z_ref)
    x_tag = ref_tag.X()
    y_tag = ref_tag.Y()

    d0 = np.sqrt(x_tag*x_tag + y_tag*y_tag)

    xz_tag = tag_trk.XZ*1e3
    yz_tag = tag_trk.YZ*1e3

    n = getN(tt, event)

    prpts = h[tt].keys()

    if "x.y" in prpts: h[tt]["x.y"][1].Fill(x_tag, y_tag, weight)
    if isWithinFiducialArea(ref_tag,
        xmin=xy_eff_range["min"]["x"],
        xmax=xy_eff_range["max"]["x"],
        ymin=xy_eff_range["min"]["y"],
        ymax=xy_eff_range["max"]["y"]
    ):
        if "x"       in prpts: h[tt]["x"][1].Fill(x_tag, weight)
        if "y"       in prpts: h[tt]["y"][1].Fill(y_tag, weight)
        if "xz"      in prpts: h[tt]['xz'][1].Fill(xz_tag, weight)
        if "yz"      in prpts: h[tt]['yz'][1].Fill(yz_tag, weight)
        if "n"       in prpts: h[tt]['n'][1].Fill(n, weight)
        if "d0"      in prpts: h[tt]['d0'][1].Fill(d0, weight)

        if "x.xz"       in prpts: h[tt]['x.xz'][1].Fill(x_tag, xz_tag, weight)
        if "x.yz"       in prpts: h[tt]['x.yz'][1].Fill(x_tag, yz_tag, weight)
        if "x.n"        in prpts: h[tt]['x.n'][1].Fill(x_tag, n, weight)

        if "y.xz"       in prpts: h[tt]['y.xz'][1].Fill(y_tag, xz_tag, weight)
        if "y.yz"       in prpts: h[tt]['y.yz'][1].Fill(y_tag, yz_tag, weight)
        if "y.n"        in prpts: h[tt]['y.n'][1].Fill(y_tag, n, weight)

        if "xz.x"       in prpts: h[tt]['xz.x'][1].Fill(xz_tag, x_tag, weight)
        if "xz.y"       in prpts: h[tt]['xz.y'][1].Fill(xz_tag, y_tag, weight)
        if "xz.yz"      in prpts: h[tt]["xz.yz"][1].Fill(xz_tag, yz_tag, weight)
        if "xz.n"       in prpts: h[tt]['xz.n'][1].Fill(xz_tag, n, weight)

        if "yz.x"       in prpts: h[tt]['yz.x'][1].Fill(yz_tag, x_tag, weight)
        if "yz.y"       in prpts: h[tt]['yz.y'][1].Fill(yz_tag, y_tag, weight)
        if "yz.n"       in prpts: h[tt]['yz.n'][1].Fill(yz_tag, n, weight)

    for trk2 in event.Reco_MuonTracks:
        trk2 = DdfTrack(Track=trk2, Event=event, IP1_Angle=ip1_angle)

        if not (
            trk2.tt==tt and
            trk2.IsGood(xz_min=-ip1_angle/1e3, xz_max=ip1_angle/1e3, yz_min=-ip1_angle/1e3, yz_max=ip1_angle/1e3) and
            flag["passed"][tt]
        ): continue

        xz_cand = 1e3*trk2.XZ
        yz_cand = 1e3*trk2.YZ
        ref2 = trk2.GetPointAtZ(z_ref)

        if "dxRef" in prpts: h[tt]['dxRef'].Fill(ref2.X()-x_tag, weight)
        if "dyRef" in prpts: h[tt]['dyRef'].Fill(ref2.Y()-y_tag, weight)
        if "dxz"   in prpts: h[tt]['dxz'].Fill(xz_cand - xz_tag, weight)
        if "dyz"   in prpts: h[tt]['dyz'].Fill(yz_cand - yz_tag, weight)

        if not areWithinAllowedDistance(ref_tag, ref2, 3.):
            continue

        if "x.y" in prpts: h[tt]["x.y"][0].Fill(x_tag, y_tag, weight)

        if not isWithinFiducialArea(ref_tag,
            xmin=xy_eff_range["min"]["x"],
            xmax=xy_eff_range["max"]["x"],
            ymin=xy_eff_range["min"]["y"],
            ymax=xy_eff_range["max"]["y"]
        ): continue

        if "x"       in prpts: h[tt]["x"][0].Fill(x_tag, weight)
        if "y"       in prpts: h[tt]["y"][0].Fill(y_tag, weight)
        if "xz"      in prpts: h[tt]['xz'][0].Fill(xz_tag, weight)
        if "yz"      in prpts: h[tt]['yz'][0].Fill(yz_tag, weight)
        if "n"       in prpts: h[tt]['n'][0].Fill(n, weight)
        if "d0"      in prpts: h[tt]['d0'][0].Fill(d0, weight)

        if "x.xz"       in prpts: h[tt]['x.xz'][0].Fill(x_tag, xz_tag, weight)
        if "x.yz"       in prpts: h[tt]['x.yz'][0].Fill(x_tag, yz_tag, weight)
        if "x.n"        in prpts: h[tt]['x.n'][0].Fill(x_tag, n, weight)

        if "y.xz"       in prpts: h[tt]['y.xz'][0].Fill(y_tag, xz_tag, weight)
        if "y.yz"       in prpts: h[tt]['y.yz'][0].Fill(y_tag, yz_tag, weight)
        if "y.n"        in prpts: h[tt]['y.n'][0].Fill(y_tag, n, weight)

        if "xz.x"       in prpts: h[tt]['xz.x'][0].Fill(xz_tag, x_tag, weight)
        if "xz.y"       in prpts: h[tt]['xz.y'][0].Fill(xz_tag, y_tag, weight)
        if "xz.yz"      in prpts: h[tt]["xz.yz"][0].Fill(xz_tag, yz_tag, weight)
        if "xz.n"       in prpts: h[tt]['xz.n'][0].Fill(xz_tag, n, weight)

        if "yz.x"       in prpts: h[tt]['yz.x'][0].Fill(yz_tag, x_tag, weight)
        if "yz.y"       in prpts: h[tt]['yz.y'][0].Fill(yz_tag, y_tag, weight)
        if "yz.n"       in prpts: h[tt]['yz.n'][0].Fill(yz_tag, n, weight)



def getEffRT(
    event:        ROOT.TChain,
    h:            dict,
    mcSet:        str,
    z_ref:        dict = {1: 430., 11: 450, 3: 430., 13: 450},
    ip1_angle:    float = 20.,
    weight:       Union[int, float] = 1,
    track_types:  tuple = (1, 11, 3, 13),
    xz_min:       float = -20.,
    yz_min:       float = -20.,
    xz_max:       float =  20.,
    yz_max:       float =  20.,
):
    _sf = sfTrackIsReconstructible(event)
    _ds = dsTrackIsReconstructible(event)

    flag = {
        "passed": {tt: False for tt in track_types},
        "total":  {tt: False for tt in track_types}
    }
    if _sf==False and _ds==False:
        return flag

    reco = {1: _sf, 11: _sf, 3: _ds, 13: _ds}

    for tt in track_types:
        if not reco[tt]:
            continue

        flag["total"][tt] = True

        for trk in event.Reco_MuonTracks:
            if trk.getTrackType() != tt:
                continue
            flag["passed"][tt] = True

            for mcTrack in event.MCTrack:
                ddfMCTrack = DdfMCTrack(mcTrack, Event=event, IP1_Angle=20.)
                if not (
                    ddfMCTrack.XZ <= xz_max/1e3 and
                    ddfMCTrack.XZ >= xz_min/1e3 and
                    ddfMCTrack.YZ <= yz_max/1e3 and
                    ddfMCTrack.YZ >= yz_min/1e3
                ): continue

                if not (ddfMCTrack.IsWithinDS3() and ddfMCTrack.IsWithinSF1()):
                    continue

                fillHistsRT(h, flag, ddfMCTrack, "muGun.rt", z_ref[tt], tt)
                break
    return flag



def getEffNRT(
    event:        ROOT.TChain,
    h:            dict,
    mcSet:        str,
    z_ref:        dict = {1: 430., 11: 450, 3: 430., 13: 450},
    ip1_angle:    float = 20.,
    weight:       Union[int, float] = 1,
    track_types:  tuple = (1, 11, 3, 13),
    xz_min:       float = -20.,
    yz_min:       float = -20.,
    xz_max:       float =  20.,
    yz_max:       float =  20.,
):
    _sf = sfTrackIsReconstructible(event)
    _ds = dsTrackIsReconstructible(event)

    flag = {
        "passed": {tt: False for tt in track_types},
        "total":  {tt: False for tt in track_types}
    }
    if _sf==False and _ds==False:
        return flag

    reco = {1: _sf, 11: _sf, 3: _ds, 13: _ds}

    for tt in track_types:
        if not reco[tt]:
            continue
        n = getN(tt, event)

        h[tt]['n'][1].Fill(n, weight)
        flag["total"][tt] = True

        for trk in event.Reco_MuonTracks:
            if trk.getTrackType() != tt:
                continue
            flag["passed"][tt] = True
            h[tt]['n'][0].Fill(n, weight)

    return flag



def getEff(
    hist: ROOT.TH2,
    xmin: float = -42.,
    xmax: float = -10.,
    ymin: float = 19.,
    ymax: float = 48.,
) -> np.ndarray:
    x_bin_edges = [hist.GetXaxis().GetBinLowEdge(i) for i in range(1, hist.GetNbinsX() + 2)]
    y_bin_edges = [hist.GetYaxis().GetBinLowEdge(i) for i in range(1, hist.GetNbinsY() + 2)]

    filtered_bins = []
    for x_bin in range(1, hist.GetNbinsX() + 1):
        x_low = hist.GetXaxis().GetBinLowEdge(x_bin)
        x_up = hist.GetXaxis().GetBinUpEdge(x_bin)

        if x_low < xmin or x_up > xmax:
            continue

        y_row = []
        for y_bin in range(1, hist.GetNbinsY() + 1):
            y_low = hist.GetYaxis().GetBinLowEdge(y_bin)
            y_up = hist.GetYaxis().GetBinUpEdge(y_bin)

            if y_low < ymin or y_up > ymax:
                continue

            y_row.append(hist.GetBinContent(x_bin, y_bin))

        if y_row:
            filtered_bins.append(y_row)

    return np.array(filtered_bins, dtype=np.float64)
