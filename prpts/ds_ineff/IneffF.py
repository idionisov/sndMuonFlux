from ROOT import TMath, TChain, TVector3
from ddf.snd.trk import sys, alg, sys_name, alg_name, get_anti_tt
from ddf.snd.mc import *

def is_us(point) -> bool:
    detID = point.GetDetectorID()
    if detID<3e4 and detID>2e4: return True
    else: return False

def is_ds(point) -> bool:
    detID = point.GetDetectorID()
    if detID>=30000: return True
    else: return False


def get_ds_points(event):
    ds_points = []
    for dsP in event.MuFilterPoint:
        if not (
            abs(dsP.PdgCode())==13 and
            dsP.GetTrackID()==0 and
            is_ds(dsP)
        ): continue
        ds_points.append(dsP)
    return ds_points

def get_sf_points(event):
    sf_points = []
    for sfP in event.ScifiPoint:
        if not (
            event.EventHeader.isIP1() and
            abs(sfP.PdgCode())==13 and
            sfP.GetTrackID()==0
        ): continue
        sf_points.append(sfP)
    return sf_points


def get_sf_dP(event):
    sf_points = get_sf_points(event)

    dP = []
    P_prev = 0
    for sfP_i, sfP in enumerate(sf_points):
        P = sfP.Momentum()    
        if sfP_i==0:
            P_prev = P
            continue

        dP.append(P - P_prev)
    return dP

def get_ds_dP(event):
    ds_points = get_ds_points(event)

    dP = []
    P_prev = 0
    for dsP_i, dsP in enumerate(ds_points):
        P = dsP.Momentum()    
        if dsP_i==0:
            P_prev = P
            continue

        dP.append(P - P_prev)
    return dP


def get_sys_dP(event):
    sfP_last  = get_sf_points(event)[-1]
    dsP_first = get_ds_points(event)[0]

    return dsP_first - sfP_last

def us_scatter(
    event,
    Aref: dict = {
        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
        'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
    },
    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.}
):
    if not (
        event.EventHeader.isIP1() and
        should_be_a_track(event)["sf"] and
        should_be_a_track(event)["ds"]
    ): return {3: False, 13: False}

    sf_points = get_sf_points(event)
    ds_points = get_ds_points(event)

    sfP_last  = sf_points[-1]
    dsP_first = ds_points[0]

    xz_sf = get_xz(sfP_last)
    yz_sf = get_yz(sfP_last)

    xz_ds = get_xz(dsP_first)
    yz_ds = get_yz(dsP_first)

    trks = {tt: False for tt in (1, 11, 3, 13)}
    us_scat = {3: False, 13: False}
    ref = {}

    for tt_sf, tt_ds in zip((1, 11), (3, 13)):
        ref[tt_sf] = get_intersection_z(sfP_last,  z_ref[tt_sf]) 
        ref[tt_ds] = get_intersection_z(dsP_first, z_ref[tt_ds]) 

        x_sf = ref[tt_sf].X()
        y_sf = ref[tt_sf].Y()
        if (
            x_sf >= Aref["sf"]["min"]["x"] and
            x_sf <= Aref["sf"]["max"]["x"] and
            y_sf >= Aref["sf"]["min"]["y"] and
            y_sf <= Aref["sf"]["max"]["y"] and
            abs(xz_sf) <= 0.02 and
            abs(yz_sf) <= 0.02
        ): trks[tt_sf] = True
        else: trks[tt_sf] = False

        x_ds = ref[tt_ds].X()
        y_ds = ref[tt_ds].Y()
        if (
            x_ds >= Aref["ds"]["min"]["x"] and
            x_ds <= Aref["ds"]["max"]["x"] and
            y_ds >= Aref["ds"]["min"]["y"] and
            y_ds <= Aref["ds"]["max"]["y"] and
            abs(xz_ds) <= 0.02 and
            abs(yz_ds) <= 0.02
        ): trks[tt_ds] = True
        else: trks[tt_ds] = False

        if trks[tt_sf] and not trks[tt_ds]:
            us_scat[tt_ds] = True
        else:
            us_scat[tt_ds] = False

    return us_scat
    






def is_within_Aref(
    mc_point, tt: int,
    Aref: dict = {
        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
        'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
    },
    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.}
) -> bool:
    
    if not (abs(mc_point.PdgCode())==13 and mc_point.GetTrackID()==0):
        return False

    z = z_ref[tt]
    x = get_intersection_z(mc_point, z).X()
    y = get_intersection_z(mc_point, z).Y()
    
    if (
        x >= Aref[sys(tt)]["min"]["x"] and
        x <= Aref[sys(tt)]["max"]["x"] and
        y >= Aref[sys(tt)]["min"]["y"] and
        y >= Aref[sys(tt)]["max"]["y"]
    ): return True
    return False


def is_ip1(
    mc_point,
    event
) -> bool:
    if not event.EventHeader.isIP1(): return False

    xz = get_xz(mc_point)
    yz = get_yz(mc_point)

    if (
        mc_point.GetPz() > 0 and
        abs(xz) <= 0.02 and
        abs(yz) <= 0.02
    ): return True
    else: return False



def get_xz(mc_point) -> float:
    px = mc_point.GetPx()
    pz = mc_point.GetPz()
    return TMath.ATan(px/pz)

def get_yz(mc_point) -> float:
    py = mc_point.GetPy()
    pz = mc_point.GetPz()
    return TMath.ATan(py/pz)


def get_intersection_z(mc_point, z) -> TVector3:
    px = mc_point.GetPx()
    py = mc_point.GetPy()
    pz = mc_point.GetPz()

    x0 = mc_point.GetX()
    y0 = mc_point.GetY()
    z0 = mc_point.GetZ()

    t = (z - z0) / pz

    x = x0 + px * t
    y = y0 + py * t

    return TVector3(x, y, z)
