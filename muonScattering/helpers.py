import scipy.stats as stats
import numpy as np
from typing import Union
from ROOT import TMath, TChain, TVector3, FairMCPoint
from sndUtils import sys, alg, system, algorithm, att
# from ddf.snd.mc import *

X0 = 13.84           # g/cm
c  = 299792458       # m/s
m  = 0.105658375523  # GeV/c²


def getPointAtZ(mcPoint, z) -> TVector3:
    px = mcPoint.GetPx()
    py = mcPoint.GetPy()
    pz = mcPoint.GetPz()

    x0 = mcPoint.GetX()
    y0 = mcPoint.GetY()
    z0 = mcPoint.GetZ()

    t = (z - z0) / pz

    x = x0 + px * t
    y = y0 + py * t

    return TVector3(x, y, z)



def isWithinAref(mcPoint,
    tt: int,
    Aref = {
        'sf': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}},
        'ds': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}}
    },
    zRef: dict = {1: 430., 11: 450., 3: 430., 13: 450.}
) -> bool:

    if not (abs(mcPoint.PdgCode())==13 and mcPoint.GetTrackID()==0):
        raise ValueError("The point is not one of a primary muon!")

    z = zRef[tt]
    x = getPointAtZ(mcPoint, z).X()
    y = getPointAtZ(mcPoint, z).Y()

    if (
        x >= Aref[sys(tt)]["min"]["x"] and
        x <= Aref[sys(tt)]["max"]["x"] and
        y >= Aref[sys(tt)]["min"]["y"] and
        y >= Aref[sys(tt)]["max"]["y"]
    ): return True
    return False


def isIP1(mcPoint, event) -> bool:
    if not event.EventHeader.isIP1():
        return False

    xz = getXZ(mcPoint)
    yz = getYZ(mcPoint)

    if (
        mcPoint.GetPz() > 0 and
        abs(xz) <= 0.02 and
        abs(yz) <= 0.02
    ): return True
    else: return False



def getXZ(mcPoint) -> float:
    px = mcPoint.GetPx()
    pz = mcPoint.GetPz()
    return TMath.ATan(px/pz)

def getYZ(mcPoint) -> float:
    py = mcPoint.GetPy()
    pz = mcPoint.GetPz()
    return TMath.ATan(py/pz)






def isUS(point) -> bool:
    detID = point.GetDetectorID()
    if detID<3e4 and detID>2e4:
        return True
    return False

def isDS(point) -> bool:
    detID = point.GetDetectorID()
    if detID>=30000:
        return True
    return False


def getPointsDS(event):
    dsPoints = []
    for dsP in event.MuFilterPoint:
        if not (
            abs(dsP.PdgCode())==13 and
            dsP.GetTrackID()==0 and
            isDS(dsP)
        ): continue
        dsPoints.append(dsP)
    return dsPoints

def getPointsSF(event):
    sfPoints = []
    for sfP in event.ScifiPoint:
        if not (
            event.EventHeader.isIP1() and
            abs(sfP.PdgCode())==13 and
            sfP.GetTrackID()==0
        ): continue
        sfPoints.append(sfP)
    return sfPoints


def getMom(point) -> float:
    px = point.GetPx()
    py = point.GetPy()
    pz = point.GetPz()
    return np.sqrt(px*px + py*py + pz*pz)

def getBeta(
    p: float,                   # [GeV]
    m: float  = 0.105658375523  # [GeV]
) -> float:
    return p/np.sqrt(m*m + p*p) # []



def getThetaRMS(
    p: float, L: float, X0: float,
    beta: Union[float, None] = None,
    z: int = 1
) -> float:
    if beta is None:
        beta = getBeta(p)

    Lz = L/X0

    K   = (0.01316) / (beta * p)
    log = np.log(Lz * (z*z) / (beta*beta))

    thetaRMS = K * z * np.sqrt(Lz) * (1 + 0.038 * log)
    return thetaRMS    # [rad]



def getProbOutAx(
    mcPoint: FairMCPoint,
    tt: int,
    Aref: dict = {
        'sf': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}},
        'ds': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}}
    },
    zRef: dict = {1: 430., 11: 450., 3: 430., 13: 450.},
    dThetaRMS: float = 0.11,
    ThetaRMS: Union[None, float] = None
):
    xi = mcPoint.GetX()
    zi = mcPoint.GetZ()

    xMin = Aref[sys(tt)]["min"]["x"]
    xMax = Aref[sys(tt)]["max"]["x"]


    xz = getXZ(mcPoint)                    # [rad]
    yz = getYZ(mcPoint)                    # [rad]

    if (
        xi > xMax or xi <  xMin or
        xz > 0.02  or xz < -0.02
    ): return 0, 0

    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
    p = getMom(mcPoint)                    # [GeV/c]
    beta = getBeta(p)                       # []

    if ThetaRMS is None:
        thetaRMS = getThetaRMS(p, dz, X0, beta)
    else:
        thetaRMS = ThetaRMS


    err = dThetaRMS * thetaRMS

    L = abs(zRef[tt] - zi)                  # [cm]
    dxMax = xMax-xi                        # [cm]
    dxMin = xMin-xi                        # [cm]

    xzMax = TMath.Tan(dxMax/L)             # [rad]
    xzMin = TMath.Tan(dxMin/L)             # [rad]

    if xzMax >  0.02: xzMax =  0.02
    if xzMin < -0.02: xzMin = -0.02

    prob = 1 - stats.norm.cdf(xzMax, loc=xz, scale=thetaRMS) + stats.norm.cdf(xzMin, loc=xz, scale=thetaRMS)

    up  = 1 - stats.norm.cdf(xzMax, loc=xz, scale=thetaRMS+err/2) + stats.norm.cdf(xzMin, loc=xz, scale=thetaRMS+err/2)
    low = 1 - stats.norm.cdf(xzMax, loc=xz, scale=thetaRMS-err/2) + stats.norm.cdf(xzMin, loc=xz, scale=thetaRMS-err/2)
    err = up - low

    return prob, err




def getProbOutAy(
    mcPoint: FairMCPoint,
    tt: int,
    Aref: dict = {
        'sf': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}},
        'ds': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}}
    },
    zRef: dict = {1: 430., 11: 450., 3: 430., 13: 450.},
    dThetaRMS: float = 0.11,
    ThetaRMS: Union[None, float] = None
):
    yi = mcPoint.GetY()
    zi = mcPoint.GetZ()

    yMin = Aref[sys(tt)]["min"]["y"]
    yMax = Aref[sys(tt)]["max"]["y"]

    xz = getXZ(mcPoint)                    # [rad]
    yz = getYZ(mcPoint)                    # [rad]

    if (
        yi > yMax or yi <  yMin or
        yz > 0.02  or xz < -0.02
    ): return 0, 0

    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
    p = getMom(mcPoint)                    # [GeV/c]
    beta = getBeta(p)                       # []

    if ThetaRMS is None:
        thetaRMS = getThetaRMS(p, dz, X0, beta)
    else:
        thetaRMS = ThetaRMS

    err = dThetaRMS * thetaRMS

    L = abs(zRef[tt] - zi)                  # [cm]
    dyMax = yMax-yi                        # [cm]
    dyMin = yMin-yi                        # [cm]

    yzMax = TMath.Tan(dyMax/L)             # [rad]
    yzMin = TMath.Tan(dyMin/L)             # [rad]

    if yzMax >  0.02:
        yzMax =  0.02
    if yzMin < -0.02:
        yzMin = -0.02

    prob = 1 - stats.norm.cdf(yzMax, loc=yz, scale=thetaRMS) + stats.norm.cdf(yzMin, loc=yz, scale=thetaRMS)
    up   = 1 - stats.norm.cdf(yzMax, loc=yz, scale=thetaRMS+err/2) + stats.norm.cdf(yzMin, loc=yz, scale=thetaRMS+err/2)
    low  = 1 - stats.norm.cdf(yzMax, loc=yz, scale=thetaRMS-err/2) + stats.norm.cdf(yzMin, loc=yz, scale=thetaRMS-err/2)
    err = up - low

    return prob, err




def getProbInAx(
    mcPoint: FairMCPoint,
    tt: int,
    Aref: dict = {
        'sf': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}},
        'ds': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}}
    },
    zRef: dict = {1: 430., 11: 450., 3: 430., 13: 450.},
    dThetaRMS: float = 0.11,
    ThetaRMS: Union[None, float] = None
):
    xi = mcPoint.GetX()
    zi = mcPoint.GetZ()

    xMin = Aref[sys(tt)]["min"]["x"]
    xMax = Aref[sys(tt)]["max"]["x"]

    xz = getXZ(mcPoint)                    # [rad]
    yz = getYZ(mcPoint)                    # [rad]

    if (
        xi <= xMax and xi >=  xMin and
        xz <= 0.02  and xz >= -0.02
    ): return 0, 0

    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
    p = getMom(mcPoint)                    # [GeV/c]
    beta = getBeta(p)                       # []

    if ThetaRMS is None:
        thetaRMS = getThetaRMS(p, dz, X0, beta)
    else:
        thetaRMS = ThetaRMS

    err = dThetaRMS * thetaRMS

    L = abs(zRef[tt] - zi)                  # [cm]
    dxMax = xMax-xi                        # [cm]
    dxMin = xMin-xi                        # [cm]

    xzMax = TMath.Tan(dxMax/L)             # [rad]
    xzMin = TMath.Tan(dxMin/L)             # [rad]

    if xzMax >  0.02: xzMax =  0.02
    if xzMin < -0.02: xzMin = -0.02

    prob = stats.norm.cdf(xzMax, loc=xz, scale=thetaRMS) - stats.norm.cdf(xzMin, loc=xz, scale=thetaRMS)
    up   = stats.norm.cdf(xzMax, loc=xz, scale=thetaRMS+err/2) - stats.norm.cdf(xzMin, loc=xz, scale=thetaRMS+err/2)
    low  = stats.norm.cdf(xzMax, loc=xz, scale=thetaRMS-err/2) - stats.norm.cdf(xzMin, loc=xz, scale=thetaRMS-err/2)
    err = up - low

    return prob, err





def getProbInAy(
    mcPoint: FairMCPoint,
    tt: int,
    Aref: dict = {
        'sf': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}},
        'ds': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}}
    },
    zRef: dict = {1: 430., 11: 450., 3: 430., 13: 450.},
    dThetaRMS: float = 0.11,
    ThetaRMS: Union[None, float] = None
):
    yi = mcPoint.GetY()
    zi = mcPoint.GetZ()

    yMin = Aref[sys(tt)]["min"]["y"]
    yMax = Aref[sys(tt)]["max"]["y"]

    xz = getXZ(mcPoint)                    # [rad]
    yz = getYZ(mcPoint)                    # [rad]

    if (
        yi <= yMax and yi >=  yMin and
        yz <= 0.02  and yz >= -0.02
    ): return 0, 0

    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
    p = getMom(mcPoint)                    # [GeV/c]
    beta = getBeta(p)                       # []

    if ThetaRMS is None:
        thetaRMS = getThetaRMS(p, dz, X0, beta)
    else:
        thetaRMS = ThetaRMS

    err = dThetaRMS * thetaRMS

    L = abs(zRef[tt] - zi)                  # [cm]
    dyMax = yMax-yi                        # [cm]
    dyMin = yMin-yi                        # [cm]

    yzMax = TMath.Tan(dyMax/L)             # [rad]
    yzMin = TMath.Tan(dyMin/L)             # [rad]

    if yzMax >  0.02: yzMax =  0.02
    if yzMin < -0.02: yzMin = -0.02

    prob = stats.norm.cdf(yzMax, loc=yz, scale=thetaRMS) - stats.norm.cdf(yzMin, loc=yz, scale=thetaRMS)
    up   = stats.norm.cdf(yzMax, loc=yz, scale=thetaRMS+err/2) - stats.norm.cdf(yzMin, loc=yz, scale=thetaRMS+err/2)
    low  = stats.norm.cdf(yzMax, loc=yz, scale=thetaRMS-err/2) - stats.norm.cdf(yzMin, loc=yz, scale=thetaRMS-err/2)
    err = up - low

    return prob, err


def getProbScatter(
    mcPoint: FairMCPoint,
    tt: int,
    Aref: dict = {
        'sf': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}},
        'ds': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}}
    },
    zRef: dict = {1: 430., 11: 450., 3: 430., 13: 450.},
    dThetaRMS: float = 0.11,
    ThetaRMS: Union[None, float] = None
):

    p_out_Ax, e_p_out_Ax = getProbOutAx(mcPoint, tt, Aref, zRef)
    p_out_Ay, e_p_out_Ay = getProbOutAy(mcPoint, tt, Aref, zRef)
    p_in_Ax,  e_p_in_Ax  = getProbInAx(mcPoint,  tt, Aref, zRef)
    p_in_Ay,  e_p_in_Ay  = getProbInAy(mcPoint,  tt, Aref, zRef)

    p_out = p_out_Ax + p_out_Ay - p_out_Ax*p_out_Ay
    e_p_out = np.sqrt(  (e_p_out_Ax + p_out_Ay*e_p_out_Ax)**2 + (e_p_out_Ay + p_out_Ax*e_p_out_Ay)**2  )

    p_in  = p_in_Ax  + p_in_Ay  - p_in_Ax*p_in_Ay
    e_p_in = np.sqrt(  (e_p_in_Ax + p_in_Ay*e_p_in_Ax)**2 + (e_p_in_Ay + p_in_Ax*e_p_in_Ay)**2  )

    prob = p_out - p_in
    err = e_p_out + e_p_in
    return prob, err




#def get_p_out_xz(mcPoint: FairMCPoint) -> float:
#    xz = getXZ(mcPoint) # [rad]
#    yz = getYZ(mcPoint) # [rad]
#
#    if xz > 0.02: return 0
#
#    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
#    p = getMom(mcPoint)                    # [GeV/c]
#    beta = getBeta(mcPoint)                # []
#
#    thetaRMS = ((13.6*1e3) / (beta*c*p)) * np.sqrt(dz/X0) * (1 + 0.038 * np.log(dz/X0))
#
#    return 1 - stats.norm.cdf(0.02, loc=xz, scale=thetaRMS) + stats.norm.cdf(-0.02, loc=xz, scale=thetaRMS)
#
#
#
#def get_p_out_yz(mcPoint: FairMCPoint) -> float:
#    xz = getXZ(mcPoint) # [rad]
#    yz = getYZ(mcPoint) # [rad]
#
#    if yz > 0.02: return 0
#
#    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
#    p = getMom(mcPoint)                    # [GeV/c]
#    beta = getBeta(mcPoint)                # []
#
#    thetaRMS = ((13.6*1e6) / (beta*c*p)) * np.sqrt(dz/X0) * (1 + 0.038 * np.log(dz/X0))
#
#    return 1 - stats.norm.cdf(20, loc=yz, scale=thetaRMS) + stats.norm.cdf(-20, loc=yz, scale=thetaRMS)
#
#
#def get_p_in_xz(mcPoint: FairMCPoint) -> float:
#    xz = getXZ(mcPoint) # [rad]
#    yz = getYZ(mcPoint) # [rad]
#
#    if xz <= 0.02: return 0
#
#    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
#    p = getMom(mcPoint)                    # [GeV/c]
#    beta = getBeta(mcPoint)                # []
#
#    thetaRMS = ((13.6*1e6) / (beta*c*p)) * np.sqrt(dz/X0) * (1 + 0.038 * np.log(dz/X0))
#
#    return stats.norm.cdf(20, loc=xz, scale=thetaRMS) - stats.norm.cdf(-20, loc=xz, scale=thetaRMS)
#
#
#
#
#def get_p_in_yz(mcPoint: FairMCPoint) -> float:
#    xz = getXZ(mcPoint) # [rad]
#    yz = getYZ(mcPoint) # [rad]
#
#    if yz <= 0.02: return 0
#
#    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
#    p = getMom(mcPoint)                    # [GeV/c]
#    beta = getBeta(mcPoint)                # []
#
#    thetaRMS = ((13.6*1e6) / (beta*c*p)) * np.sqrt(dz/X0) * (1 + 0.038 * np.log(dz/X0))
#
#    return stats.norm.cdf(20, loc=yz, scale=thetaRMS) - stats.norm.cdf(-20, loc=yz, scale=thetaRMS)
#
#
#
#def get_p_ang(mcPoint: FairMCPoint) -> float:
#    p_out_xz = get_p_out_xz(mcPoint)
#    p_out_yz = get_p_out_yz(mcPoint)
#    p_in_xz  = get_p_in_xz(mcPoint)
#    p_in_yz  = get_p_in_yz(mcPoint)
#
#    p_out = p_out_xz + p_out_yz - p_out_xz*p_out_yz
#    p_in  = p_in_xz  + p_in_yz  - p_in_xz*p_in_yz
#
#    return p_out-p_in
#
#
#def get_p_scat(
#    mcPoint: FairMCPoint,
#    tt: int,
#    Aref: dict = {
#        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
#        'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
#    },
#    zRef: dict = {1: 430., 11: 450., 3: 430., 13: 450.}
#) -> float:
#    p_out_A   = get_p_A(mcPoint, tt, Aref, zRef)
#    p_out_ang = get_p_ang(mcPoint)
#
#    return p_out_A + p_out_ang - p_out_A*p_out_ang
