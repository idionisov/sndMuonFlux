import scipy.stats as stats
import numpy as np
from typing import Union
from ROOT import TMath, TChain, TVector3, FairMCPoint
from ddf.snd.trk import sys, alg, sys_name, alg_name, get_anti_tt
from ddf.snd.mc import *

X0 = 13.84           # g/cm
c  = 299792458       # m/s
m  = 0.105658375523  # GeV/c²


def is_within_Aref(
    mc_point, tt: int,
    Aref: dict = {
        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
        'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
    },
    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.}
) -> bool:
    
    if not (abs(mc_point.PdgCode())==13 and mc_point.GetTrackID()==0):
        raise ValueError("The point is not one of a primary muon!")

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


def get_mom(point) -> float:
    px = point.GetPx()
    py = point.GetPy()
    pz = point.GetPz()
    return np.sqrt(px*px + py*py + pz*pz)
    
def get_beta(
    p: float,                   # [GeV]
    m: float  = 0.105658375523  # [GeV]
) -> float:
    return p/np.sqrt(m*m + p*p) # []



def get_theta_rms(
    p: float, L: float, X0: float,
    beta: Union[float, None] = None,
    z: int = 1
) -> float:
    if beta is None: beta = get_beta(p)

    Lz = L/X0

    K   = (0.01316) / (beta * p)
    log = np.log(Lz * (z*z) / (beta*beta))

    theta_rms = K * z * np.sqrt(Lz) * (1 + 0.038 * log)
    return theta_rms    # [rad]



def get_p_out_Ax(
    mc_point: FairMCPoint,
    tt: int,
    Aref: dict = {
        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
        'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
    },
    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.},
    dThetaRMS: float = 0.11,
    ThetaRMS: Union[None, float] = None
):
    xi = mc_point.GetX()
    zi = mc_point.GetZ()

    x_min = Aref[sys(tt)]["min"]["x"]
    x_max = Aref[sys(tt)]["max"]["x"]
   
    
    xz = get_xz(mc_point)                    # [rad]
    yz = get_yz(mc_point)                    # [rad]

    if (
        xi > x_max or xi <  x_min or
        xz > 0.02  or xz < -0.02
    ): return 0, 0

    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
    p = get_mom(mc_point)                    # [GeV/c]
    beta = get_beta(p)                       # []

    if ThetaRMS is None:
        theta_rms = get_theta_rms(p, dz, X0, beta)
    else:
        theta_rms = ThetaRMS


    err = dThetaRMS * theta_rms

    L = abs(z_ref[tt] - zi)                  # [cm]
    dx_max = x_max-xi                        # [cm]
    dx_min = x_min-xi                        # [cm]

    xz_max = TMath.Tan(dx_max/L)             # [rad]
    xz_min = TMath.Tan(dx_min/L)             # [rad]

    if xz_max >  0.02: xz_max =  0.02,
    if xz_min < -0.02: xz_min = -0.02

    prob = 1 - stats.norm.cdf(xz_max, loc=xz, scale=theta_rms) + stats.norm.cdf(xz_min, loc=xz, scale=theta_rms) 
    
    up  = 1 - stats.norm.cdf(xz_max, loc=xz, scale=theta_rms+err/2) + stats.norm.cdf(xz_min, loc=xz, scale=theta_rms+err/2)
    low = 1 - stats.norm.cdf(xz_max, loc=xz, scale=theta_rms-err/2) + stats.norm.cdf(xz_min, loc=xz, scale=theta_rms-err/2)
    err = up - low

    return prob, err




def get_p_out_Ay(
    mc_point: FairMCPoint,
    tt: int,
    Aref: dict = {
        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
        'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
    },
    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.},
    dThetaRMS: float = 0.11,
    ThetaRMS: Union[None, float] = None
):
    yi = mc_point.GetY()
    zi = mc_point.GetZ()

    y_min = Aref[sys(tt)]["min"]["y"]
    y_max = Aref[sys(tt)]["max"]["y"]
   
    xz = get_xz(mc_point)                    # [rad]
    yz = get_yz(mc_point)                    # [rad]

    if (
        yi > y_max or yi <  y_min or
        yz > 0.02  or xz < -0.02
    ): return 0, 0

    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
    p = get_mom(mc_point)                    # [GeV/c]
    beta = get_beta(p)                       # []

    if ThetaRMS is None:
        theta_rms = get_theta_rms(p, dz, X0, beta)
    else:
        theta_rms = ThetaRMS

    err = dThetaRMS * theta_rms

    L = abs(z_ref[tt] - zi)                  # [cm]
    dy_max = y_max-yi                        # [cm]
    dy_min = y_min-yi                        # [cm]

    yz_max = TMath.Tan(dy_max/L)             # [rad]
    yz_min = TMath.Tan(dy_min/L)             # [rad]

    if yz_max >  0.02: yz_max =  0.02
    if yz_min < -0.02: yz_min = -0.02

    prob = 1 - stats.norm.cdf(yz_max, loc=yz, scale=theta_rms) + stats.norm.cdf(yz_min, loc=yz, scale=theta_rms) 
    up   = 1 - stats.norm.cdf(yz_max, loc=yz, scale=theta_rms+err/2) + stats.norm.cdf(yz_min, loc=yz, scale=theta_rms+err/2) 
    low  = 1 - stats.norm.cdf(yz_max, loc=yz, scale=theta_rms-err/2) + stats.norm.cdf(yz_min, loc=yz, scale=theta_rms-err/2) 
    err = up - low

    return prob, err




def get_p_in_Ax(
    mc_point: FairMCPoint,
    tt: int,
    Aref: dict = {
        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
        'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
    },
    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.},
    dThetaRMS: float = 0.11,
    ThetaRMS: Union[None, float] = None
):
    xi = mc_point.GetX()
    zi = mc_point.GetZ()

    x_min = Aref[sys(tt)]["min"]["x"]
    x_max = Aref[sys(tt)]["max"]["x"]
   
    xz = get_xz(mc_point)                    # [rad]
    yz = get_yz(mc_point)                    # [rad]

    if (
        xi <= x_max and xi >=  x_min and
        xz <= 0.02  and xz >= -0.02
    ): return 0, 0

    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
    p = get_mom(mc_point)                    # [GeV/c]
    beta = get_beta(p)                       # []

    if ThetaRMS is None:
        theta_rms = get_theta_rms(p, dz, X0, beta)
    else:
        theta_rms = ThetaRMS

    err = dThetaRMS * theta_rms

    L = abs(z_ref[tt] - zi)                  # [cm]
    dx_max = x_max-xi                        # [cm]
    dx_min = x_min-xi                        # [cm]

    xz_max = TMath.Tan(dx_max/L)             # [rad]
    xz_min = TMath.Tan(dx_min/L)             # [rad]

    if xz_max >  0.02: xz_max =  0.02
    if xz_min < -0.02: xz_min = -0.02

    prob = stats.norm.cdf(xz_max, loc=xz, scale=theta_rms) - stats.norm.cdf(xz_min, loc=xz, scale=theta_rms)
    up   = stats.norm.cdf(xz_max, loc=xz, scale=theta_rms+err/2) - stats.norm.cdf(xz_min, loc=xz, scale=theta_rms+err/2)
    low  = stats.norm.cdf(xz_max, loc=xz, scale=theta_rms-err/2) - stats.norm.cdf(xz_min, loc=xz, scale=theta_rms-err/2)
    err = up - low 
 
    return prob, err





def get_p_in_Ay(
    mc_point: FairMCPoint,
    tt: int,
    Aref: dict = {
        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
        'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
    },
    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.},
    dThetaRMS: float = 0.11,
    ThetaRMS: Union[None, float] = None
):
    yi = mc_point.GetY()
    zi = mc_point.GetZ()

    y_min = Aref[sys(tt)]["min"]["y"]
    y_max = Aref[sys(tt)]["max"]["y"]
   
    xz = get_xz(mc_point)                    # [rad]
    yz = get_yz(mc_point)                    # [rad]

    if (
        yi <= y_max and yi >=  y_min and
        yz <= 0.02  and yz >= -0.02
    ): return 0, 0

    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
    p = get_mom(mc_point)                    # [GeV/c]
    beta = get_beta(p)                       # []

    if ThetaRMS is None:
        theta_rms = get_theta_rms(p, dz, X0, beta)
    else:
        theta_rms = ThetaRMS

    err = dThetaRMS * theta_rms

    L = abs(z_ref[tt] - zi)                  # [cm]
    dy_max = y_max-yi                        # [cm]
    dy_min = y_min-yi                        # [cm]

    yz_max = TMath.Tan(dy_max/L)             # [rad]
    yz_min = TMath.Tan(dy_min/L)             # [rad]

    if yz_max >  0.02: yz_max =  0.02
    if yz_min < -0.02: yz_min = -0.02

    prob = stats.norm.cdf(yz_max, loc=yz, scale=theta_rms) - stats.norm.cdf(yz_min, loc=yz, scale=theta_rms)
    up   = stats.norm.cdf(yz_max, loc=yz, scale=theta_rms+err/2) - stats.norm.cdf(yz_min, loc=yz, scale=theta_rms+err/2)
    low  = stats.norm.cdf(yz_max, loc=yz, scale=theta_rms-err/2) - stats.norm.cdf(yz_min, loc=yz, scale=theta_rms-err/2)
    err = up - low

    return prob, err


def get_p_scat(
    mc_point: FairMCPoint,
    tt: int,
    Aref: dict = {
        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
        'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
    },
    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.},
    dThetaRMS: float = 0.11,
    ThetaRMS: Union[None, float] = None
):
    
    p_out_Ax, e_p_out_Ax = get_p_out_Ax(mc_point, tt, Aref, z_ref)
    p_out_Ay, e_p_out_Ay = get_p_out_Ay(mc_point, tt, Aref, z_ref)
    p_in_Ax,  e_p_in_Ax  = get_p_in_Ax(mc_point,  tt, Aref, z_ref)
    p_in_Ay,  e_p_in_Ay  = get_p_in_Ay(mc_point,  tt, Aref, z_ref)

    p_out = p_out_Ax + p_out_Ay - p_out_Ax*p_out_Ay
    e_p_out = np.sqrt(  (e_p_out_Ax + p_out_Ay*e_p_out_Ax)**2 + (e_p_out_Ay + p_out_Ax*e_p_out_Ay)**2  )

    p_in  = p_in_Ax  + p_in_Ay  - p_in_Ax*p_in_Ay
    e_p_in = np.sqrt(  (e_p_in_Ax + p_in_Ay*e_p_in_Ax)**2 + (e_p_in_Ay + p_in_Ax*e_p_in_Ay)**2  )

    prob = p_out - p_in
    err = e_p_out + e_p_in
    return prob, err




#def get_p_out_xz(mc_point: FairMCPoint) -> float:
#    xz = get_xz(mc_point) # [rad]
#    yz = get_yz(mc_point) # [rad]
#
#    if xz > 0.02: return 0
#
#    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
#    p = get_mom(mc_point)                    # [GeV/c]
#    beta = get_beta(mc_point)                # []
#
#    theta_rms = ((13.6*1e3) / (beta*c*p)) * np.sqrt(dz/X0) * (1 + 0.038 * np.log(dz/X0))
#
#    return 1 - stats.norm.cdf(0.02, loc=xz, scale=theta_rms) + stats.norm.cdf(-0.02, loc=xz, scale=theta_rms)
#
#
#
#def get_p_out_yz(mc_point: FairMCPoint) -> float:
#    xz = get_xz(mc_point) # [rad]
#    yz = get_yz(mc_point) # [rad]
#
#    if yz > 0.02: return 0
#
#    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
#    p = get_mom(mc_point)                    # [GeV/c]
#    beta = get_beta(mc_point)                # []
#
#    theta_rms = ((13.6*1e6) / (beta*c*p)) * np.sqrt(dz/X0) * (1 + 0.038 * np.log(dz/X0))
#
#    return 1 - stats.norm.cdf(20, loc=yz, scale=theta_rms) + stats.norm.cdf(-20, loc=yz, scale=theta_rms)
#
#
#def get_p_in_xz(mc_point: FairMCPoint) -> float:
#    xz = get_xz(mc_point) # [rad]
#    yz = get_yz(mc_point) # [rad]
#
#    if xz <= 0.02: return 0
#
#    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
#    p = get_mom(mc_point)                    # [GeV/c]
#    beta = get_beta(mc_point)                # []
#
#    theta_rms = ((13.6*1e6) / (beta*c*p)) * np.sqrt(dz/X0) * (1 + 0.038 * np.log(dz/X0))
#
#    return stats.norm.cdf(20, loc=xz, scale=theta_rms) - stats.norm.cdf(-20, loc=xz, scale=theta_rms)
#
#
#
#
#def get_p_in_yz(mc_point: FairMCPoint) -> float:
#    xz = get_xz(mc_point) # [rad]
#    yz = get_yz(mc_point) # [rad]
#
#    if yz <= 0.02: return 0
#
#    dz = 160/(TMath.Cos(xz) * TMath.Cos(yz)) # [cm]
#    p = get_mom(mc_point)                    # [GeV/c]
#    beta = get_beta(mc_point)                # []
#
#    theta_rms = ((13.6*1e6) / (beta*c*p)) * np.sqrt(dz/X0) * (1 + 0.038 * np.log(dz/X0))
#
#    return stats.norm.cdf(20, loc=yz, scale=theta_rms) - stats.norm.cdf(-20, loc=yz, scale=theta_rms)
#
#
#
#def get_p_ang(mc_point: FairMCPoint) -> float:
#    p_out_xz = get_p_out_xz(mc_point)
#    p_out_yz = get_p_out_yz(mc_point)
#    p_in_xz  = get_p_in_xz(mc_point)
#    p_in_yz  = get_p_in_yz(mc_point)
#
#    p_out = p_out_xz + p_out_yz - p_out_xz*p_out_yz
#    p_in  = p_in_xz  + p_in_yz  - p_in_xz*p_in_yz
#
#    return p_out-p_in
#
#
#def get_p_scat(
#    mc_point: FairMCPoint,
#    tt: int,
#    Aref: dict = {
#        'sf': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}},
#        'ds': {'min': {'x': -42., 'y': 18.}, 'max': {'x': -11., 'y': 49.}}
#    },
#    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.}
#) -> float:
#    p_out_A   = get_p_A(mc_point, tt, Aref, z_ref)
#    p_out_ang = get_p_ang(mc_point)
#
#    return p_out_A + p_out_ang - p_out_A*p_out_ang
