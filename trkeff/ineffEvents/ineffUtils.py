import ROOT

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

def oppositeTrackExists(trk0, event):
    _att = trk0.tt
    for trk1 in event.Reco_MuonTracks:
        tt = trk1.getTrackType()
        if (
            (tt==1  and _att==3) or
            (tt==3  and _att==1) or
            (tt==11 and _att==13) or
            (tt==13 and _att==11)
        ): return True
    return False

def tracksArePaired(trk0, trk1, zRef, tt):
    zRef0 = trk0.GetPointAtZ(zRef[tt])
    x0 = zRef0.X()
    y0 = zRef0.Y()

    zRef1 = trk1.GetPointAtZ(zRef[tt])
    x1 = zRef1.X()
    y1 = zRef1.Y()

    if abs(x0-x1)<=3 and abs(y0-y1)<=3:
        return True
    else:
        return False

def getMaxSfHitsPerPlane(event):
    sfHits = event.Digi_ScifiHits
    nHits={1:0, 2:0, 3:0, 4:0, 5:0}
    for sfHit in sfHits:
        plane = sfHit.GetStation()
        nHits[plane] += 1
    return max(nHits.values())

def getMaxDsHitsPerPlane(event):
    mfHits = event.Digi_MuFilterHits
    nHits={1:0, 2:0, 3:0, 4:0}
    for dsHit in mfHits:
        if dsHit.GetSystem()!=3:
            continue

        plane = dsHit.GetPlane()+1
        nHits[plane] += 1
    return max(nHits.values())
