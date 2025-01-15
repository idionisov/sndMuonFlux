import ROOT
from typing import Union
from sndUtils import DdfTrack, SndData

def isSelected_1(
    tagTrack: DdfTrack,
    sndData: Union[SndData, ROOT.MuFilter],
    eventOrMfHits: Union[ROOT.TChain, ROOT.TClonesArray, None] = None,
    xz_min: float = -20.,
    xz_max: float =  20.,
    yz_min: float = -20.,
    yz_max: float =  20.
) -> bool:
    if eventOrMfHits is None:
        if tagTrack.Event is not None:
            mfHits = tagTrack.Event.Digi_MuFilterHits
        else:
            raise ValueError("Neither an event nor MuFilterHits are provided!")
    else:
        mfHits = eventOrMfHits
    del eventOrMfHits

    if isinstance(sndData, ROOT.MuFilter):
        mufi = sndData
    else:
        mufi = sndData.Mufi

    if not (
        tagTrack.tt == 1 and
        tagTrack.IsGood(xz_min=xz_min, xz_max=xz_max, yz_min=yz_min, yz_max=yz_max) and
        tagTrack.IsWithinUS5Bar(mufi, mfHits) and
        tagTrack.IsWithinDS3()
    ): return False
    else: return True

def isSelected_11(
    tagTrack: DdfTrack,
    sndData: Union[SndData, ROOT.MuFilter],
    eventOrMfHits: Union[ROOT.TChain, ROOT.TClonesArray, None] = None,
    xz_min: float = -20.,
    xz_max: float =  20.,
    yz_min: float = -20.,
    yz_max: float =  20.
) -> bool:
    if eventOrMfHits is None:
        if tagTrack.Event is not None:
            mfHits = tagTrack.Event.Digi_MuFilterHits
        else:
            raise ValueError("Neither an event nor MuFilterHits are provided!")
    else:
        mfHits = eventOrMfHits
    del eventOrMfHits

    if isinstance(sndData, ROOT.MuFilter):
        mufi = sndData
    else:
        mufi = sndData.Mufi

    if not (
        tagTrack.tt == 11 and
        tagTrack.IsGood(xz_min=xz_min, xz_max=xz_max, yz_min=yz_min, yz_max=yz_max) and
        tagTrack.IsWithinUS5Bar(mufi, mfHits) and
        tagTrack.IsWithinDS3()
    ): return False
    else: return True

def isSelected_3(
    tagTrack: DdfTrack,
    sndData: Union[SndData, ROOT.MuFilter],
    eventOrMfHits: Union[ROOT.TChain, ROOT.TClonesArray, None] = None,
    xz_min: float = -20.,
    xz_max: float =  20.,
    yz_min: float = -20.,
    yz_max: float =  20.
) -> bool:
    if eventOrMfHits is None:
        if tagTrack.Event is not None:
            mfHits = tagTrack.Event.Digi_MuFilterHits
        else:
            raise ValueError("Neither an event nor MuFilterHits are provided!")
    else:
        mfHits = eventOrMfHits
    del eventOrMfHits


    if isinstance(sndData, ROOT.MuFilter):
        mufi = sndData
    else:
        mufi = sndData.Mufi

    if not (
        tagTrack.tt == 3 and
        tagTrack.IsGood(xz_min=xz_min, xz_max=xz_max, yz_min=yz_min, yz_max=yz_max) and
        tagTrack.IsWithinVetoBar(mufi, mfHits)
    ): return False
    else: return True


def isSelected_13(
    tagTrack: DdfTrack,
    sndData: Union[SndData, ROOT.MuFilter],
    eventOrMfHits: Union[ROOT.TChain, ROOT.TClonesArray, None] = None,
    xz_min: float = -20.,
    xz_max: float =  20.,
    yz_min: float = -20.,
    yz_max: float =  20.
) -> bool:
    if eventOrMfHits is None:
        if tagTrack.Event is not None:
            mfHits = tagTrack.Event.Digi_MuFilterHits
        else:
            raise ValueError("Neither an event nor MuFilterHits are provided!")
    else:
        mfHits = eventOrMfHits
    del eventOrMfHits


    if isinstance(sndData, ROOT.MuFilter):
        mufi = sndData
    else:
        mufi = sndData.Mufi

    if not (
        tagTrack.tt == 13 and
        tagTrack.IsGood(xz_min=xz_min, xz_max=xz_max, yz_min=yz_min, yz_max=yz_max) and
        tagTrack.IsWithinVetoBar(mufi, mfHits)
    ): return False
    else: return True



def isSelected(
    tt: int,
    tagTrack: DdfTrack,
    sndData: Union[SndData, ROOT.MuFilter],
    eventOrMfHits: Union[ROOT.TChain, ROOT.TClonesArray, None] = None,
    xz_min: float = -20.,
    xz_max: float =  20.,
    yz_min: float = -20.,
    yz_max: float =  20.
) -> bool:
    if tagTrack.att() != tt:
        return False
    if eventOrMfHits is None:
        if tagTrack.Event is not None:
            mfHits = tagTrack.Event.Digi_MuFilterHits
        else:
            raise ValueError("Neither an event nor MuFilterHits are provided!")
    else:
        mfHits = eventOrMfHits
    del eventOrMfHits

    if isinstance(sndData, ROOT.MuFilter):
        mufi = sndData
    else:
        mufi = sndData.Mufi

    if tt==11:
        return isSelected_13(tagTrack, mufi, mfHits,
            xz_min=xz_min/1e3,
            xz_max=xz_max/1e3,
            yz_min=yz_min/1e3,
            yz_max=yz_max/1e3
        )

    elif tt==1:
        return isSelected_3(tagTrack, mufi, mfHits,
            xz_min=xz_min/1e3,
            xz_max=xz_max/1e3,
            yz_min=yz_min/1e3,
            yz_max=yz_max/1e3
        )

    elif tt==13:
        return isSelected_11(tagTrack, mufi, mfHits,
            xz_min=xz_min/1e3,
            xz_max=xz_max/1e3,
            yz_min=yz_min/1e3,
            yz_max=yz_max/1e3
        )

    elif tt==3:
        return isSelected_1(tagTrack, mufi, mfHits,
            xz_min=xz_min/1e3,
            xz_max=xz_max/1e3,
            yz_min=yz_min/1e3,
            yz_max=yz_max/1e3
        )

    else: return False
