import ROOT
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import datetime
import dateutil.parser as dp
import time
import os, csv, json
import uproot
from typing import Union
from sndUtils import getFill, getLumi, sys, alg
from ddfUtils import getSubDirPath, getAllFiles



def getEffWithErr(run: int, tt: int) -> tuple:
    fIn = uproot.open("/eos/user/i/idioniso/mfout/MuonFlux.root")
    df = fIn["TrackingEfficiency/trkeff"].arrays(library="pd")


    eff    = df.query(f"run == {run}")[f"trkeff_{tt}"].values[0]
    effErr = df.query(f"run == {run}")[f"trkeffErr_{tt}"].values[0]

    return eff, effErr

def getFluxWithErr(run: int, tt: int, area: float, lumi: str = "eos"):
    eff, effErr = getEffWithErr(7080, tt)

    eos="/eos/user/i/idioniso"
    mfout=f"{eos}/mfout"

    fIn = uproot.open("/eos/user/i/idioniso/mfout/MuonFlux.root")
    df = fIn["nTracks"].arrays(library="pd")

    N     = df.query(f"run == {run}")[f"ntracks_{tt}_IP1"].values[0]
    scale = df.query(f"run == {run}")[f"scale"].values[0]
    dN    = np.sqrt(N)

    L = getLumi(run, option=lumi)
    if tt==1 and run==7080:
        print(f"Run {run} , Lumi [1/ub] ({lumi}) = {L*1e3}")
    dL = 0.02 * L


    flux = (N * scale) / (area * eff * L)

    dPhi_dN = 1 / (area*eff*L)
    dPhi_dL = (-N)/(area*eff*(L*L))
    dPhi_dEps = (-N)/(area*(eff*eff)*L)
    deltaPhiSquared = dPhi_dN**2 * dN**2 + dPhi_dL**2 * dL**2 + dPhi_dEps**2 * effErr**2
    fluxErr = scale * np.sqrt(deltaPhiSquared)

    if tt==1 and run==7080:
        print(f"Flux [nb/cm^2] ({lumi}) = {flux}")
    return flux, fluxErr



def updateNtracks(rootFile: str = "/eos/user/i/idioniso/mfout/MuonFlux.root"):
    df = pd.read_csv("/eos/user/i/idioniso/mfout/nTracks.csv")

    with uproot.update(rootFile) as f:
        f["nTracks"] = df


def main():
    eos="/eos/user/i/idioniso"
    mfout=f"{eos}/mfout"

    A = {
        'min': {'x': -42., 'y': 19.},
        'max': {'x': -10., 'y': 48.}
    }
    area = abs(A["max"]["x"] - A["min"]["x"]) * abs(A["max"]["y"] - A["min"]["y"])


    with uproot.open("/eos/user/i/idioniso/mfout/MuonFlux.root") as fRead:
        nTracks = fRead["nTracks"].arrays(library="pd")

    with uproot.update("/eos/user/i/idioniso/mfout/MuonFlux.root") as fUpdate:
        nTracks = nTracks.loc[nTracks["run"].apply(getFill) != 9301]
        dfMuonFluxEos = nTracks.copy()
        dfMuonFluxSt  = nTracks.copy()

        for tt in (1, 11, 3, 13):
            dfMuonFluxEos[f"flux_{tt}"], dfMuonFluxEos[f"fluxErr_{tt}"] = zip(
                *dfMuonFluxEos["run"].apply(
                    lambda r: getFluxWithErr(r, tt, area, lumi="eos")
                )
            )

            dfMuonFluxSt[f"flux_{tt}"], dfMuonFluxSt[f"fluxErr_{tt}"] = zip(
                *dfMuonFluxSt["run"].apply(
                    lambda r: getFluxWithErr(r, tt, area, lumi="st")
                )
            )

        dfMuonFluxEos.to_csv(f"~/MuonFlux.csv")
        dfMuonFluxSt.to_csv(f"~/MuonFluxSt.csv")

        fUpdate["MuonFlux"] = dfMuonFluxEos
        fUpdate["MuonFluxSt"] = dfMuonFluxSt


if __name__=="__main__":
    main()
