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
from ddfUtils import getSubDirPath, getAllFiles


def getFill(run: int, rootDir: str = "/eos/experiment/sndlhc/convertedData/physics/2023_reprocess") -> Union[int, None]:
    dataDir = getSubDirPath(TopDir=f"run_{run:06d}", RootDir=rootDir)

    try:
        file_ = getAllFiles(dataDir, "*.root")[0]
        tfile = ROOT.TFile.Open(file_)
    except Exception as e:
        print(f"Error opening file: {e}")
        return None
    if tfile.Get("cbmsim"):
        ttree = tfile.Get("cbmsim")
    elif tfile.Get("rawConv"):
        ttree = tfile.Get("rawConv")
    else:
        return None

    try:
        ttree.GetEntry(0)
        return int(ttree.EventHeader.GetFillNumber())
    except Exception as e:
        print(f"Error accessing data: {e}")
        return None

def getLumi(run: int) -> float:
    def makeUnixTime(year, month, day, hour, minute, second) :
        dt = datetime.datetime(year, month, day, hour, minute, second)
        return time.mktime(dt.timetuple())

    atlas_online_lumi = ROOT.TChain("atlas_lumi")

    fill = getFill(run)

    input_dir = "/eos/experiment/sndlhc/atlas_lumi"
    atlas_online_lumi.Add(f"{input_dir}/fill_{fill:06d}.root")

    delivered_inst_lumi = []
    delivered_unix_timestamp = []
    delivered_run_number = []
    delivered_fill_number = []
    fill = 0

    for entry in atlas_online_lumi :
        delivered_inst_lumi.append(entry.var)
        delivered_unix_timestamp.append(entry.unix_timestamp)

    recorded_mask = np.array(True)
    delivered_inst_lumi = np.array(delivered_inst_lumi)
    delivered_unix_timestamp = np.array(delivered_unix_timestamp)

    delivered_deltas = delivered_unix_timestamp[1:] - delivered_unix_timestamp[:-1]
    delivered_mask = delivered_deltas < 600

    delivered_run = np.logical_and(delivered_unix_timestamp[1:] > fill, delivered_mask)


    return np.cumsum(
        np.multiply(
            delivered_deltas[delivered_run], delivered_inst_lumi[1:][delivered_run]
        )
    )[-1]/1e3




def getEffWithErr(run: int, tt: int) -> tuple:
    i = uproot.open(f"/eos/user/i/idioniso/mfout/trkeffDataTC.root:Run{run}")[f"eff_{tt}_data.tc"]

    eff    = i["eff"].array(library="np")[0]
    effErr = i["effErr"].array(library="np")[0]

    return eff, effErr

def getFluxWithErr(run: int, tt: int, area: float):
    eff, effErr = getEffWithErr(7080, tt)

    eos="/eos/user/i/idioniso"
    mfout=f"{eos}/mfout"
    nTracks = pd.read_csv(f"{mfout}/nTracks.csv", index_col="run")

    N = nTracks.at[run, f"ntracks_{tt}_IP1"]
    dN = np.sqrt(N)

    L = getLumi(run)
    dL = 0.02 * L


    flux = N / (area * eff * L)

    dPhi_dN = 1 / (area*eff*L)
    dPhi_dL = (-N)/(area*eff*(L*L))
    dPhi_dEps = (-N)/(area*(eff*eff)*L)
    deltaPhiSquared = dPhi_dN**2 * dN**2 + dPhi_dL**2 * dL**2 + dPhi_dEps**2 * effErr**2
    fluxErr = np.sqrt(deltaPhiSquared)

    return flux, fluxErr


def saveToCsv(header, data, fout: str):
    foutExists = os.path.isfile(fout)

    with open(fout, mode='a', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)

        if not foutExists or os.stat(fout).st_size == 0:
            writer.writerow(header)

        writer.writerows(data)
    print(f"Output: {fout}")


def getFormattedData(
    inputFile: str = "/eos/user/i/idioniso/mfout/muonFlux.csv"
):
    dfMuonFlux = pd.read_csv(f"/eos/user/i/idioniso/mfout/muonFlux.csv", index_col="run")
    runs = dfMuonFlux.index.tolist()

    trackTypes = {
        1:  "Scifi simple\ntracking",
        11: "Scifi Hough\ntransform",
        3:  "DS simple\ntracking",
        13: "DS Hough\ntransform"
    }

    muonFlux = []
    for run in runs:
        for tt, TT in trackTypes.items():
            muonFlux.append((TT, f"Run {run}", dfMuonFlux.at[run, f"flux_{tt}"], dfMuonFlux.at[run, f"fluxErr_{tt}"]))

    columns=['Track type', 'Dataset', 'Muon flux', 'Errors']
    return pd.DataFrame(muonFlux, columns=columns)








def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run', type=int, default=7080)
    args = parser.parse_args()
    run = args.run

    L = getLumi(run)
    dL = 0.02 * L

    eos="/eos/user/i/idioniso"
    mfout=f"{eos}/mfout"

    A = {
        'min': {'x': -42., 'y': 19.},
        'max': {'x': -10., 'y': 48.}
    }
    area = abs(A["max"]["x"] - A["min"]["x"]) * abs(A["max"]["y"] - A["min"]["y"])

    nTracks = pd.read_csv(f"{mfout}/nTracks.csv", index_col="run")

    eff = {}; N = {}; flux = {}
    header = ["run"]
    data = [[run]]
    for tt in (1, 11, 3, 13):
        eff[tt] = {}
        eff[tt]["v"], eff[tt]["e"] = getEffWithErr(7080, tt)

        N[tt] = {}
        N[tt]["v"] = nTracks.at[run, f"ntracks_{tt}_IP1"]
        N[tt]["e"] = np.sqrt(N[tt]["v"])

        flux[tt] = {}
        flux[tt]["v"], flux[tt]["e"] = getFluxWithErr(run, tt, area)

        print(f">> {tt}:\t{flux[tt]['v']/1e4:.03f} ± {flux[tt]['e']/1e4:.03f}    \033[1;39m×10⁴ [nb/cm²]\033[0m")


        header.append(f"flux_{tt}")
        data[0].append(flux[tt]["v"])

    for tt in (1, 11, 3, 13):
        header.append(f"fluxErr_{tt}")
        data[0].append(flux[tt]["e"])

    saveToCsv(header, data, f"{mfout}/muonFlux.csv")



if __name__=="__main__":
    main()
