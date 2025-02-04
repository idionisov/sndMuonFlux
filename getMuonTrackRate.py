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



def getRunDuration(
    run: int,
):
    tchain = ROOT.TChain("cbmsim")
    tchain.Add(f"/eos/experiment/sndlhc/convertedData/physics/2023_reprocess/run_{run:06d}/*.root")

    tchain.GetEntry(0)
    tStart = tchain.EventHeader.GetUTCtimestamp()
    tchain.GetEntry(tchain.GetEntries()-1)
    tEnd = tchain.EventHeader.GetUTCtimestamp()

    return tEnd-tStart

def getEffWithErr(run: int, tt: int) -> tuple:
    i = uproot.open(f"/eos/user/i/idioniso/mfout/trkeffDataTC.root:Run{run}")[f"eff_{tt}_data.tc"]

    eff    = i["eff"].array(library="np")[0]
    effErr = i["effErr"].array(library="np")[0]

    return eff, effErr

def getTrkRateWithErr(run: int, tt: int, area: float) -> tuple:
    eos="/eos/user/i/idioniso"
    mfout=f"{eos}/mfout"
    nTracks = pd.read_csv(f"{mfout}/nTracks.csv", index_col="run")

    N = nTracks.at[run, f"ntracks_{tt}_all"]
    dN = np.sqrt(N)

    T = getRunDuration(run)
    dT = 0

    rate = (N * 961) / (area * T)

    dR_dN = 961 / (area*T)
    dR_dT = (N*961) / (area * T**2)
    dR_squared = (dR_dN * dN)**2 + (dR_dT * dT)**2
    rateErr = np.sqrt(dR_squared)

    return rate, rateErr


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
    dfMuonTrkRate = pd.read_csv(f"/eos/user/i/idioniso/mfout/muonTrkRate.csv", index_col="run")
    runs = dfMuonTrkRate.index.tolist()

    trackTypes = {
        1:  "Scifi simple\ntracking",
        11: "Scifi Hough\ntransform",
        3:  "DS simple\ntracking",
        13: "DS Hough\ntransform"
    }

    muonTrkRate = []
    for run in runs:
        for tt, TT in trackTypes.items():
            muonTrkRate.append((TT, f"Run {run}", dfMuonTrkRate.at[run, f"flux_{tt}"], dfMuonTrkRate.at[run, f"fluxErr_{tt}"]))

    columns=['Track type', 'Dataset', 'Muon track rate', 'Errors']
    return pd.DataFrame(muonTrkRate, columns=columns)






def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run', type=int, default=7080)
    args = parser.parse_args()
    run = args.run

    T = getRunDuration(run)
    dT = 0.02 * T

    eos="/eos/user/i/idioniso"
    mfout=f"{eos}/mfout"

    A = {
        'min': {'x': -42., 'y': 19.},
        'max': {'x': -10., 'y': 48.}
    }
    area = abs(A["max"]["x"] - A["min"]["x"]) * abs(A["max"]["y"] - A["min"]["y"])

    nTracks = pd.read_csv(f"{mfout}/nTracks.csv", index_col="run")

    N = {}; rate = {}
    header = ["run"]
    data = [[run]]
    for tt in (1, 11, 3, 13):
        N[tt] = {}
        N[tt]["v"] = nTracks.at[run, f"ntracks_{tt}_all"]
        N[tt]["e"] = np.sqrt(N[tt]["v"])

        rate[tt] = {}
        rate[tt]["v"], rate[tt]["e"] = getTrkRateWithErr(run, tt, area)

        print(f">> {tt}:\t{rate[tt]['v']:.03f} ± {rate[tt]['e']:.03f}    \033[1;39m [Hz]\033[0m")


        header.append(f"trkRate_{tt}")
        data[0].append(rate[tt]["v"])

    for tt in (1, 11, 3, 13):
        header.append(f"trkRateErr_{tt}")
        data[0].append(rate[tt]["e"])

    saveToCsv(header, data, f"{mfout}/muonTrkRate.csv")



if __name__=="__main__":
    main()
