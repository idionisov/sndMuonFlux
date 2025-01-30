import uproot
import os, csv
import numpy as np
from ddfRoot import getNumpyFromTH2, getNumpyFromUprootTH2

def saveToCsv(header, data, fout: str):
    foutExists = os.path.isfile(fout)

    with open(fout, mode='a', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)

        if not foutExists or os.stat(fout).st_size == 0:
            writer.writerow(header)

        writer.writerows(data)
    print(f"Output: {fout}")


def saveVarChi2ndfToCsv(
    fout: str = "/eos/user/i/idioniso/mfout/trkeff/trkeff-chi2ndfMax/trkeff-chi2ndf.csv",
    inputDir: str = "/eos/user/i/idioniso/mfout/trkeff/trkeff-chi2ndfMax",
    chi2ndf: dict = {
        1:  [2.5, 5, 7.5, 10, 15, 20, 30, 40],
        11: [2, 3, 4, 5, 7.5, 10, 15, 20],
        3:  [2.5, 5, 10, 15, 30, 50, 80, 120],
        13: [5, 7.5, 10, 25, 50, 80, 125, 250]
    }
):
    TTs = (1, 11, 3, 13)
    mfout = "/eos/user/i/idioniso/mfout"
    header, data = [], []

    for tt in TTs: header.append(f"chi2ndf_{tt}")
    for tt in TTs: header.append(f"trkeff_{tt}")
    for tt in TTs: header.append(f"trkeffErr_{tt}")

    for i_chi, chi in enumerate(chi2ndf[1]):
        data.append([])
        for tt in TTs:
            data[-1].append(chi2ndf[tt][i_chi])

        with uproot.open(f"{mfout}/trkeff/trkeff-chi2ndfMax/trkeff_chi2ndf.{chi}_tc.root") as fin:
            eff, effErr = {}, {}

            for tt in TTs:
                eff[tt]    = fin[f"eff_{tt}_data.tc/eff"].array()[0]
                effErr[tt] = fin[f"eff_{tt}_data.tc/effErr"].array()[0]

                data[-1].append(eff[tt])

            for tt in TTs:
                data[-1].append(effErr[tt])

    saveToCsv(header, data, fout)


def saveTrkeffAsCsv(
    fout: str = "/eos/user/i/idioniso/mfout/trkeff.csv",
    input: str = "/eos/user/i/idioniso/mfout/trkeff.root",
    runs: tuple = (7080, 7268, 7329, 7346)
):
    if not fout.endswith(".csv"):
        fout = f"{fout}.csv"

    if not input.endswith(".root"):
        input = f"{input}.root"

    with uproot.open(input) as fin:
        runs = (7080, 7268, 7329, 7346)
        header = ["run"]
        data = []

        eff = {}
        effErr = {}
        for tt in (1, 11, 3, 13):
            eff[tt] = {}
            effErr[tt] = {}

            header.append(f"trkeff_{tt}")
            header.append(f"trkeffErr_{tt}")

        for i_run, run in enumerate(runs):
            data.append([run])
            for i_tt, tt in enumerate((1, 11, 3, 13)):

                eff[tt][run] = fin[f"Run{run}/eff_{tt}_data.tc/eff"].array()[0]
                effErr[tt][run] = fin[f"Run{run}/eff_{tt}_data.tc/effErr"].array()[0]

                data[i_run].append(eff[tt][run])
                data[i_run].append(effErr[tt][run])

        data.append([-2])
        effErrUp, effErrLow = {}, {}
        for i_tt, tt in enumerate((1, 11, 3, 13)):
            effErrUp[tt], effErrLow[tt] = {}, {}

            eff[tt][-2] = fin[f"MonteCarlo-RT/eff_{tt}_sim.rt"]["eff"].array()[0]
            effErrUp[tt][-2] = fin[f"MonteCarlo-RT/eff_{tt}_sim.rt/effErrUp"].array()[0]
            effErrLow[tt][-2] = fin[f"MonteCarlo-RT/eff_{tt}_sim.rt/effErrLow"].array()[0]

            data[len(runs)].append(eff[tt][-1])
            data[len(runs)].append(np.mean([effErrLow[tt][-1], effErrUp[tt][-1]]))


        saveToCsv(header, data, fout)




def saveVarEToCsv(
    fout: str = "/eos/user/i/idioniso/mfout/trkeff/trkeff-E/trkeff-E.csv",
    inputDir: str = "/eos/user/i/idioniso/mfout/trkeff/trkeff-E",
    energies: list = [10, 12.5, 15, 17.5, 20, 25, 30, 35, 55, 100, 200, 300, 600, 800, 1010]
):
    TTs = (1, 11, 3, 13)
    mfout = "/eos/user/i/idioniso/mfout"
    header, data = [], []

    header.append(f"e")
    for tt in TTs: header.append(f"trkeff_{tt}_tc")
    for tt in TTs: header.append(f"trkeffErr_{tt}_tc")
    for tt in TTs: header.append(f"trkeff_{tt}_rt")
    for tt in TTs: header.append(f"trkeffErr_{tt}_rt")

    for i_e, e in enumerate(energies):
        data.append([e])

        with uproot.open(f"{inputDir}/trkeff_muGun.{e}GeV_tc.root") as fin:
            eff, effErrLow, effErrUp, effErr = {}, {}, {}, {}

            for tt in TTs:
                eff[tt]       = fin[f"eff_{tt}_muGun.tc/eff"].array()[0]
                effErrLow[tt] = fin[f"eff_{tt}_muGun.tc/effErrLow"].array()[0]
                effErrUp[tt]  = fin[f"eff_{tt}_muGun.tc/effErrUp"].array()[0]
                effErr[tt]    = np.mean([effErrLow[tt], effErrUp[tt]])

            for tt in TTs: data[-1].append(eff[tt])
            for tt in TTs: data[-1].append(effErr[tt])


        with uproot.open(f"{inputDir}/trkeff_muGun.{e}GeV_rt.root") as fin:
            eff, effErrLow, effErrUp, effErr = {}, {}, {}, {}

            for tt in TTs:
                eff[tt]       = fin[f"eff_{tt}_muGun.rt/eff"].array()[0]
                effErrLow[tt] = fin[f"eff_{tt}_muGun.rt/effErrLow"].array()[0]
                effErrUp[tt]  = fin[f"eff_{tt}_muGun.rt/effErrUp"].array()[0]
                effErr[tt]    = np.mean([effErrLow[tt], effErrUp[tt]])

            for tt in TTs: data[-1].append(eff[tt])
            for tt in TTs: data[-1].append(effErr[tt])


    saveToCsv(header, data, fout)



def saveVarAngToCsv(
    fout: str = "/eos/user/i/idioniso/mfout/trkeff/trkeff-ang/trkeff-ang.csv",
    inputDir: str = "/eos/user/i/idioniso/mfout/trkeff/trkeff-ang",
    angles: list = [20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
):
    TTs = (1, 11, 3, 13)
    header, data = [], []

    header.append(f"ang")
    for tt in TTs: header.append(f"trkeff_{tt}_tc")
    for tt in TTs: header.append(f"trkeffErr_{tt}_tc")

    for i_ang, ang in enumerate(angles):
        data.append([ang])

        with uproot.open(f"/eos/user/i/idioniso/mfout/trkeff/trkeff-ang/_trkeff_angLt{ang}_Run7080.root") as fin:
            eff, effErr = {}, {}

            for tt in TTs:
                eff[tt]    = np.mean(getNumpyFromUprootTH2(fin[f"eff_x.y_{tt}_7080.tc"])[0]
                effErr[tt] = np.std(getNumpyFromUprootTH2(fin[f"eff_x.y_{tt}_7080.tc"])[0]

            for tt in TTs: data[-1].append(eff[tt])
            for tt in TTs: data[-1].append(effErr[tt])


    saveToCsv(header, data, fout)
