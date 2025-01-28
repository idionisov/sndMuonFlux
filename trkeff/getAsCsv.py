import uproot
import argparse
import os, csv
import numpy as np


def saveToCsv(header, data, fout: str):
    foutExists = os.path.isfile(fout)

    with open(fout, mode='a', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)

        if not foutExists or os.stat(fout).st_size == 0:
            writer.writerow(header)

        writer.writerows(data)
    print(f"Output: {fout}")



parser = argparse.ArgumentParser()

parser.add_argument('-o', '--fout', type=str, default="trkeff.csv")

args = parser.parse_args()
fout = f"/eos/user/i/idioniso/mfout/{args.fout}"

if not fout.endswith(".csv"):
    fout = f"{fout}.csv"


with uproot.open("/eos/user/i/idioniso/mfout/trkeff.root") as fin:
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
