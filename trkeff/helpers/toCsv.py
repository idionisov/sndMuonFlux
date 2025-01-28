import uproot
import os, csv

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
        1:  [5,   10,  15,  20, 40,  60, 100],
        11: [8,   12,  16,  25, 50,  80, 120],
        3:  [2,   3,   4,   5,  7.5, 10, 15 ],
        13: [3.5, 4.5, 5.5, 7,  9,   12, 15 ]
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
