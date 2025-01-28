import ROOT
import uproot
import os, csv
import numpy as np
import pandas as pd
from time import time
from sndUtils import DdfTrack, SndMCData, sys, alg



def getEffWithErr(csv: str, tt: int) -> tuple:
    df = pd.read_csv(csv, index_col="run")
    eff    = df.at[-1, f"trkeff_{tt}"]
    effErr = df.at[-1, f"trkeffErr_{tt}"]

    return eff, effErr


def saveToCsv(header, data, fout: str):
    foutExists = os.path.isfile(fout)

    with open(fout, mode='a', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)

        if not foutExists or os.stat(fout).st_size == 0:
            writer.writerow(header)

        writer.writerows(data)
    print(f"Output: {fout}")



def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-dir', type=str, default="/eos/user/i/idioniso/1_Data/Monte_Carlo")
    parser.add_argument('--trkeff', type=str, default="/eos/user/i/idioniso/mfout/trkeff.csv")

    args = parser.parse_args()

    mfout = "/eos/user/i/idioniso/mfout"
    track_types = (1, 11, 3, 13)
    xy = {
        'sf': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}},
        'ds': {'min': {'x': -42., 'y': 19.}, 'max': {'x': -10., 'y': 48.}}
    }
    z_ref = {1: 430., 11: 450., 3: 430., 13: 450.}
    area = {
        sys: abs(xy[sys]['max']['x']-xy[sys]['min']['x']) * abs(xy[sys]['max']['y']-xy[sys]['min']['y']) for sys in ("sf", "ds")
    }



    data = {
        a: SndMCData(InputDir=args.input_dir, Files=f"*{a}*PbPb*.root") for a in ("EMD", "NI")
    }
    data["EMD"].Print()
    data["NI"].Print()


    L_LHC = 6.4e-6
    sigma = {"NI": 7.8e9, "EMD": 4.5e11}
    N_rate = {"NI": 1e5, "EMD": 4e5}
    L_MC = {i: N_rate[i]/sigma[i] for i in ("EMD", "NI")}

    factor = {i: L_MC[i]/L_LHC for i in ("EMD", "NI")}



    eps = {}
    for tt in track_types:
        eps[tt] = {}
        eps[tt]['v'], eps[tt]['e'] = getEffWithErr(args.trkeff, tt)


    eventMCmu = {i: {} for i in ("EMD", "NI")}
    Nrate = {
        i: {tt: 0 for tt in track_types} for i in ("EMD", "NI")
    }


    totalEvents = data["EMD"].Tree.GetEntries() + data["NI"].Tree.GetEntries()
    for i in ("EMD", "NI"):
        for i_event, event in enumerate(data[i].Tree):
            if i=="EMD": current_event = i_event
            else:        current_event = i_event + data["EMD"].Tree.GetEntries()

            if current_event%10000==0:
                print(f"{current_event*100/totalEvents:.02f}%")

            for mctrack in event.MCTrack:
                if mctrack.GetMotherId()==-1:
                    eventMCmu[i][i_event] = mctrack.GetWeight()/factor[i]

            for trk in event.Reco_MuonTracks:
                trk = DdfTrack(Track=trk, Event=event)

                tt = trk.tt
                if not (
                    trk.IsIP1() and
                    trk.IsWithinAref(
                        Zref=z_ref[tt],
                        xmin=xy[sys(tt)]["min"]["x"],
                        xmax=xy[sys(tt)]["max"]["x"],
                        ymin=xy[sys(tt)]["min"]["y"],
                        ymax=xy[sys(tt)]["max"]["y"]
                    )
                ): continue

                ref = trk.GetPointAtZ(Z=z_ref[trk.tt])

                Nrate[i][tt] += eventMCmu[i][i_event]


    mf = {}
    for i in ("EMD", "NI"):
        mf[i] = {}

        for tt in track_types:
            A = area[sys(tt)]
            Nr = Nrate[i][tt]

            mf[i][tt] = (1e6*Nr)/(6.4*A*eps[tt]['v'])


    dNrate = {
        i: {
            tt: np.sqrt(Nrate[i][tt]) for tt in track_types
        } for i in ("EMD", "NI")
    }

    dmf = {}
    for i in ("EMD", "NI"):
        dmf[i] = {}

        for tt in track_types:
            A = area[sys(tt)]
            Nr = Nrate[i][tt]
            dNr = dNrate[i][tt]

            dmf[i][tt] = (1e6/(A*eps[tt]['v']*6.4)) * np.sqrt(
                (Nr**2 * eps[tt]['e']**2)/eps[tt]['v']**2 + Nr
            )


    header = ["run"]
    data = [[-1]]
    for tt in track_types:
        header.append(f"flux_{tt}")
        data[0].append(mf["EMD"][tt] + mf["NI"][tt])

    for tt in (1, 11, 3, 13):
        flux = mf["EMD"][tt] + mf["NI"][tt]
        fluxErr = dmf["EMD"][tt] + dmf["NI"][tt]

        header.append(f"fluxErr_{tt}")
        data[0].append(fluxErr)

        print(f" > {tt}:\t{flux/1e4} ± {fluxErr/1e4}    \033[1;39m[10⁴ nb/cm²]\033[0m")

    saveToCsv(header, data, f"{mfout}/muonFlux.csv")


if __name__=="__main__":
    main()
