import argparse

from sndUtils import *
from ddfUtils import *
from ddfRoot import *
from time import time

import ineffUtils
import ineffHists

parser = argparse.ArgumentParser()

parser.add_argument('-r', '--runs', nargs="+", type=int, default=[10337])
parser.add_argument('-f', '--files', type=str, default='*f10*.root')
parser.add_argument('-z', '--z_ref', nargs="+", type=float, default=[430., 450., 430., 450.])
parser.add_argument('-g', '--geofile', type=str, default="/eos/experiment/sndlhc/convertedData/physics/2024/run_2412/geofile_sndlhc_TI18_V12_2024.root")

args = parser.parse_args()

zRef = {1: args.z_ref[0], 11: args.z_ref[1], 3: args.z_ref[2], 13: args.z_ref[3]}
runs = [run for run in args.runs if getSubDirPath(RootDir="/eos/user/i/idioniso/1_Data/Tracks/", TopDir=f"run_{run:06d}")]

trackTypes=(1, 11, 3, 13)
files = args.files
geo = args.geofile


if len(runs)==1:
    h = {
        "x.y":  ineffHists.getIneffXYHists(runs[0]),
        "n":    ineffHists.getIneffNHists(runs[0]),
        "maxH": ineffHists.getIneffMaxPlaneHitsHists(runs[0])
    }
else:
    h = {
        "x.y":  ineffHists.getIneffXYHists(2024),
        "n":    ineffHists.getIneffNHists(2024),
        "maxH": ineffHists.getIneffMaxPlaneHitsHists(2024)
    }

if len(runs)==1:
    if runs[0]==7080:
        geo="/eos/experiment/sndlhc/convertedData/physics/2023/geofile_sndlhc_TI18_V4_2023.root"
    elif runs[0]==10337:
        geo="/eos/experiment/sndlhc/convertedData/physics/2024/run_2412/geofile_sndlhc_TI18_V12_2024.root"
    else:
        raise ValueError("A geofile needs to be provided!")

if not files.endswith(".root"):
    files = f"{files}.root"

data = {}
for run in runs:
    data[run] = SndData(
        Run=run,
        InputDir="/eos/user/i/idioniso/1_Data/Tracks/",
        Files=files,
        Geofile=geo
    )
print(data)



counts = {
    a: {
        tt: 0 for tt in trackTypes
    } for a in ("all", "eff", "ineff", "notBuilt", "notPaired", "cand>20mrad")
}

ineffEvents = {
    "run": [],
    "eventNum": [],
    "tt": [],
    "built": [],
    "paired": [],
    "cand>20mrad": []
}
t0 = time()
count = 0
trackTypes = (1, 11, 3, 13)

trkPairs = {}
for run in runs:
    if run not in data.keys():
        print(f"Skipping run {run}")
        continue

    data[run].InitGeo()
    nEntries = data[run].Tree.GetEntries()

    for i_event, event in enumerate(data[run].Tree):
        eventNumber = event.EventHeader.GetEventNumber()
        count = printStatus(i_event, nEntries, t0, count)

        if not (event.EventHeader.isIP1() and event.Reco_MuonTracks):
            continue

        for tt in (1, 11, 3, 13):
            flag = False
            built = False
            paired = False
            lt20mrad = False

            for trk0 in event.Reco_MuonTracks:
                trk0 = DdfTrack(Track=trk0, Event=event, IP1_Angle=0.02)
                if not (
                    trk0.tt == att(tt) and
                    trk0.Mom.Z() > 0 and
                    abs(trk0.XZ) <= 0.08 and
                    abs(trk0.YZ) <= 0.08
                ):
                    continue
                else:
                    if trk0.tt==1 or trk0.tt==11:
                        if not (
                            trk0.IsWithinUS5Bar(data[run].Mufi, event.Digi_MuFilterHits) and
                            trk0.IsWithinDS3()
                        ): continue

                    elif (trk0.tt==3 or trk0.tt==13):
                        if not (
                            trk0.IsWithinVetoBar(data[run].Mufi, event.Digi_MuFilterHits)
                        ): continue
                    else: continue


                zRef0 = trk0.GetPointAtZ(zRef[tt])
                x0 = zRef0.X()
                y0 = zRef0.Y()

                if not (
                    x0 > -42. and x0 < -10. and y0 > 19. and y0 < 48.
                ):
                    continue

                flag = True

                N = getN(tt, event)
                if tt in (1, 11):
                    nMaxPlane = ineffUtils.getMaxSfHitsPerPlane(event)
                else:
                    nMaxPlane = ineffUtils.getMaxDsHitsPerPlane(event)


                counts["all"][tt] += 1

                if not ineffUtils.oppositeTrackExists(trk0, event):
                    built = False
                    paired = False
                    lt20mrad = False

                    counts["ineff"][tt] += 1
                    counts["notBuilt"][tt] += 1

                else:
                    built = True
                    paired = False
                    lt20mrad = False

                for trk1 in event.Reco_MuonTracks:
                    if not (
                        trk1.getTrackType() == tt and
                        trk1.getTrackMom().Z() > 0 and
                        abs(trk1.getAngleXZ()) <= 0.08 and
                        abs(trk1.getAngleYZ()) <= 0.08
                    ):
                        continue

                    else:
                        trk1 = DdfTrack(
                            Track=trk1,
                            Event=event,
                            IP1_Angle=0.02,
                        )


                        if ineffUtils.tracksArePaired(trk0, trk1, zRef, tt):
                            paired = True

                            if abs(trk1.XZ) <= 0.02 and abs(trk1.YZ) <= 0.02:
                                lt20mrad = True
                            else:
                                lt20mrad = False
                            break

                        else:
                            paired = False
                            lt20mrad = False

                if not paired:
                    counts["ineff"][tt] += 1
                    counts["notPaired"][tt] += 1

                else:
                    if lt20mrad:
                        counts["eff"][tt] += 1
                    else:
                        counts["cand>20mrad"][tt] += 1

            if not flag:
                continue

            ineffEvents["run"].append(run)
            ineffEvents["eventNum"].append(eventNumber)
            ineffEvents["tt"].append(tt)

            h["x.y"][tt]["all"].Fill(x0, y0)
            h["n"][tt]["all"].Fill(N)
            h["maxH"][tt]["all"].Fill(nMaxPlane)

            if not built:
                ineffEvents["built"].append(False)
                ineffEvents["paired"].append(False)
                ineffEvents["cand>20mrad"].append(False)

                h["x.y"][tt]["ineff"].Fill(x0, y0)
                h["n"][tt]["ineff"].Fill(N)
                h["maxH"][tt]["ineff"].Fill(nMaxPlane)

                h["x.y"][tt]["notBuilt"].Fill(x0, y0)
                h["n"][tt]["notBuilt"].Fill(N)
                h["maxH"][tt]["notBuilt"].Fill(nMaxPlane)
            else:
                ineffEvents["built"].append(False)
            
                if paired:
                    ineffEvents["paired"].append(True)

                    h["x.y"][tt]["eff"].Fill(x0, y0)
                    h["n"][tt]["eff"].Fill(N)
                    h["maxH"][tt]["eff"].Fill(nMaxPlane)

                    if lt20mrad:
                        ineffEvents["cand>20mrad"].append(False)
                    else:
                        ineffEvents["cand>20mrad"].append(True)

                        h["x.y"][tt]["cand>20mrad"].Fill(x0, y0)
                        h["n"][tt]["cand>20mrad"].Fill(N)
                        h["maxH"][tt]["cand>20mrad"].Fill(nMaxPlane)
                else:
                    ineffEvents["paired"].append(False)
                    ineffEvents["cand>20mrad"].append(False)

                    h["x.y"][tt]["ineff"].Fill(x0, y0)
                    h["n"][tt]["ineff"].Fill(N)
                    h["maxH"][tt]["ineff"].Fill(nMaxPlane)

                    h["x.y"][tt]["notPaired"].Fill(x0, y0)
                    h["n"][tt]["notPaired"].Fill(N)
                    h["maxH"][tt]["notPaired"].Fill(nMaxPlane)

ineffEvents = pd.DataFrame(ineffEvents)


def saveHists(fout: str):
    h_div_nb = {"x.y": {}, "n": {}, "maxH": {}}
    h_div_np = {"x.y": {}, "n": {}, "maxH": {}}
    for tt in trackTypes:

        h_div_nb["n"][tt] = getTEff(
            passed = h["n"][tt]["notBuilt"],
            total  = h["n"][tt]["all"],
            statOption = "kfcp",
            name = f'{h["n"][tt]["notBuilt"].GetName()}.DivBy.All'
        )
        h_div_nb["n"][tt] = getGraphFromTEff(h_div_nb["n"][tt])

        h_div_np["n"][tt] = getTEff(
            passed = h["n"][tt]["notPaired"],
            total  = h["n"][tt]["all"],
            statOption = "kfcp",
            name = f'{h["n"][tt]["notPaired"].GetName()}.DivBy.All'
        )
        h_div_np["n"][tt] = getGraphFromTEff(h_div_np["n"][tt])

        h_div_nb["x.y"][tt] = getTEff(
            passed = h["x.y"][tt]["notBuilt"],
            total  = h["x.y"][tt]["all"],
            statOption = "kfcp",
            name = f'{h["x.y"][tt]["notBuilt"].GetName()}.DivBy.All'
        )
        h_div_nb["x.y"][tt] = getHistFromTEff2D(h_div_nb["x.y"][tt])

        h_div_np["x.y"][tt] = getTEff(
            passed = h["x.y"][tt]["notPaired"],
            total  = h["x.y"][tt]["all"],
            statOption = "kfcp",
            name = f'{h["x.y"][tt]["notPaired"].GetName()}.DivBy.All'
        )
        h_div_np["x.y"][tt] = getHistFromTEff2D(h_div_np["x.y"][tt])

        h_div_nb["maxH"][tt] = getTEff(
            passed = h["maxH"][tt]["notBuilt"],
            total  = h["maxH"][tt]["all"],
            statOption = "kfcp",
            name = f'{h["maxH"][tt]["notBuilt"].GetName()}.DivBy.All'
        )
        h_div_nb["maxH"][tt] = getGraphFromTEff(h_div_nb["maxH"][tt])

        h_div_np["maxH"][tt] = getTEff(
            passed = h["maxH"][tt]["notPaired"],
            total  = h["maxH"][tt]["all"],
            statOption = "kfcp",
            name = f'{h["maxH"][tt]["notPaired"].GetName()}.DivBy.All'
        )
        h_div_np["maxH"][tt] = getGraphFromTEff(h_div_np["maxH"][tt])

    saveToRoot(h, h_div_nb, h_div_np, fout=fout)

# hh = h["n"].copy()
# for k1 in hh:
#     for k2 in hh[k1]:
#         hh[k1][k2] = getAsPandas(hh[k1][k2])

# def plot_Ineff(tt: int):
#     plt.figure(figsize=(8,6))

#     plt.plot(hh[tt]["ineff"]["x"], hh[tt]["ineff"]["y"], label="ineff")
#     plt.plot(hh[tt]["notPaired"]["x"], hh[tt]["notPaired"]["y"], label="notPaired")
#     plt.plot(hh[tt]["notBuilt"]["x"], hh[tt]["notBuilt"]["y"], label="notBuilt")

#     plt.tight_layout()
#     plt.legend()

#     plt.show()
