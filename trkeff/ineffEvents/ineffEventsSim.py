import argparse

from sndUtils import *
from ddfUtils import *
from ddfRoot import *
from time import time

import ineffUtils
import ineffHists

parser = argparse.ArgumentParser()

parser.add_argument('-z', '--z_ref', nargs="+", type=float, default=[430., 450., 430., 450.])
parser.add_argument('-g', '--geofile', type=str, default="/eos/user/i/idioniso/1_Data/Monte_Carlo/geofile_MC.PbPb.root")

args = parser.parse_args()

zRef = {1: args.z_ref[0], 11: args.z_ref[1], 3: args.z_ref[2], 13: args.z_ref[3]}
trackTypes=(1, 11, 3, 13)
geo = args.geofile


h = {
    "x.y": ineffHists.getIneffXYHists("sim"),
    "n":   ineffHists.getIneffNHists("sim")
}


data = {
    "EMD": SndMCData(
        InputDir = "/eos/user/i/idioniso/1_Data/Monte_Carlo",
        Files = "muonReco_MC-EMD_PbPb.root",
        Geofile = geo
    ),
    "NI": SndMCData(
        InputDir = "/eos/user/i/idioniso/1_Data/Monte_Carlo",
        Files = "muonReco_MC-NI_PbPb.root",
        Geofile = geo
    ),
}

mcFactor = {
    'EMD': 0.1388888888888889,
    'NI':  2.003205128205128
}

ineffEvents = {
    "EMD": {
        "eventNum": [],
        "tt": [],
        "built": [],
        "paired": []
    },
    "NI": {
        "eventNum": [],
        "tt": [],
        "built": [],
        "paired": []
    }
}
t0 = time()
count = 0
trackTypes = (1, 11, 3, 13)
nEntries = data["EMD"].Tree.GetEntries() + data["NI"].Tree.GetEntries()

trkPairs = {}
for mcSet in data:
    data[mcSet].InitGeo()
    for i_event, event in enumerate(data[mcSet].Tree):
        eventNumber = event.EventHeader.GetEventNumber()
        count = printStatus(i_event, nEntries, t0, count)

        if not (event.EventHeader.isIP1() and event.Reco_MuonTracks):
            continue

        w = mcFactor[mcSet] * event.MCTrack[0].GetWeight()
        for tt in (1, 11, 3, 13):
            for trk0 in event.Reco_MuonTracks:
                if not (
                    trk0.getTrackType() == att(tt) and
                    trk0.getTrackMom().Z() > 0 and
                    trk0.getAngleXZ() <= 0.02 and
                    trk0.getAngleYZ() <= 0.02
                ):
                    continue
                else:
                    trk0 = DdfTrack(Track=trk0, Event=event, IP1_Angle=0.02)

                    if trk0.tt==1 or trk0.tt==11:
                        if not (
                            trk0.IsWithinUS5Bar(data[mcSet].Mufi, event.Digi_MuFilterHits) and
                            trk0.IsWithinDS3()
                        ): continue

                    elif (trk0.tt==3 or trk0.tt==13):
                        if not (
                            trk0.IsWithinVetoBar(data[mcSet].Mufi, event.Digi_MuFilterHits)
                        ): continue
                    else: continue


                zRef0 = trk0.GetPointAtZ(zRef[tt])
                x0 = zRef0.X()
                y0 = zRef0.Y()

                if not (
                    x0 > -42. and x0 < -10. and y0 > 19. and y0 < 48.
                ): continue

                N = getN(tt, event)
                h["x.y"][tt]["all"].Fill(x0, y0, w)
                h["n"][tt]["all"].Fill(N, w)

                ineffEvents[mcSet]["eventNum"].append(eventNumber)
                ineffEvents[mcSet]["tt"].append(att(trk0.tt))

                if not ineffUtils.oppositeTrackExists(trk0, event):
                    ineffEvents[mcSet]["built"].append(False)
                    ineffEvents[mcSet]["paired"].append(False)

                    h["x.y"][tt]["notBuilt"].Fill(x0, y0, w)
                    h["x.y"][tt]["ineff"].Fill(x0, y0, w)
                    h["n"][tt]["notBuilt"].Fill(N, w)
                    h["n"][tt]["ineff"].Fill(N, w)
                    continue
                else:
                    ineffEvents[mcSet]["built"].append(True)


                pair = False
                for trk1 in event.Reco_MuonTracks:
                    if not (
                        trk1.getTrackType() == tt and
                        trk1.getTrackMom().Z() > 0
                    ):
                        continue
                    else:
                        trk1 = DdfTrack(
                            Track=trk1,
                            Event=event,
                            IP1_Angle=0.02,
                        )


                        if ineffUtils.tracksArePaired(trk0, trk1, zRef, tt):
                            pair = True
                            break

                ineffEvents[mcSet]["paired"].append(pair)

                if not pair:
                    h["x.y"][tt]["notPaired"].Fill(x0, y0, w)
                    h["x.y"][tt]["ineff"].Fill(x0, y0, w)
                    h["n"][tt]["notPaired"].Fill(N, w)
                    h["n"][tt]["ineff"].Fill(N, w)
                else:
                    h["x.y"][tt]["eff"].Fill(x0, y0, w)
                    h["n"][tt]["eff"].Fill(N, w)

ineffEvents = {
    k: pd.DataFrame(ineffEvents[k]) for k in data
}



h_div_nb = {"x.y": {}, "n": {}}
h_div_np = {"x.y": {}, "n": {}}
for tt in trackTypes:

    h_div_nb["n"][tt] = h["n"][tt]["notBuilt"].Clone()
    h_div_nb["n"][tt].Divide(h["n"][tt]["all"])
    h_div_nb["n"][tt].SetName(f"{h_div_nb['n'][tt].GetName()}.DivBy.All")

    h_div_np["n"][tt] = h["n"][tt]["notPaired"].Clone()
    h_div_np["n"][tt].Divide(h["n"][tt]["all"])
    h_div_np["n"][tt].SetName(f"{h_div_np['n'][tt].GetName()}.DivBy.All")

    h_div_nb["x.y"][tt] = h["x.y"][tt]["notBuilt"].Clone()
    h_div_nb["x.y"][tt].Divide(h["x.y"][tt]["all"])
    h_div_nb["x.y"][tt].SetName(f"{h_div_nb['n'][tt].GetName()}.DivBy.All")

    h_div_np["x.y"][tt] = h["x.y"][tt]["notPaired"].Clone()
    h_div_np["x.y"][tt].Divide(h["x.y"][tt]["all"])
    h_div_np["x.y"][tt].SetName(f"{h_div_np['n'][tt].GetName()}.DivBy.All")


saveToRoot(h, h_div_nb, h_div_np, fout="/home/idioniso/ineffEvents/_ineff-sim.root")
