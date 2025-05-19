import argparse
from time import time
from SndlhcGeo import GeoInterface
from ROOT import TH1F, TFile, TChain
from sndUtils import *
from ddfRoot import *
from ddfUtils import *
from hists import get_hists, create_h

parser = argparse.ArgumentParser()

parser.add_argument('--run', type=int, default=10241)
parser.add_argument('--tteff', type=int, default=11)
parser.add_argument('--geofile', type=str, default="/eos/experiment/sndlhc/convertedData/physics/2024/run_2412/geofile_sndlhc_TI18_V12_2024.root")

args = parser.parse_args()
geofile = args.geofile
run = args.run
tt = args.tteff
mfout = "/eos/user/i/idioniso/mfout"
z_ref = {1: 430., 11: 450., 3: 430., 13: 450.}

if tt == 1:
    att = 3
elif tt == 11:
    att = 13
elif tt == 3:
    att = 1
elif tt==13:
    att = 11
else:
    exit()

cbmsim = TChain("cbmsim")
input_files = f"/eos/user/i/idioniso/1_Data/Tracks/run_{run:06d}/*f10*.root"
print(f"Input files:\t{input_files}")
cbmsim.Add(input_files)
print(f"Total events:\t{cbmsim.GetEntries()}")

geo   = GeoInterface(geofile)
scifi = geo.modules['Scifi']
mufi  = geo.modules['MuFilter']

start_time = time()

cbmsim.GetEntry(0)
scifi.InitEvent(cbmsim.EventHeader)
mufi.InitEvent(cbmsim.EventHeader)

bad_events = {tt: []}

nTotalEvents = 0
nCandNotBuilt = 0
nCandNotMatched = 0
nBadEvents = 0

nentries = cbmsim.GetEntries()
already_flagged = set()

for i_event, event in enumerate(cbmsim):
    if nBadEvents >= 150:
        break

    if not event.EventHeader.isIP1():
        continue

    eventNum = event.EventHeader.GetEventNumber()
    if eventNum in already_flagged:
        continue

    mf_hits = event.Digi_MuFilterHits

    for tag_trk in event.Reco_MuonTracks:
        tag_trk = DdfTrack(Track=tag_trk, Event=event)

        if not (tag_trk.tt != att and tag_trk.IsGood(xz_min=-0.02, xz_max=0.02, yz_min=-0.02, yz_max=0.02)):
            continue

        ref_tag = tag_trk.GetPointAtZ(z_ref[tt])
        x_tag = ref_tag.X()
        y_tag = ref_tag.Y()

        if not (
            x_tag < -10. and x_tag > -42. and
            y_tag <  48. and y_tag > 19.
        ):
            continue

        nTotalEvents += 1
        trk2_is_built = False
        trk2_is_matched = False
        for trk2 in event.Reco_MuonTracks:
            trk2 = DdfTrack(Track=trk2, Event=event)

            if not (trk2.tt == tt and trk2.Mom.Z() > 0):
                continue
            trk2_is_built = True

            ref2 = trk2.GetPointAtZ(z_ref[tt])
            x2 = ref2.X()
            y2 = ref2.Y()

            if (
                abs(x2 - x_tag) <= 3. and
                abs(y2 - y_tag) <= 3.
            ):
                trk2_is_matched = True

        if trk2_is_built:
            if eventNum not in already_flagged:
                if not trk2_is_matched:
                    print(nBadEvents, eventNum, "\tTracks not paired!")
                    bad_events[tt].append(i_event)
                    nBadEvents += 1
                    nCandNotMatched += 1
                    already_flagged.add(eventNum)  # Avoid duplicate logs

        else:
            if eventNum not in already_flagged:
                print(nBadEvents, eventNum, "\tCandidate track was not built!")
                bad_events[tt].append(i_event)
                nBadEvents += 1
                nCandNotBuilt += 1
                already_flagged.add(eventNum)  # Avoid duplicate logs
print(tt)
print(f"> Inefficient events:  \t{nBadEvents} ({nBadEvents*100/nTotalEvents:.02f}% of all tagged events)")
print(f"> Not built:           \t{nCandNotBuilt} ({nCandNotBuilt*100/nBadEvents:.02f}% of all inefficient events)")
print(f"> Not paired:          \t{nCandNotMatched} ({nCandNotMatched*100/nBadEvents:.02f}% of all inefficient events)")
