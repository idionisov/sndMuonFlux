import argparse
from time import time
from SndlhcGeo import GeoInterface
from ROOT import TH1F, TFile, TChain
from ddf.root.misc import save_to_root
from ddf.root.teff import get_teff_dict
from ddf.snd.trk import is_good, is_within_ds3, is_within_us5_bar, get_intersection, sys, alg, sys_name, alg_name
from ddf.snd.trkeff import ref1_and_ref2_are_within_allowed_distance, ref1_is_within_eff_area, xy_eff_range
from ddf.pyfuncs import print_status, get_current_dmy, get_cl_sigma
from hists import get_hists, create_h
from process import loop, get_fit_eq

parser = argparse.ArgumentParser()

parser.add_argument('--run', type=int, default=7080)
parser.add_argument('--tteff', type=int, default=13)

args = parser.parse_args()
run = args.run
tt = args.tteff
mfout = "/eos/user/i/idioniso/mfout"
cl = get_cl_sigma(1)
z_ref = {3: 460., 13: 475.}

if   tt==3:  att=1
elif tt==13: att=11
else: exit()

cbmsim = TChain("cbmsim")
input_files = f"/eos/user/i/idioniso/1_Data/Tracks/run_{run:06d}/*.root"
print(f"Input files:\t{input_files}")
cbmsim.Add(input_files)
print(f"Total events:\t{cbmsim.GetEntries()}")

geo   = GeoInterface("$SND_DATA/2023/geofile_sndlhc_TI18_V4_2023.root")
scifi = geo.modules['Scifi']
mufi  = geo.modules['MuFilter']

start_time = time()

cbmsim.GetEntry(0)
scifi.InitEvent(cbmsim.EventHeader)
mufi.InitEvent(cbmsim.EventHeader)

bad_events = {tt: []}

badEventsNr = 0
nentries = cbmsim.GetEntries()
for i_event, event in enumerate(cbmsim):
    if badEventsNr >= 500: break


    if not event.EventHeader.isIP1(): continue
    mf_hits = event.Digi_MuFilterHits

    for tag_trk in event.Reco_MuonTracks:
        if not (
            tag_trk.getTrackType() != att and
            is_good(
                tag_trk,
                xz_ang_min  = -0.02,
                xz_ang_max  =  0.02,
                yz_ang_min  = -0.02,
                yz_ang_max  =  0.02
            ) and
            is_within_us5_bar(tag_trk, mf_hits, mufi) and
            is_within_ds3(tag_trk)
        ): continue

        ref_tag = get_intersection(tag_trk, z_ref[tt])
        x_tag = ref_tag.X()
        y_tag = ref_tag.Y()

        if not ref1_is_within_eff_area(ref_tag,
            x_min=xy_eff_range["min"]["x"],
            x_max=xy_eff_range["max"]["x"],
            y_min=xy_eff_range["min"]["y"],
            y_max=xy_eff_range["max"]["y"]
        ): continue

        track2_list: list = []
        for trk2 in event.Reco_MuonTracks:
            if (trk2.getTrackType()==tt and is_good(trk2)): track2_list.append(trk2)

        if track2_list:
            for trk2 in track2_list:
                if not trk2.getTrackMom().Z(): continue
                ref2 = get_intersection(trk2, z_ref[tt])

                if not ref1_and_ref2_are_within_allowed_distance(ref_tag, ref2, 3.):
                    print(badEventsNr, i_event, "\tTracks not paired!")
                    continue

        print(badEventsNr, i_event)
        bad_events[tt].append(i_event)
        badEventsNr += 1
        break
