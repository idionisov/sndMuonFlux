from typing import Union
from time import time
from ROOT import TChain, sndRecoTrack, MuFilter, MuFilterHit, TClonesArray, Scifi, SNDLHCEventHeader, TF1
import numpy as np
from ddf.pyfuncs import print_status
from ddf.snd.snd import get_N
from ddf.snd.trk import xy_eff_range, get_intersection, is_good, is_within_ds3, is_within_us5_bar, is_within_veto_bar
from ddf.snd.trkeff import ref1_and_ref2_are_within_allowed_distance, ref1_is_within_eff_area
from selection import is_selected


def loop(
    hists: dict, tree,
    scifi: Scifi,
    mufi:  MuFilter,
    run:   int = 7080,
    track_types = (1, 11, 3, 13),
    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.},
    xz_ang_min: dict = {1: -0.02, 11: -0.02, 3: -0.02, 13: -0.02},
    xz_ang_max: dict = {1:  0.02, 11:  0.02, 3:  0.02, 13:  0.02},
    yz_ang_min: dict = {1: -0.02, 11: -0.02, 3: -0.02, 13: -0.02},
    yz_ang_max: dict = {1:  0.02, 11:  0.02, 3:  0.02, 13:  0.02}
):
    start_time = time()
    count = 0

    tree.GetEntry(0)
    scifi.InitEvent(tree.EventHeader)
    mufi.InitEvent(tree.EventHeader)

    nentries = tree.GetEntries()
    for i_event, event in enumerate(tree):
        count = print_status(i_event, nentries, start_time, count)
        if not event.EventHeader.isIP1(): continue
 
        for tag_trk in event.Reco_MuonTracks:
            if(
                tag_track.getTrackType() != 1 or
                not is_good(
                    tag_track,
                    xz_ang_min  = xz_ang_min,
                    xz_ang_max  = xz_ang_max,
                    yz_ang_min  = yz_ang_min,
                    yz_ang_max  = yz_ang_max
                ) or
                not is_within_us5_bar(tag_track, mf_hits, mf) or
                not is_within_ds3(tag_track)
            ): continue

            ref_tag = get_intersection(tag_trk, z_ref)
            x_tag = ref_tag.X()
            y_tag = ref_tag.Y()

            if not ref1_is_within_eff_area(ref_tag,
                x_min=xy_eff_range["min"]["x"],
                x_max=xy_eff_range["max"]["x"],
                y_min=xy_eff_range["min"]["y"],
                y_max=xy_eff_range["max"]["y"]
            ): return False

            ds_st_list: list = []
            for trk2 in event.Reco_MuonTracks:
                if (trk2.getTrackType()==tt and is_good(trk2)): track2_list.append(trk2)

            if track2_list: continue

            print(i_event)
