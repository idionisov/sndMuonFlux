from ROOT import sndRecoTrack, TClonesArray, MuFilter
from ddf.snd.trk import is_good
from ddf.snd.trkeff import is_within_us5_bar, is_within_ds3, is_within_veto_bar


def is_selected_1(
    tag_track:   sndRecoTrack,
    mf_hits:     TClonesArray,
    mf:          MuFilter,
    xz_ang_min:  float = -0.08,
    xz_ang_max:  float =  0.08,
    yz_ang_min:  float = -0.08,
    yz_ang_max:  float =  0.08
) -> bool:
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
    ): return False
    else: return True



def is_selected_11(
    tag_track:   sndRecoTrack,
    mf_hits:     TClonesArray,
    mf:          MuFilter,
    xz_ang_min:  float = -0.08,
    xz_ang_max:  float =  0.08,
    yz_ang_min:  float = -0.08,
    yz_ang_max:  float =  0.08
) -> bool:
    if (
        tag_track.getTrackType() != 11 or
        not is_good(
            tag_track,
            xz_ang_min  = xz_ang_min,
            xz_ang_max  = xz_ang_max,
            yz_ang_min  = yz_ang_min,
            yz_ang_max  = yz_ang_max
        ) or
        not is_within_us5_bar(tag_track, mf_hits, mf) or
        not is_within_ds3(tag_track)
    ): return False
    else: return True



def is_selected_3(
    tag_track:   sndRecoTrack,
    mf_hits:     TClonesArray,
    mf:          MuFilter,
    xz_ang_min:  float = -0.08,
    xz_ang_max:  float =  0.08,
    yz_ang_min:  float = -0.08,
    yz_ang_max:  float =  0.08
) -> bool:
    if (
        tag_track.getTrackType() != 3 or
        not is_good(
            tag_track,
            xz_ang_min  = xz_ang_min,
            xz_ang_max  = xz_ang_max,
            yz_ang_min  = yz_ang_min,
            yz_ang_max  = yz_ang_max
        ) or
        not is_within_veto_bar(tag_track, mf_hits, mf)
    ): return False
    else: return True



def is_selected_13(
    tag_track:   sndRecoTrack,
    mf_hits:     TClonesArray,
    mf:          MuFilter,
    xz_ang_min:  float = -0.08,
    xz_ang_max:  float =  0.08,
    yz_ang_min:  float = -0.08,
    yz_ang_max:  float =  0.08
) -> bool:
    if (
        tag_track.getTrackType() != 13 or
        not is_good(
            tag_track,
            xz_ang_min  = xz_ang_min,
            xz_ang_max  = xz_ang_max,
            yz_ang_min  = yz_ang_min,
            yz_ang_max  = yz_ang_max
        ) or
        not is_within_veto_bar(tag_track, mf_hits, mf)
    ): return False
    else: return True


def is_selected(
    track_type:  int,
    tag_track:   sndRecoTrack,
    mf_hits:     TClonesArray,
    mufi:        MuFilter,
    xz_ang_min:  dict = {1: -0.08,  11: -0.08,   3: -0.08,   13: -0.08},
    xz_ang_max:  dict = {1:  0.08,  11:  0.08,   3:  0.08,   13:  0.08},
    yz_ang_min:  dict = {1: -0.08,  11: -0.08,   3: -0.08,   13: -0.08},
    yz_ang_max:  dict = {1:  0.08,  11:  0.08,   3:  0.08,   13:  0.08}
) -> bool:

    if track_type==11:
        return is_selected_13(tag_track, mf_hits, mufi,
            xz_ang_min  = xz_ang_min[11],
            xz_ang_max  = xz_ang_max[11],
            yz_ang_min  = yz_ang_min[11],
            yz_ang_max  = yz_ang_max[11]
        )

    elif track_type==1:
        return is_selected_3(tag_track, mf_hits, mufi,
            xz_ang_min  = xz_ang_min[1],
            xz_ang_max  = xz_ang_max[1],
            yz_ang_min  = yz_ang_min[1],
            yz_ang_max  = yz_ang_max[1]
        )

    elif track_type==13:
        return is_selected_11(tag_track, mf_hits, mufi,
            xz_ang_min  = xz_ang_min[13],
            xz_ang_max  = xz_ang_max[13],
            yz_ang_min  = yz_ang_min[13],
            yz_ang_max  = yz_ang_max[13]
        )

    elif track_type==3:
        return is_selected_1(tag_track, mf_hits, mufi,
            xz_ang_min  = xz_ang_min[3],
            xz_ang_max  = xz_ang_max[3],
            yz_ang_min  = yz_ang_min[3],
            yz_ang_max  = yz_ang_max[3]
        )

    else: return False
