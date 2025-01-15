from typing import Union
from time import time
from ROOT import TChain, sndRecoTrack, MuFilter, MuFilterHit, TClonesArray, Scifi, SNDLHCEventHeader, TF1
import numpy as np
from ddf.pyfuncs import print_status
from ddf.snd.snd import get_N
from ddf.snd.trk import xy_eff_range, get_intersection, is_good, is_within_ds3, is_within_us5_bar, is_within_veto_bar
from ddf.snd.trkeff import ref1_and_ref2_are_within_allowed_distance, ref1_is_within_eff_area
from selection import is_selected

def ref1_is_within_eff_area(
    ref1: TVector3,
    x_min: float = -43.,
    x_max: float = -10.,
    y_min: float = 18.,
    y_max: float = 50.
) -> bool:
    if (
        ref1.X() >= x_min and
        ref1.X() <= x_max and
        ref1.Y() >= y_min and
        ref1.Y() <= y_max
    ): return True
    else: return False


def fill_hists(
    h:            dict,
    event:        TChain,
    tag_trk:      sndRecoTrack,
    run_or_mcSet: Union[int, str] = 7080,
    tt:           int = 1,
    z_ref:        float = 490.,
    weight:       Union[int, float] = 1
):
    ref_tag = get_intersection(tag_trk, z_ref)
    x_tag = ref_tag.X()
    y_tag = ref_tag.Y()

    d0 = np.sqrt(x_tag*x_tag + y_tag*y_tag)

    xz_tag = tag_trk.getAngleXZ()*1e3
    yz_tag = tag_trk.getAngleYZ()*1e3

    chi2_tag = tag_trk.getChi2()
    chi2ndf_tag = tag_trk.getChi2Ndf()

    trkP_tag = tag_trk.getTrackPoints().size()
    n_tag = get_N(tt, event)

    prpts = h[tt].keys()

    if "x.y" in prpts: h[tt]["x.y"][1].Fill(x_tag, y_tag, weight)
    if ref1_is_within_eff_area(ref_tag,
        x_min=xy_eff_range["min"]["x"],
        x_max=xy_eff_range["max"]["x"],
        y_min=xy_eff_range["min"]["y"],
        y_max=xy_eff_range["max"]["y"]
    ):
        if "x"       in prpts: h[tt]["x"][1].Fill(x_tag, weight)
        if "y"       in prpts: h[tt]["y"][1].Fill(y_tag, weight)
        if "xz"      in prpts: h[tt]['xz'][1].Fill(xz_tag, weight)
        if "yz"      in prpts: h[tt]['yz'][1].Fill(yz_tag, weight)
        if "chi2"    in prpts: h[tt]['chi2'][1].Fill(chi2_tag, weight)
        if "chi2ndf" in prpts: h[tt]['chi2ndf'][1].Fill(chi2ndf_tag, weight)
        if "n"       in prpts: h[tt]['n'][1].Fill(n_tag, weight)
        if "trkP"    in prpts: h[tt]['trkP'][1].Fill(trkP_tag, weight)
        if "d0"      in prpts: h[tt]['d0'][1].Fill(d0, weight)

        if "x.xz"       in prpts: h[tt]['x.xz'][1].Fill(x_tag, xz_tag, weight)
        if "x.yz"       in prpts: h[tt]['x.yz'][1].Fill(x_tag, yz_tag, weight)
        if "x.chi2"     in prpts: h[tt]['x.chi2'][1].Fill(x_tag, chi2_tag, weight)
        if "x.chi2ndf"  in prpts: h[tt]['x.chi2ndf'][1].Fill(x_tag, chi2ndf_tag, weight)
        if "x.trkP"     in prpts: h[tt]['x.trkP'][1].Fill(x_tag, trkP_tag, weight)
        if "x.n"        in prpts: h[tt]['x.n'][1].Fill(x_tag, n_tag, weight)

        if "y.xz"       in prpts: h[tt]['y.xz'][1].Fill(y_tag, xz_tag, weight)
        if "y.yz"       in prpts: h[tt]['y.yz'][1].Fill(y_tag, yz_tag, weight)
        if "y.chi2"     in prpts: h[tt]['y.chi2'][1].Fill(y_tag, chi2_tag, weight)
        if "y.chi2ndf"  in prpts: h[tt]['y.chi2ndf'][1].Fill(y_tag, chi2ndf_tag, weight)
        if "y.trkP"     in prpts: h[tt]['y.trkP'][1].Fill(y_tag, trkP_tag, weight)
        if "y.n"        in prpts: h[tt]['y.n'][1].Fill(y_tag, n_tag, weight)

        if "xz.x"       in prpts: h[tt]['xz.x'][1].Fill(xz_tag, x_tag, weight)
        if "xz.y"       in prpts: h[tt]['xz.y'][1].Fill(xz_tag, y_tag, weight)
        if "xz.yz"      in prpts: h[tt]["xz.yz"][1].Fill(xz_tag, yz_tag, weight)
        if "xz.chi2"    in prpts: h[tt]['xz.chi2'][1].Fill(xz_tag, chi2_tag, weight)
        if "xz.chi2ndf" in prpts: h[tt]['xz.chi2ndf'][1].Fill(xz_tag, chi2ndf_tag, weight)
        if "xz.trkP"    in prpts: h[tt]['xz.trkP'][1].Fill(xz_tag, trkP_tag, weight)
        if "xz.n"       in prpts: h[tt]['xz.n'][1].Fill(xz_tag, n_tag, weight)

        if "yz.x"       in prpts: h[tt]['yz.x'][1].Fill(yz_tag, x_tag, weight)
        if "yz.y"       in prpts: h[tt]['yz.y'][1].Fill(yz_tag, y_tag, weight)
        if "yz.chi2"    in prpts: h[tt]['yz.chi2'][1].Fill(yz_tag, chi2_tag, weight)
        if "yz.chi2ndf" in prpts: h[tt]['yz.chi2ndf'][1].Fill(yz_tag, chi2ndf_tag, weight)
        if "yz.trkP"    in prpts: h[tt]['yz.trkP'][1].Fill(yz_tag, trkP_tag, weight)
        if "yz.n"       in prpts: h[tt]['yz.n'][1].Fill(yz_tag, n_tag, weight)


        track2_list: list = []
        for trk2 in event.Reco_MuonTracks:
            if (
                trk2.getTrackType()==tt and
                is_good(trk2, xz_ang_min=-0.02, xz_ang_max=0.02, yz_ang_min=-0.02, yz_ang_max=0.02)
            ): track2_list.append(trk2)

        if track2_list:
            for trk2 in track2_list:
                xz_cand = 1e3*trk2.getAngleXZ()
                yz_cand = 1e3*trk2.getAngleYZ()
                ref2 = get_intersection(trk2, z_ref)

                if "dxRef" in prpts: h[tt]['dxRef'].Fill(ref2.X()-x_tag, weight)
                if "dyRef" in prpts: h[tt]['dyRef'].Fill(ref2.Y()-y_tag, weight)
                if "dxz"   in prpts: h[tt]['dxz'].Fill(xz_cand - xz_tag, weight)
                if "dyz"   in prpts: h[tt]['dyz'].Fill(yz_cand - yz_tag, weight)

                if not ref1_and_ref2_are_within_allowed_distance(ref_tag, ref2, 3.): continue

                if "x.y" in prpts: h[tt]["x.y"][0].Fill(x_tag, y_tag, weight)

                if not ref1_is_within_eff_area(ref_tag): continue
                if "x"       in prpts: h[tt]["x"][0].Fill(x_tag, weight)
                if "y"       in prpts: h[tt]["y"][0].Fill(y_tag, weight)
                if "xz"      in prpts: h[tt]['xz'][0].Fill(xz_tag, weight)
                if "yz"      in prpts: h[tt]['yz'][0].Fill(yz_tag, weight)
                if "chi2"    in prpts: h[tt]['chi2'][0].Fill(chi2_tag, weight)
                if "chi2ndf" in prpts: h[tt]['chi2ndf'][0].Fill(chi2ndf_tag, weight)
                if "n"       in prpts: h[tt]['n'][0].Fill(n_tag, weight)
                if "trkP"    in prpts: h[tt]['trkP'][0].Fill(trkP_tag, weight)
                if "d0"      in prpts: h[tt]['d0'][0].Fill(d0, weight)

                if "x.xz"       in prpts: h[tt]['x.xz'][0].Fill(x_tag, xz_tag, weight)
                if "x.yz"       in prpts: h[tt]['x.yz'][0].Fill(x_tag, yz_tag, weight)
                if "x.chi2"     in prpts: h[tt]['x.chi2'][0].Fill(x_tag, chi2_tag, weight)
                if "x.chi2ndf"  in prpts: h[tt]['x.chi2ndf'][0].Fill(x_tag, chi2ndf_tag, weight)
                if "x.trkP"     in prpts: h[tt]['x.trkP'][0].Fill(x_tag, trkP_tag, weight)
                if "x.n"        in prpts: h[tt]['x.n'][0].Fill(x_tag, n_tag, weight)

                if "y.xz"       in prpts: h[tt]['y.xz'][0].Fill(y_tag, xz_tag, weight)
                if "y.yz"       in prpts: h[tt]['y.yz'][0].Fill(y_tag, yz_tag, weight)
                if "y.chi2"     in prpts: h[tt]['y.chi2'][0].Fill(y_tag, chi2_tag, weight)
                if "y.chi2ndf"  in prpts: h[tt]['y.chi2ndf'][0].Fill(y_tag, chi2ndf_tag, weight)
                if "y.trkP"     in prpts: h[tt]['y.trkP'][0].Fill(y_tag, trkP_tag, weight)
                if "y.n"        in prpts: h[tt]['y.n'][0].Fill(y_tag, n_tag, weight)

                if "xz.x"       in prpts: h[tt]['xz.x'][0].Fill(xz_tag, x_tag, weight)
                if "xz.y"       in prpts: h[tt]['xz.y'][0].Fill(xz_tag, y_tag, weight)
                if "xz.yz"      in prpts: h[tt]["xz.yz"][0].Fill(xz_tag, yz_tag, weight)
                if "xz.chi2"    in prpts: h[tt]['xz.chi2'][0].Fill(xz_tag, chi2_tag, weight)
                if "xz.chi2ndf" in prpts: h[tt]['xz.chi2ndf'][0].Fill(xz_tag, chi2ndf_tag, weight)
                if "xz.trkP"    in prpts: h[tt]['xz.trkP'][0].Fill(xz_tag, trkP_tag, weight)
                if "xz.n"       in prpts: h[tt]['xz.n'][0].Fill(xz_tag, n_tag, weight)

                if "yz.x"       in prpts: h[tt]['yz.x'][0].Fill(yz_tag, x_tag, weight)
                if "yz.y"       in prpts: h[tt]['yz.y'][0].Fill(yz_tag, y_tag, weight)
                if "yz.chi2"    in prpts: h[tt]['yz.chi2'][0].Fill(yz_tag, chi2_tag, weight)
                if "yz.chi2ndf" in prpts: h[tt]['yz.chi2ndf'][0].Fill(yz_tag, chi2ndf_tag, weight)
                if "yz.trkP"    in prpts: h[tt]['yz.trkP'][0].Fill(yz_tag, trkP_tag, weight)
                if "yz.n"       in prpts: h[tt]['yz.n'][0].Fill(yz_tag, n_tag, weight)




def loop(
    hists: dict, tree,
    scifi: Scifi,
    mufi:  MuFilter,
    run:   int = 7080,
    track_types = (1, 11, 3, 13),
    z_ref: dict = {1: 440., 11: 435., 3: 460., 13: 475.},
    xz_min: dict = {1: -0.02, 11: -0.02, 3: -0.02, 13: -0.02},
    xz_max: dict = {1:  0.02, 11:  0.02, 3:  0.02, 13:  0.02},
    yz_min: dict = {1: -0.02, 11: -0.02, 3: -0.02, 13: -0.02},
    yz_max: dict = {1:  0.02, 11:  0.02, 3:  0.02, 13:  0.02}
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

        for tt in track_types:
            for tag_trk in event.Reco_MuonTracks:

                if not is_selected(
                    track_type = tt,
                    tag_track = tag_trk,
                    mf_hits = event.Digi_MuFilterHits,
                    mufi = mufi,
                    xz_ang_min = xz_min,
                    xz_ang_max = xz_max,
                    yz_ang_min = yz_min,
                    yz_ang_max = yz_max
                ): continue

                fill_hists(
                    h = hists,
                    event = event,
                    tag_trk = tag_trk,
                    run_or_mcSet = run,
                    tt = tt,
                    z_ref = z_ref[tt]
                )


def get_fit_eq(
    teff: dict, run: int=7080, track_types=(1, 11, 3, 13)
):
    eq = {
        tt: {
            xy: TF1(f"fit_{xy}_{tt}_{run}", "[0]", xy_eff_range["min"][xy], xy_eff_range["max"][xy])
            for xy in ["x", "y"]
        } for tt in track_types
    }

    for tt in track_types:
        for xy in ["x", "y"]:
            eq[tt][xy].SetParameters(0.9)
            teff[tt][xy].Fit(eq[tt][xy], "SRQ0+")
    return eq
