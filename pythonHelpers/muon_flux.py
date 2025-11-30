import csv
import os
from array import array

import numpy as np
import ROOT


def compute_muon_flux(
    n_tracks: int, lumi: float, area: float,
    eff: float = 1.,
    eff_err: float = 0.,
    lumi_err: float = 0.,
    scale: int = 1,
    bunch_correction_factor: float = 1.,
):
    K = scale * bunch_correction_factor

    muon_flux = K * n_tracks / (lumi * area * eff)

    rel_N = 1.0 / np.sqrt(n_tracks) if n_tracks > 0 else 0.0
    rel_L = lumi_err / lumi
    rel_eff = eff_err / eff

    rel_total = np.sqrt(rel_N**2 + rel_L**2 + rel_eff**2)

    muon_flux_err = muon_flux * rel_total

    return muon_flux, muon_flux_err



def compute_fluxes_per_track_type(
    track_counts: dict,
    lumi: float,
    area: float,
    effs: dict,
    lumi_err: float,
    scale: int,
    correction: dict
) -> dict:
    muon_flux = {}
    for tt in track_counts["IP1"]:
        muon_flux[tt] = compute_muon_flux(
            track_counts["IP1"][tt],
            lumi, area,
            effs[tt][0], effs[tt][1],
            lumi_err, scale, correction[tt]
        )
    return muon_flux


def save_all_to_csv(
    filename: str, run: int, fill: int, scale: int,
    nTracks_dict: dict,
    effs_dict: dict,
    muon_flux_dict: dict,
    cols: list[str] = [
        "Run", "Fill", "scale",
        "nTracks1", "nTracks11", "nTracks3", "nTracks13",
        "trkeff1", "trkeff11", "trkeff3", "trkeff13",
        "trkeffErr1", "trkeffErr11", "trkeffErr3", "trkeffErr13",
        "flux1", "flux11", "flux3", "flux13",
        "fluxErr1", "fluxErr11", "fluxErr3", "fluxErr13"
    ]
):
    if not filename.endswith(".csv"):
        filename += ".csv"


    row = {
        "Run": run,
        "Fill": fill,
        "scale": scale,
        "nTracks1":      nTracks_dict.get(1,  0),
        "nTracks11":     nTracks_dict.get(11, 0),
        "nTracks3":      nTracks_dict.get(3,  0),
        "nTracks13":     nTracks_dict.get(13, 0),
        "trkeff1":       effs_dict.get(1,  0)[0],
        "trkeff11":      effs_dict.get(11, 0)[0],
        "trkeff3":       effs_dict.get(3,  0)[0],
        "trkeff13":      effs_dict.get(13, 0)[0],
        "trkeffErr1":    effs_dict.get(1,  0)[1],
        "trkeffErr11":   effs_dict.get(11, 0)[1],
        "trkeffErr3":    effs_dict.get(3,  0)[1],
        "trkeffErr13":   effs_dict.get(13, 0)[1],
        "flux1":         muon_flux_dict.get(1,  (0, 0))[0],
        "flux11":        muon_flux_dict.get(11, (0, 0))[0],
        "flux3":         muon_flux_dict.get(3,  (0, 0))[0],
        "flux13":        muon_flux_dict.get(13, (0, 0))[0],
        "fluxErr1":      muon_flux_dict.get(1,  (0, 0))[1],
        "fluxErr11":     muon_flux_dict.get(11, (0, 0))[1],
        "fluxErr3":      muon_flux_dict.get(3,  (0, 0))[1],
        "fluxErr13":     muon_flux_dict.get(13, (0, 0))[1]
    }

    file_exists = os.path.exists(filename)
    with open(filename, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=cols)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)




def save_muonFlux_to_root(
    filename: str,
    run: int,
    fill: int,
    muon_flux_dict: dict,
    cols: list[str] = ["Run", "Fill", "flux1", "flux11", "flux3", "flux13", "fluxErr1", "fluxErr11", "fluxErr3", "fluxErr13"],
    tree_name: str = "muonFlux"
):

    if not filename.endswith(".root"):
        filename = filename + ".root"

    def _entry_to_row(tree):
        return {
            "Run":       int(getattr(tree, cols[0])),
            "Fill":      int(getattr(tree, cols[1])),
            "flux1":     int(getattr(tree, cols[2])),
            "flux11":    int(getattr(tree, cols[3])),
            "flux3":     int(getattr(tree, cols[4])),
            "flux13":    int(getattr(tree, cols[5])),
            "fluxErr1":  int(getattr(tree, cols[6])),
            "fluxErr11": int(getattr(tree, cols[7])),
            "fluxErr3":  int(getattr(tree, cols[8])),
            "fluxErr13": int(getattr(tree, cols[9]))
        }

    existing_rows = []

    if os.path.exists(filename):
        f = ROOT.TFile(filename, "READ")
        existing_tree = f.Get(tree_name)
        if existing_tree:
            n_entries = int(existing_tree.GetEntries())
            for i in range(n_entries):
                existing_tree.GetEntry(i)
                existing_rows.append(_entry_to_row(existing_tree))
        f.Close()

    new_row = {
        "Run":       int(run),
        "Fill":      int(fill),
        "flux1":     int(muon_flux_dict.get(1,  (0, 0))[0]),
        "flux11":    int(muon_flux_dict.get(11, (0, 0))[0]),
        "flux3":     int(muon_flux_dict.get(3,  (0, 0))[0]),
        "flux13":    int(muon_flux_dict.get(13, (0, 0))[0]),
        "fluxErr1":  int(muon_flux_dict.get(1,  (0, 0))[1]),
        "fluxErr11": int(muon_flux_dict.get(11, (0, 0))[1]),
        "fluxErr3":  int(muon_flux_dict.get(3,  (0, 0))[1]),
        "fluxErr13": int(muon_flux_dict.get(13, (0, 0))[1])
    }
    existing_rows.append(new_row)

    outfile = ROOT.TFile(filename, "UPDATE")
    old_tree = outfile.Get(tree_name)
    if old_tree:
        outfile.Delete(f"{tree_name};*")
    tree = ROOT.TTree(tree_name, "Muon flux tree")

    run_buf     = array('i', [0])
    fill_buf    = array('i', [0])
    mf1_buf     = array('i', [0])
    mf11_buf    = array('i', [0])
    mf3_buf     = array('i', [0])
    mf13_buf    = array('i', [0])
    mferr1_buf  = array('i', [0])
    mferr11_buf = array('i', [0])
    mferr3_buf  = array('i', [0])
    mferr13_buf = array('i', [0])

    tree.Branch(cols[0], run_buf,  f"{cols[0]}/I")
    tree.Branch(cols[1], fill_buf, f"{cols[1]}/I")
    tree.Branch(cols[2], mf1_buf,  f"{cols[2]}/I")
    tree.Branch(cols[3], mf11_buf, f"{cols[3]}/I")
    tree.Branch(cols[4], mf3_buf,  f"{cols[4]}/I")
    tree.Branch(cols[5], mf13_buf, f"{cols[5]}/I")

    for row in existing_rows:
        run_buf[0]     = int(row["Run"])
        fill_buf[0]    = int(row["Fill"])
        mf1_buf[0]     = int(row["flux1"])
        mf11_buf[0]    = int(row["flux11"])
        mf3_buf[0]     = int(row["flux3"])
        mf13_buf[0]    = int(row["flux13"])
        mferr1_buf[0]  = int(row["fluxErr1"])
        mferr11_buf[0] = int(row["fluxErr11"])
        mferr3_buf[0]  = int(row["fluxErr3"])
        mferr13_buf[0] = int(row["fluxErr13"])
        tree.Fill()

    tree.Write()
    outfile.Close()

def write_output(
    filename: str,
    run: int, fill: int, scale: int,
    track_counts: dict,
    effs: dict,
    muon_flux: dict,
    tree_name: str = ""
):
    if not filename:
        return

    if filename.endswith(".csv"):
        save_all_to_csv(
            filename, run, fill, scale,
            track_counts["IP1"], effs, muon_flux
        )

    else:
        if not filename.endswith(".root"):
            filename += ".root"

        save_muonFlux_to_root(
            filename, run, fill, muon_flux, tree_name=tree_name
        )
