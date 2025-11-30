
import csv
import glob
import os
from array import array
from typing import Union

import ROOT


def has_tree(input_str: str, tree_name: str) -> bool:
    tmp_f = ROOT.TFile.Open(glob.glob(input_str)[0])
    if tmp_f.GetListOfKeys().Contains(tree_name):
        return True
    else:
        return False


def load_snd_TChain(input_str: str) -> ROOT.TChain:
    for tree_name in ("rawConv", "cbmsim"):
        if has_tree(input_str, tree_name):
            chain = ROOT.TChain(tree_name)
            chain.Add(input_str)
            return chain
    raise ValueError(f"No rawConv or cbmsim tree found in: {input_str}")


def get_snd_run(input_str: str) -> int:
    chain = load_snd_TChain(input_str)
    chain.GetEntry(0)
    return chain.EventHeader.GetRunId()

def get_lhc_fill(input_str: str) -> int:
    chain = load_snd_TChain(input_str)
    chain.GetEntry(0)
    return chain.EventHeader.GetFillNumber()

def load_run_info(input_str: str):
    ch = load_snd_TChain(input_str)
    ch.GetEntry(0)
    run = ch.EventHeader.GetRunId()
    fill = ch.EventHeader.GetFillNumber()
    acc_mode = ch.EventHeader.GetAccMode()
    return ch, run, fill, acc_mode


def compute_area(
    x_range: tuple[float, float],
    y_range: tuple[float, float]
) -> float:
    return abs(x_range[1] - x_range[0]) * abs(y_range[1] - y_range[0])

def get_outfiles(o: str = ""):
    if not o:
        return "", ""

    parts = o.split()

    if len(parts) > 2:
        raise ValueError("More than two output files were provided.")

    out_csv = ""
    out_root = ""

    for p in parts:
        if p.endswith(".csv"):
            if out_csv:
                raise ValueError("Multiple CSV output files provided.")
            out_csv = p

        elif p.endswith(".root"):
            if out_root:
                raise ValueError("Multiple ROOT output files provided.")
            out_root = p

        else:
            # If it's a single non-suffix file, assume it's .root
            if len(parts) == 1:
                out_root = p + ".root"
            else:
                raise ValueError(f"Unrecognized output file: {p}")

    return out_root, out_csv












def save_trkeff_to_csv(
    filename: str, run: int, fill: int, effs: dict,
    cols: list[str] = ["Run", "Fill", "trkeff1", "trkeffErr1", "trkeff11", "trkeffErr11", "trkeff3", "trkeffErr3", "trkeff13", "trkeffErr13"]
):
    if not filename.endswith(".csv"):
        filename += ".csv"


    row = {
        "Run": run,
        "Fill": fill,
        "trkeff1":      effs.get(1, 0)[0],
        "trkeffErr1":   effs.get(1, 0)[1],
        "trkeff11":     effs.get(11, 0)[0],
        "trkeffErr11":  effs.get(11, 0)[1],
        "trkeff3":      effs.get(3, 0)[0],
        "trkeffErr3":   effs.get(3, 0)[1],
        "trkeff13":     effs.get(13, 0)[0],
        "trkeffErr13":  effs.get(13, 0)[1]
    }

    file_exists = os.path.exists(filename)
    with open(filename, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=cols)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)



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
