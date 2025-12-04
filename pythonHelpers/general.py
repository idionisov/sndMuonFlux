
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
