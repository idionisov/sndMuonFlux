
import glob

import numpy as np
import ROOT

import pythonHelpers.general


def get_fill_lumi_path(input_files: str) -> str:
    files = glob.glob(input_files)
    if not files:
        raise ValueError(f"No files found matching: {input_files}")

    tmp_f = ROOT.TFile.Open(files[0])
    if tmp_f.GetListOfKeys().Contains("atlas_lumi"):
        return files[0]
    del tmp_f

    ch = pythonHelpers.general.load_snd_TChain(input_files)
    ch.GetEntry(0)
    fill = ch.EventHeader.GetFillNumber()
    return f"/eos/experiment/sndlhc/atlas_lumi/fill_{fill:06d}.root"




def get_lumi_eos(input_files: str) -> float:
    files = glob.glob(input_files)
    if not files:
        raise ValueError(f"No files found matching: {input_files}")

    atlas_lumi_path = get_fill_lumi_path(files[0])
    atlas_lumi = ROOT.TChain("atlas_lumi")
    atlas_lumi.Add(atlas_lumi_path)

    delivered_inst_lumi = []
    delivered_unix_timestamp = []
    for entry in atlas_lumi:
        delivered_inst_lumi.append(entry.var)
        delivered_unix_timestamp.append(entry.unix_timestamp)

    delivered_inst_lumi = np.array(delivered_inst_lumi)
    delivered_unix_timestamp = np.array(delivered_unix_timestamp)

    delivered_deltas = delivered_unix_timestamp[1:] - delivered_unix_timestamp[:-1]
    delivered_mask = delivered_deltas < 600
    delivered_run = delivered_mask

    return np.cumsum(delivered_deltas[delivered_run] * delivered_inst_lumi[1:][delivered_run])[-1] / 1e3
