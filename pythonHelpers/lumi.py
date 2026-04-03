import glob
import os
from typing import Optional, Union

import numpy as np
import ROOT

sndsw_path = os.environ["SNDSW_ROOT"]
ROOT.gInterpreter.ProcessLine(
    f'#include "{sndsw_path}/analysis/tools/sndTchainGetter.h"'
)

import pythonHelpers.general


def get_fill_lumi_path(
    arg: Union[ROOT.TTree, ROOT.TChain, int, str, None] = None,
    run: Optional[int] = None,
) -> str:
    if arg is None and run is None:
        raise ValueError("Provide a fill, input string, TTree or run number!")

    def get_fill_lumi_path_from_fill(fill: int) -> str:
        return f"/eos/experiment/sndlhc/atlas_lumi/fill_{fill:06d}.root"

    def get_fill_lumi_path_from_tchain(ch: Union[ROOT.TTree, ROOT.TChain]) -> str:
        if ch.GetEntries() == 0:
            raise ValueError(
                "The provided TChain/TTree is empty. Cannot extract fill number."
            )

        ch.GetEntry(0)
        fill = ch.EventHeader.GetFillNumber()
        return get_fill_lumi_path_from_fill(fill)

    def get_fill_lumi_path_from_run(run: int) -> str:
        ch = ROOT.snd.analysis_tools.GetTChain(run)
        return get_fill_lumi_path_from_tchain(ch)

    def get_fill_lumi_path_from_files(input_files: str) -> str:
        files = glob.glob(input_files)
        if not files:
            raise ValueError(f"No files found matching: {input_files}")

        tmp_f = ROOT.TFile.Open(files[0])
        if not tmp_f or tmp_f.IsZombie():
            raise IOError(f"Failed to open ROOT file: {files[0]}")

        has_lumi_branch = tmp_f.GetListOfKeys().Contains("atlas_lumi")
        tmp_f.Close()
        if has_lumi_branch:
            return files[0]

        ch = pythonHelpers.general.load_snd_TChain(input_files)
        return get_fill_lumi_path_from_tchain(ch)

    if run is not None:
        return get_fill_lumi_path_from_run(run)

    if isinstance(arg, ROOT.TTree):
        return get_fill_lumi_path_from_tchain(arg)
    elif isinstance(arg, str):
        return get_fill_lumi_path_from_files(arg)
    elif isinstance(arg, int):
        return get_fill_lumi_path_from_fill(arg)
    else:
        raise TypeError(f"Unsupported argument type for arg: {type(arg)}")


def get_lumi_eos(
    arg: Optional[Union[str, int]] = None,
    run: Optional[int] = None,
    start_ts: Optional[float] = None,
    end_ts: Optional[float] = None,
) -> float:
    atlas_lumi_path = get_fill_lumi_path(arg=arg, run=run)
    atlas_lumi = ROOT.TChain("atlas_lumi")
    atlas_lumi.Add(atlas_lumi_path)

    delivered_inst_lumi = []
    delivered_unix_timestamp = []

    for entry in atlas_lumi:
        delivered_inst_lumi.append(entry.var)
        delivered_unix_timestamp.append(entry.unix_timestamp)

    delivered_inst_lumi = np.array(delivered_inst_lumi)
    delivered_unix_timestamp = np.array(delivered_unix_timestamp)

    if start_ts is not None:
        mask_start = delivered_unix_timestamp >= start_ts
        delivered_inst_lumi = delivered_inst_lumi[mask_start]
        delivered_unix_timestamp = delivered_unix_timestamp[mask_start]

    if end_ts is not None:
        mask_end = delivered_unix_timestamp <= end_ts
        delivered_inst_lumi = delivered_inst_lumi[mask_end]
        delivered_unix_timestamp = delivered_unix_timestamp[mask_end]

    delivered_deltas = delivered_unix_timestamp[1:] - delivered_unix_timestamp[:-1]
    delivered_mask = delivered_deltas < 600
    delivered_run = delivered_mask

    if not np.any(delivered_run):
        return 0.0

    return (
        np.cumsum(
            delivered_deltas[delivered_run] * delivered_inst_lumi[1:][delivered_run]
        )[-1]
        / 1e3
    )
