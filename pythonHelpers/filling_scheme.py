import os
import sys
from datetime import datetime, timezone

parts = os.environ["SNDSW_ROOT"].split('/')
sw_index = parts.index('sw')
base = "/".join(parts[:sw_index])
package = parts[sw_index + 2]
fs_source_path = os.path.join(base, package)
script_path = os.path.join(fs_source_path, 'shipLHC', 'scripts')
sys.path.insert(0, fs_source_path)

import FillingScheme as FS
import ROOT

import pythonHelpers.general

year_configs = {
    2022: {"rmin": 4360,  "has_run_subfolder": False},
    2023: {"rmin": 5412,  "has_run_subfolder": False},
    2024: {"rmin": 7648,  "has_run_subfolder": True},
    2025: {"rmin": 10918, "has_run_subfolder": True},
}

class FSOptions:
    def __init__(self,
        path: str = "./",
        raw_data_path: str = f"{os.environ['EOSSHIP']}/eos/experiment/sndlhc/raw_data/physics",
        www_path: str = f"{os.environ['EOSSHIP']}/eos/experiment/sndlhc/www/",
        lumi_version: str = "offline"
    ):
        self.path = "./"               # Required by the fillingScheme class
        self.rawData = raw_data_path
        self.www = www_path

        # These are needed by FS.Init, even if not directly used by getBunchStructureDict
        self.fillNumbers = ""
        self.runNumbers = ""
        self.command = ""
        self.lumiversion = lumi_version
        self.withIP2 = True
        self.nMin = 100000
        self.batch = True


def get_bunch_slots(
    input_files: str,
) -> dict:
    ch = pythonHelpers.general.load_snd_TChain(input_files)
    ch.GetEntry(0)
    fill = ch.EventHeader.GetFillNumber()
    run = ch.EventHeader.GetRunId()

    ts = ch.EventHeader.GetUTCtimestamp()
    ts_date = datetime.fromtimestamp(ts, tz=timezone.utc)
    year = ts_date.year

    eos_prefix = os.environ['EOSSHIP']
    raw_data_path = f"{eos_prefix}/eos/experiment/sndlhc/raw_data/physics/{year}"
    www_path = f"{eos_prefix}/eos/experiment/sndlhc/www/"

    options = FSOptions(raw_data_path=raw_data_path, www_path=www_path)
    config = year_configs.get(year)
    if config:
        options.rmin = config["rmin"]

        run_suffix = ""
        if config["has_run_subfolder"]:
            idx = options.rawData.find("run_")
            run_suffix = options.rawData[idx:] if idx != -1 else ""

        options.convpath = f"/eos/experiment/sndlhc/convertedData/physics/{year}/{run_suffix}"

    filling_scheme = FS.fillingScheme()
    filling_scheme.Init(options)

    bunch_slots = filling_scheme.getBunchStructureDict(fill)
    if not bunch_slots:
        print(f"Could not retrieve bunch dictionary for run {run} (fill {fill}).")
        return None

    return bunch_slots
