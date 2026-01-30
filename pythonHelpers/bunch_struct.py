import glob
import re
from datetime import datetime, timezone
from typing import Union

import FillingScheme as FS
import ROOT

import pythonHelpers.general
from pythonHelpers.filling_scheme import get_bunch_slots


def get_bunch_spacing(
    arg: Union[int, str, None] = None
) -> int:
    """
    Extract the bunch spacing for
    - A given filling scheme string OR
    - Given input root files containing 'rawConv' or 'cbmsim' tree with an SNDLHCEventHeader OR
    - A given SNDLHCEventHeader acc mode (integer)
    """
    def from_acc_mode(acc_mode: int) -> int:
        acc_modes = {
            # Constants taken from https://snd-lhc.github.io/sndsw/SNDLHCEventHeaderConst_8h_source.html
            # IonPhysics: 12 and ProtonPhysics: 11
            11: 25, 12: 50,
        }

        if acc_mode in acc_modes:
            return acc_modes[acc_mode]
        else:
            raise ValueError(f"Cannot extract bunch spacing for acc mode {arg}.")


    if isinstance(arg, int):
        return from_acc_mode(arg)

    elif isinstance(arg, str):
        if "*" in arg or "?" in arg:
            files = glob.glob(arg)
            if not files:
                raise ValueError(f"No files match wildcard: {arg}")
            arg = files[0]

        if arg.endswith(".root"):
            ch = pythonHelpers.general.load_snd_TChain(arg)
            ch.GetEntry(0)
            acc_mode = ch.EventHeader.GetAccMode()
            return from_acc_mode(acc_mode)

        else:
            m = re.search(r'(\d+)ns', arg)
            if m:
                return int(m.group(1))
            else:
                raise ValueError(f"Cannot extract bunch spacing from scheme {arg}.")

    else:
        raise ValueError("No argument provided.")



def get_slot_numbers(
    arg: Union[int, str, None] = None,
    harmonic_number: int = 35640,
    bucket_spacing: float = 2.5
) -> int:
    """
    Extract the bunch spacing for
    - A given filling scheme string OR
    - Given input root files containing 'rawConv' or 'cbmsim' tree with an SNDLHCEventHeader OR
    - A given SNDLHCEventHeader acc mode (integer)
    """
    # LHC has a harmonic number h = 35640
    # (number of oscillations of the RF electromagnetic field that fit
    # exactly into one full turn of the ring)
    # That means we have one bucket every T = T_rev / h ~ 2.5 ns
    # Nominal bunch spacing is a multiple of that 2.5 ns spacing
    # For protons the scaling factor is S=10 to get 10*2.5=25 ns
    # For ions it is S=20 to get 20*2.5=50 ns
    # Number of slots are N = (h / nominal_spacing * S)
    bunch_spacing = get_bunch_spacing(arg)
    return int(harmonic_number * bucket_spacing / bunch_spacing)


def extract_bunch_struct(
    input_files: str,
    fout: str = "",
    n_break: int = 1e6
) -> dict:

    print("Loading TChain...", end="")
    ch = pythonHelpers.general.load_snd_TChain(input_files)
    n_entries = ch.GetEntries()
    print(f" {n_entries} entries.")

    ch.GetEntry(0)
    run = ch.EventHeader.GetRunId()
    fill = ch.EventHeader.GetFillNumber()
    acc_mode = ch.EventHeader.GetAccMode()
    bunch_slots = get_slot_numbers(acc_mode)

    if fout:
        fout = ROOT.TFile(fout, "update")
        run_dir = fout.GetDirectory(f"Run{run}")
        if not run_dir:
            run_dir = fout.mkdir(f"Run{run}")
        run_dir.cd()

    print(f"Starting bunch structure extraction for run {run} (fill {fill})")
    h = {
        'IP1': ROOT.TH1D(f"h_IP1_run{run}", "IP1;Bunch number;", bunch_slots, 0, bunch_slots),
        'IP2': ROOT.TH1D(f"h_IP2_run{run}", "IP2;Bunch number;", bunch_slots, 0, bunch_slots),
        'B1':  ROOT.TH1D(f"h_B1_run{run}",  "B1;Bunch number;",  bunch_slots, 0, bunch_slots),
        'B2':  ROOT.TH1D(f"h_B2_run{run}",  "B2;Bunch number;",  bunch_slots, 0, bunch_slots),

        'B1Only':  ROOT.TH1D(f"h_B1Only_run{run}",  "B1Only;Bunch number;",  bunch_slots, 0, bunch_slots),
        'B2noB1':  ROOT.TH1D(f"h_B2noB1_run{run}",  "B2noB1;Bunch number;",  bunch_slots, 0, bunch_slots),

        'IP1&B1':  ROOT.TH1D(f"h_IP1&B1_run{run}",  "IP1&B1;Bunch number;",  bunch_slots, 0, bunch_slots),
        'IP1&B2':  ROOT.TH1D(f"h_IP1&B2_run{run}",  "IP1&B2;Bunch number;",  bunch_slots, 0, bunch_slots),
    }
    h['IP1'].SetLineColor(ROOT.kBlue)
    h['B1'].SetLineColor(ROOT.kRed)
    h['B2'].SetLineColor(ROOT.kCyan+1)
    h['B2'].SetLineWidth(2)
    h['IP2'].SetLineColor(ROOT.kOrange)
    h['IP2'].SetLineWidth(2)

    h["B1Only"].SetLineColor(ROOT.kYellow+3)
    h["B2noB1"].SetLineColor(ROOT.kRed-5)
    h["B1Only"].SetLineWidth(2)
    h["B2noB1"].SetLineWidth(2)

    h["IP1&B1"].SetLineColor(ROOT.kMagenta+3)
    h["IP1&B2"].SetLineColor(ROOT.kGreen-6)
    h["IP1&B1"].SetLineWidth(2)
    h["IP1&B2"].SetLineWidth(2)


    b = {'IP1': [], 'IP2': [], 'B1': [], 'B2': [], "B1Only": [], "B2noB1": [], "IP1&B1": [], "IP1&B2": []}
    # b = {'IP1': [], 'IP2': [], 'B1': [], 'B2': []}

    h_trackRate = None
    if ch.GetListOfBranches().FindObject("Reco_MuonTracks"):
        h_trackRate = ROOT.TH1D("h_trackRate", f"Track rate (run {run});Bunch number;", bunch_slots+1, 0, bunch_slots)
        h_trackRate.SetLineColor(1)
        h_trackRate.SetFillColor(14)
        h_trackRate.SetLineWidth(1)

    h_eventRate = ROOT.TH1D("h_eventRate", f"Event rate (run {run});Bunch number;", bunch_slots+1, 0, bunch_slots)
    h_eventRate.SetLineColor(1)
    h_eventRate.SetFillColor(14)
    h_eventRate.SetLineWidth(1)
    print("Successfully created histograms.")

    print("Entering event loop...")
    n_break = min(n_break, n_entries)
    next_print = 0
    for i_event, event in enumerate(ch):
        # Print progress every 5%
        progress = (i_event * 100) // n_break
        if progress >= next_print:
            print(f"{progress:.0f} %")
            next_print += 5

        # Breakpoint to prevent looping through too many events
        if i_event >= n_break:
            break


        eh = event.EventHeader
        eventTime = eh.GetEventTime()
        bunch_phase = eventTime % (8*bunch_slots)/8+0.5
        bunch_nr = int( bunch_phase )

        h_eventRate.Fill(bunch_phase)
        if eh.isIP1() and bunch_nr not in b['IP1']: b['IP1'].append(bunch_nr)
        if eh.isIP2() and bunch_nr not in b['IP2']: b['IP2'].append(bunch_nr)
        if eh.isB1()  and bunch_nr not in b['B1']:  b['B1'].append(bunch_nr)
        if eh.isB2()  and bunch_nr not in b['B2']:  b['B2'].append(bunch_nr)

        if eh.isB1Only() and bunch_nr not in b['B1Only']:  b['B1Only'].append(bunch_nr)
        if eh.isB2noB1() and bunch_nr not in b['B2noB1']:  b['B2noB1'].append(bunch_nr)

        if eh.isIP1() and eh.isB1() and bunch_nr not in b['IP1&B1']: b['IP1&B1'].append(bunch_nr)
        if eh.isIP1() and eh.isB2() and bunch_nr not in b['IP1&B2']: b['IP1&B2'].append(bunch_nr)

        if h_trackRate and event.Reco_MuonTracks:
            for trk in event.Reco_MuonTracks:
                if not (
                    trk.getTrackFlag() and
                    trk.getTrackMom().Z()>0
                ): continue

                h_trackRate.Fill(bunch_phase)
                break


    for i in "IP1", "IP2", "B1", "B2", "B1Only", "B2noB1", "IP1&B1", "IP1&B2":
        for bunch_nr in range(bunch_slots):
            if bunch_nr in b[i]:
                h[i].Fill(bunch_nr)

        # Delete the ones needed just to fill histograms
        # Keeps output cleaner
        if i in ("B1Only", "B2noB1", "IP1&B1", "IP1&B2"):
            del b[i]

    er_max  = h_eventRate.GetMaximum()
    ip1_max = 2.10*er_max
    b1_max  = 1.40*er_max
    b2_max  = 0.70*er_max
    ip2_max = 0.35*er_max

    h['IP1'].Scale(ip1_max)
    h['IP2'].Scale(ip2_max)
    h['B1'].Scale(b1_max)
    h['B1Only'].Scale(1.05*er_max)
    h['B2'].Scale(b2_max)
    h['B2noB1'].Scale(0.85*er_max)
    h['IP1&B1'].Scale(1.65*er_max)
    h['IP1&B2'].Scale(1.25*(er_max))

    if fout:
        c_events = ROOT.TCanvas(f"c_BunchStruct_events_Run{run}", f"Bunch Structure of Run {run}", 900, 700)
        h_eventRate.SetMaximum(1.025*ip1_max)
        h_eventRate.Draw("hist")
        for i in "IP1", "B1", "B2", "IP2":
            h[i].Draw("hist same")
        c_events.Write()
        h_eventRate.Write()
        if h_trackRate:
            h_trackRate.Write()
            c_tracks = ROOT.TCanvas(f"c_BunchStruct_tracks_Run{run}", f"Bunch Structure of Run {run}", 900, 700)
            h_trackRate.SetMaximum(1.025*ip1_max)
            h_trackRate.Draw("hist")
            for i in "IP1", "B1", "B2", "IP2":
                h[i].Draw("hist same")
                h[i].Write()
            c_tracks.Write()

        for i in "B1Only", "B2noB1", "IP1&B1", "IP1&B2":
            h[i].Write()

        fout.Close()
        print(f"Saved bunch structured in: {fout.GetName()}.")
    return b



def get_bunch_slots(
    input_files: str,
) -> dict:
    ch = pythonHelpers.general.load_snd_TChain(input_files)
    fill = ch.EventHeader.GetFillNumber()
    run = ch.EventHeader.GetRunId()
    ch.GetEntry(0)

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


def get_bunch_counts(
    bunch_slots: dict,
    include: tuple = (),
    exclude: tuple = ()
):
    if not include:
        return 0

    result = set(bunch_slots[include[0]])
    for key in include[1:]:
        result &= set(bunch_slots[key])

    for key in exclude:
        if key in bunch_slots:
            result -= set(bunch_slots[key])

    return len(result)
