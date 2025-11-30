import glob
import re
from typing import Union

import ROOT

import pythonHelpers.general


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



def get_bunch_slots(
    arg: Union[int, str, None] = None
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
    return int(35640*2.5/bunch_spacing)


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
    bunch_slots = get_bunch_slots(acc_mode)

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


def get_bunch_counts(
    input_files: str,
    output_file: str,
    n_break: int = 1e7,
) -> dict:
    bunch_slots = extract_bunch_struct(
        input_files,
        output_file,
        n_break
    )


    return {
        'IP1': len(bunch_slots["IP1"]),
        'IP2': len(bunch_slots["IP2"]),
        'B1': len(bunch_slots["B1"]),
        'B2': len(bunch_slots["B2"]),

        # IP2, B1 and B2 almost never appear individually in ion runs
        # ---> Run 7080 isolated bunches: IP2=3, B1=3, B2=0
        #      IP2 and B1 and B2 = 125
        # Since they're always together their contributions to non-IP1 muons are added together
        'IP2&B1&B2': len(
            set(bunch_slots["B2"])
            & set(bunch_slots["B1"])
            & set(bunch_slots["B2"])
            - set(bunch_slots["IP1"])
        )
    }
