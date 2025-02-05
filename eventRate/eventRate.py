import argparse
import ROOT
from time import time
from ddfUtils import printStatus
from sndUtils import SndData


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run', type=int, default=7080)
    parser.add_argument('-f', '--files', type=str, default="*.root")
    args = parser.parse_args()

    run = args.run
    files = args.files
    mfout = "/eos/user/i/idioniso/mfout"

    data = SndData(Run=args.run, InputDir="/eos/user/i/idioniso/1_Data/Tracks", Files=files)
    data.Print()

    data.Tree.GetEntry(0)
    tStart = data.EventHeader.GetUTCtimestamp()
    data.Tree.GetEntry(data.Tree.GetEntries()-1)
    tEnd = data.EventHeader.GetUTCtimestamp()
    runDuration = (tEnd - tStart)/3600
    print(f" >> Run duration: {runDuration} h")


    fout = ROOT.TFile(f"{mfout}/eventRateRun{run}.root", "recreate")

    h = ROOT.TH1D(
        f"h_eventRate_Run{run}",
        f"Event Rate (run {run});Time in seconds after start DAQ run;Number of events",
        1022, 0, runDuration
    )

    t0 = time()
    count = 0
    for i_event, event in enumerate(data.Tree):
        count = printStatus(i_event, event, t0, count)

        t = event.EventHeader.GetUTCtimestamp() - tStart
        h.Fill(t)

    fout.Write()
    fout.Close()


if __name__=="__main__":
    main()
