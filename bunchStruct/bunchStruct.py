import argparse
from ROOT import TChain, TFile, TCanvas, TH1F, kBlue, kRed, kCyan, kOrange

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--run", type=int, default=7268)
    parser.add_argument("-f", "--files", type=str, default="sndsw_raw*-000*")
    args = parser.parse_args()

    run = args.run
    selection = args.selection

    t = TChain("rawConv")
    t.Add(f"$SND_DATA/2023_reprocess/run_{run:06d}/{selection}.root")
    nevents=t.GetEntries()
    print(f"Events: {nevents:,}")


    fout = TFile(f"/eos/user/i/idioniso/mfout/bunch_struct_{run}.root", "recreate")

    h = {
        'IP1': TH1F(f"h_IP1_run{run}", "IP1;Bunch number;", 1782, 0, 1782),
        'IP2': TH1F(f"h_IP2_run{run}", "IP2;Bunch number;", 1782, 0, 1782),
        'B1':  TH1F(f"h_B1_run{run}",  "B1;Bunch number;",  1782, 0, 1782),
        'B2':  TH1F(f"h_B2_run{run}",  "B2;Bunch number;",  1782, 0, 1782),
    }
    h['IP1'].SetLineColor(kBlue)
    h['B1'].SetLineColor(kRed)
    h['B2'].SetLineColor(kCyan+1)
    h['B2'].SetLineWidth(2)
    h['IP2'].SetLineColor(kOrange)
    h['IP2'].SetLineWidth(2)


    b = {'IP1': [], 'IP2': [], 'B1': [], 'B2': []}

    h_eventRate = TH1F("h_eventRate", f"Event rate (run {run});Bunch number;", 1782, 0, 1782)
    h_eventRate.SetLineColor(1)
    h_eventRate.SetFillColor(14)
    h_eventRate.SetLineWidth(1)

    for i_event, event in enumerate(t):
        if i_event%100000==0: print(f"Run {run}:\t{i_event:,}/{nevents:,}")
        eh = event.EventHeader
        eventTime = eh.GetEventTime()
        _eRate = eventTime % (8*1782)/8+0.5
        bunchNr = int( _eRate )

        h_eventRate.Fill(_eRate)

        if eh.isIP1() and bunchNr not in b['IP1']: b['IP1'].append(bunchNr)
        if eh.isIP2() and bunchNr not in b['IP2']: b['IP2'].append(bunchNr)
        if eh.isB1()  and bunchNr not in b['B1']:  b['B1'].append(bunchNr)
        if eh.isB2()  and bunchNr not in b['B2']:  b['B2'].append(bunchNr)

    for i in "IP1", "IP2", "B1", "B2":
        for bunchNr in range(1782):
            if bunchNr in b[i]: h[i].Fill(bunchNr)

    er_max  = h_eventRate.GetMaximum()
    ip1_max = 2.10*er_max
    b1_max  = 1.40*er_max
    b2_max  = 0.70*er_max
    ip2_max = 0.35*er_max

    h['IP1'].Scale(ip1_max)
    h['IP2'].Scale(ip2_max)
    h['B1'].Scale(b1_max)
    h['B2'].Scale(b2_max)


    c = TCanvas(f"c_BunchStruct_Run{run}", f"Bunch Structure of Run {run}", 900, 700)
    h_eventRate.SetMaximum(1.025*ip1_max)
    h_eventRate.Draw("hist")
    for i in "IP1", "B1", "B2", "IP2":
        h[i].Draw("hist same")

    c.Write()

    fout.Write()
    fout.Close()


if __name__=="__main__":
    main()
