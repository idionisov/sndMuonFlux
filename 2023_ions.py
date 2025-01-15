import os, re, ROOT


DATA = "/eos/experiment/sndlhc/convertedData/physics/2023"

print("SND@LHC run\tLHC fill\tEntries\tBunch Type")
for root, run_dirs, files in os.walk(DATA):
    for run_dir in run_dirs:
        #print(f"Subdirectory: {os.path.join(root, run_dir)}")
        
        dir = os.path.join(root, run_dir)
        tree = ROOT.TChain("cbmsim")
        tree.Add(f"{dir}/*.root")
        nentries = tree.GetEntries()
    
        runNr = int(re.search(r"run_(\d+)$", run_dir).group(1))
        tree.GetEntry(0)
        acc = tree.EventHeader.GetAccMode()
        fillNr = tree.EventHeader.GetFillNumber()
        
        if acc == 12:
            print(f"\033[1;32m{runNr}\t{fillNr}\t{nentries}\tPbPb\033[0m")
        else:
            print(f"\033[1;31m{runNr}\t{fillNr}\t{nentries}\tpp\033[0m")
