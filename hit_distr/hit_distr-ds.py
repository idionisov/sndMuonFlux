import argparse
from typing import Union
from ROOT import TFile, TH1I, TChain, TFile, MuFilterHit
from time import time
from ddf.root.misc import save_to_root
from ddf.snd.snd import get_N
from ddf.snd.trk import get_intersection
from ddf.pyfuncs import get_current_dmy, print_status

parser = argparse.ArgumentParser()
parser.add_argument('--run', type=int, default=7080)
parser.add_argument('--selection', type=str, default="*f10*")
parser.add_argument('--input_dir', type=str, default="")
parser.add_argument('--output_dir', type=str, default="")

args = parser.parse_args()

run = args.run
selection = args.selection

if not args.input_dir:
    input_dir = f"/eos/user/i/idioniso/1_Data/Tracks/run_{run:06d}"
else:
    input_dir = args.input_dir

if not args.output_dir:
    output_dir = "/eos/user/i/idioniso/mfout"
else:
    output_dir = args.input_dir

input_files = f"{input_dir}/{selection}.root"
fout_name = f"hitDistr_run{run}.root"

cbmsim = TChain("cbmsim")
cbmsim.Add(input_files)
nevents = cbmsim.GetEntries()

fout = TFile(f"{output_dir}/{fout_name}", "recreate")


h = {
    'v': {pl: TH1I(f"h_vDS{pl}-hits_run{run}", f"DS hit distribution, vDS{pl} (run {run});Vertical DS bar;", 60, 0, 61) for pl in range(1, 5)},
    'h': {pl: TH1I(f"h_hDS{pl}-hits_run{run}", f"DS hit distribution, hDS{pl} (run {run});Horizontal DS bar;", 60, 0, 61) for pl in range(1, 4)}
}
for i in h.values():
    for hist in i.values():
        hist.SetLineColor(1)


def get_dsBar(mf_hit: MuFilterHit) -> int:
    if mf_hit.GetSystem() != 3: raise ValueError("Not a DS MuFilterHit!")

    if mf_hit.isVertical():
        return mf_hit.GetDetectorID()%1000-60
    else:
        return mf_hit.GetDetectorID()%1000

start_time = time()
count = 0

for i_event, event in enumerate(cbmsim):
    count = print_status(i_event, nevents, start_time, count)

    for mf_hit in event.Digi_MuFilterHits:
        if mf_hit.GetSystem() != 3: continue
        
        bar = get_dsBar(mf_hit)
        pl = mf_hit.GetPlane()+1
        
        if mf_hit.isVertical():
            h['v'][pl].Fill(bar)
        else:
            h['h'][pl].Fill(bar)

save_to_root(h, tfile=fout)
fout.Close()
