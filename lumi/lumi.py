import argparse
from sndUtils import SndData
import matplotlib.pyplot as plt
import datetime
import dateutil.parser as dp
import time
import ROOT
import numpy as np
import json

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--run', type=int, default=7080)
args = parser.parse_args()

run = args.run

def makeUnixTime(year, month, day, hour, minute, second) :
    dt = datetime.datetime(year, month, day, hour, minute, second)
    return time.mktime(dt.timetuple())

atlas_online_lumi = ROOT.TChain("atlas_lumi")
data = SndData(Run=run)

plt.style.use("root")
# path to luminosity file
# one can extract the LHC fill number from the SNDLHCEventHeader.
input_dir = ""
atlas_online_lumi.Add(f"/eos/experiment/sndlhc/atlas_lumi/fill_{data.Fill:06d}.root")

delivered_inst_lumi = []
delivered_unix_timestamp = []
delivered_run_number = []
delivered_fill_number = []
# start timestamp - in case our snd run starts during Stable beams
# must be 0 for almost all runs
fill = 0

for entry in atlas_online_lumi :
    delivered_inst_lumi.append(entry.var)
    delivered_unix_timestamp.append(entry.unix_timestamp)

recorded_mask = np.array(True)
delivered_inst_lumi = np.array(delivered_inst_lumi)
delivered_unix_timestamp = np.array(delivered_unix_timestamp)
delivered_deltas = delivered_unix_timestamp[1:] - delivered_unix_timestamp[:-1]
delivered_mask = delivered_deltas < 600 # related to sampling of tstamps in the L input file-> typically don't touch
delivered_run = np.logical_and(delivered_unix_timestamp[1:] > fill, delivered_mask)



plt.figure(figsize=(10,6))
plt.plot(delivered_unix_timestamp, delivered_inst_lumi)
plt.xlabel('Unix Timestamp')
plt.ylabel('Delivered Instantaneous Luminosity')
plt.title('ATLAS Luminosity vs Time')
plt.grid(True)
plt.show()


print(len(delivered_unix_timestamp))
