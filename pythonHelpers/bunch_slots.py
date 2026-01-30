import os
import pprint
import sys

import ROOT

# --- Configuration (Adjust these paths and values as needed) ---
# Ensure the directory containing FillingScheme.py is in your Python path
# If not, you might need to add it, e.g.:
sys.path.insert(0, '/afs/cern.ch/user/i/idioniso/snd_FS/sndsw/shipLHC/scripts')

# Set EOSSHIP environment variable if not already set in your session
if 'EOSSHIP' not in os.environ:
    os.environ['EOSSHIP'] = '/eos/experiment/sndlhc' # Adjust to your actual EOSSHIP path

# The fill number you want to analyze
TARGET_FILL_NUMBER = "8786"
# Path where intermediate ROOT files (e.g., fillingScheme-XXXX.root) will be stored
OUTPUT_PATH = "./" # Current directory, or specify a writable path
# Raw data path, used to infer convpath, rmin, etc.
RAW_DATA_PATH = "/eos/experiment/sndlhc/raw_data/physics/2023" # Adjust for the correct year
# Path to the WWW directory, used by some functions to find offline data
WWW_PATH = os.environ['EOSSHIP'] + "/eos/experiment/sndlhc/www/"
# -----------------------------------------------------------------

# Import the FillingScheme module
import FillingScheme as FS


# Create a dummy options object to mimic argparse arguments
class Options:
    def __init__(self):
        self.path = OUTPUT_PATH
        self.rawData = RAW_DATA_PATH
        self.www = WWW_PATH
        # These are needed by FS.Init, even if not directly used by getBunchStructureDict
        self.fillNumbers = ""
        self.runNumbers = ""
        self.command = ""
        self.lumiversion = "offline"
        self.withIP2 = True
        self.nMin = 100000
        self.batch = False # Set to True if you don't want ROOT canvases to pop up

# Instantiate options
options = Options()

# Replicate the logic from FillingScheme.py's __main__ block to set
# options.convpath and options.rmin based on rawData
if options.rawData.find('2022') > 0:
    options.convpath = "/eos/experiment/sndlhc/convertedData/physics/2022/"
    options.rmin = 4361 - 1
elif options.rawData.find('2023') > 0:
    options.convpath = "/eos/experiment/sndlhc/convertedData/physics/2023/"
    options.rmin = 5413 - 1
elif options.rawData.find('2024') > 0:
    em_run_start_idx = options.rawData.find("run_")
    em_run = options.rawData[em_run_start_idx:] if em_run_start_idx != -1 else ""
    options.convpath = "/eos/experiment/sndlhc/convertedData/physics/2024/" + em_run
    options.rmin = 7649 - 1
elif options.rawData.find('2025') > 0:
    em_run_start_idx = options.rawData.find("run_")
    em_run = options.rawData[em_run_start_idx:] if em_run_start_idx != -1 else ""
    options.convpath = "/eos/experiment/sndlhc/convertedData/physics/2025/" + em_run
    options.rmin = 10919 - 1

# Create a fillingScheme object
filling_scheme_analyzer = FS.fillingScheme()
filling_scheme_analyzer.Init(options)

# Get the bunch dictionary
bunch_dictionary = filling_scheme_analyzer.getBunchStructureDict(TARGET_FILL_NUMBER)

if bunch_dictionary:
    print(f"\nBunch structure for fill {TARGET_FILL_NUMBER}:")
    pprint.pprint(bunch_dictionary)
    # You can now use bunch_dictionary for further analysis
    # For example:
    # print(f"\nMax bunch number in IP1: {max(bunch_dictionary['IP1'])}")
else:
    print(f"Could not retrieve bunch dictionary for fill {TARGET_FILL_NUMBER}.")
