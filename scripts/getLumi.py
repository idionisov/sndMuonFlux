import os
import sys

# Add project root to path to allow importing from pythonHelpers
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pythonHelpers.lumi

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input_files", type=str, help="Input files for luminosity extraction.")
    args = parser.parse_args()

    lumi = pythonHelpers.lumi.get_lumi_eos(args.input_files)
    print(lumi)
