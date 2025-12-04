import argparse
import time

import ROOT

import pythonHelpers.bunch_struct

if __name__=="__main__":
    start = time.time()

    parser = argparse.ArgumentParser(
        description="Extract bunch structure from input files. Adapted to work with both 25 ns (protons) and 50 ns (heavy-ions) schemes."
    )
    parser.add_argument("-i", "--input-files", type=str, required=True, help="Input files.")
    parser.add_argument("-o", "--fout", type=str, default="", help="Optional output file to store bunch structure histograms.")
    parser.add_argument("--n-break", type=int, default=1e6, help="Breakpoint for event loop.")
    args = parser.parse_args()



    bunch_counts = pythonHelpers.bunch_struct.get_bunch_counts(
        args.input_files,
        args.fout,
        args.n_break
    )
    print(bunch_counts)

    elapsed = time.time() - start
    m, s = divmod(elapsed, 60)
    h, m = divmod(m, 60)
    print(f"Runtime: {int(h)}h {int(m)}m {int(s)}s  ({elapsed:.02f} seconds)")
