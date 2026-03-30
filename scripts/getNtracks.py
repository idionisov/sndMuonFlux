import argparse
import os
import sys

# Add project root to path to allow importing from pythonHelpers
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import ROOT

import pythonHelpers.bunch_struct
import pythonHelpers.general
import pythonHelpers.nTracks

this_dir = os.path.dirname(os.path.abspath(__file__))


def get_nTracks_pipeline(args):
    print("Starting nTracks extraction")

    # Assumes the analysis is compiled in build directory
    pythonHelpers.general.load_cpp_extension("libnTracks.so", "nTracks")

    if len(args.fout) == 0:
        print("No output file specified. Track counts will be printed without saving.")
        outnames = args.fout
    else:
        outnames = " ".join(args.fout)
    outfile_root, outfile_csv = pythonHelpers.general.get_outfiles(outnames)
    # print(outfile_root)
    # print(outfile_csv)

    # run  = pythonHelpers.general.get_snd_run(args.input_files)
    # fill = pythonHelpers.general.get_lhc_fill(args.input_files)
    _, run, fill, acc_mode = pythonHelpers.general.load_run_info(args.input_files)

    start_ts = args.t_range[0] if args.t_range[0] is not None else -1.0
    end_ts   = args.t_range[1] if args.t_range[1] is not None else -1.0

    vec = ROOT.getNTracks(
        args.input_files,
        *args.x_range, *args.y_range,
        *args.xz_range, *args.yz_range,
        *args.z_ref,
        start_ts,
        end_ts
    )
    bunches = ("IP1", "IP2", "B1Only", "B2noB1", "IP2B1B2")
    counts = {
        b: pythonHelpers.nTracks.get_track_counts(vec, b) for b in bunches
    }

    pythonHelpers.nTracks.write_output(outfile_root, run, fill, counts["IP1"], args.scale)
    pythonHelpers.nTracks.write_output(outfile_csv,  run, fill, counts, args.scale)

    return counts


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script for getting the number of tracks of a given run.")

    parser.add_argument('-i', '--input-files', type=str, required=True, help="Regex pattern for input ROOT files with reconstructed tracks, e.g., '/path/to/files*.root'.")
    parser.add_argument('-t', '--t-range', nargs=2, type=float, default=[None, None], help="UTC timestamp range [start, end] in seconds for data selection.")
    parser.add_argument('-o', '--fout', type=str, nargs="+", default=[], help="Optional output files (root by default, but could be csv -- both formats simultaneously are supported).")
    parser.add_argument('-z', '--z-ref', nargs=4, type=float, default=[430., 430., 450., 450.], help="Reference z-plane coordinates for each track type (types 1, 11, 3, 13). Provide 4 numbers: zRef1 zRef11 zRef3 zRef13.")
    parser.add_argument('-x', '--x-range', nargs=2, type=float, default=[-42., -10.], help="Fiducial x-coordinate range [xmin, xmax] in cm for tracks to be counted.")
    parser.add_argument('-y', '--y-range', nargs=2, type=float, default=[19., 48.], help="Fiducial y-coordinate range [ymin, ymax] in cm for tracks to be counted.")
    parser.add_argument('-xz', '--xz-range', nargs=2, type=float, default=[-1e12, 1e12], help="Allowed track angle in XZ plane [xzMin, xzMax] in mrad. Tracks outside this range are ignored.")
    parser.add_argument('-yz', '--yz-range', nargs=2, type=float, default=[-1e12, 1e12], help="Allowed track angle in YZ plane [yzMin, yzMax] in mrad. Tracks outside this range are ignored.")
    parser.add_argument('-sc', '--scale', type=int, default=1, help="Scaling used for track reconstruction")

    args = parser.parse_args()


    counts = get_nTracks_pipeline(args)
    print(counts)
