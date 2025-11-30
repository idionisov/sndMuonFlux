import argparse
import os

import ROOT

import pythonHelpers.general
import pythonHelpers.nTracks

this_dir = os.path.dirname(os.path.abspath(__file__))


def get_nTracks_pipeline(args):
    # Assumes the analysis is compiled in build directory
    path_to_ntrackslib = os.path.join(this_dir, ".", "build", "nTracks", "libnTracks.so")
    path_to_ntrackslib = os.path.abspath(path_to_ntrackslib)
    ROOT.gSystem.Load(path_to_ntrackslib)

    if len(args.fout) == 0:
        print("No output file specified. Tracking efficiency results will be printed without saving.")
        outnames = args.fout
    else:
        outnames = " ".join(args.fout)
    outfile_root, outfile_csv = pythonHelpers.general.get_outfiles(outnames)

    run  = pythonHelpers.general.get_snd_run(args.input_files)
    fill = pythonHelpers.general.get_lhc_fill(args.input_files)

    vec = ROOT.getNTracks(
        args.input_files,
        *args.x_range, *args.y_range,
        *args.xz_range, *args.yz_range,
        *args.z_ref
    )
    counts = pythonHelpers.nTracks.get_track_counts(vec, "IP1")
    counts.update({
        "Run": run, "Fill": fill, "scale": args.scale
    })
    print(counts)

    pythonHelpers.nTracks.write_output(outfile_root, run, fill, counts, args.scale)
    pythonHelpers.nTracks.write_output(outfile_csv,  run, fill, counts, args.scale)



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script for getting the number of tracks of a given run.")

    parser.add_argument('-i', '--input-files', type=str, required=True, help="Regex pattern for input ROOT files with reconstructed tracks, e.g., '/path/to/files*.root'.")
    parser.add_argument('-o', '--fout', type=str, nargs="+", default=[], help="Optional output files (root by default, but could be csv -- both formats simultaneously are supported).")
    parser.add_argument('-z', '--z-ref', nargs="+", type=float, default=[430., 430., 450., 450.], help="Reference z-plane coordinates for each track type (types 1, 11, 3, 13). Provide 4 numbers: zRef1 zRef11 zRef3 zRef13.")
    parser.add_argument('-x', '--x-range', nargs=2, type=float, default=[-42., -10.], help="Fiducial x-coordinate range [xmin, xmax] in cm for tracks to be counted.")
    parser.add_argument('-y', '--y-range', nargs=2, type=float, default=[19., 48.], help="Fiducial y-coordinate range [ymin, ymax] in cm for tracks to be counted.")
    parser.add_argument('-xz', '--xz-range', nargs=2, type=float, default=[-1e12, 1e12], help="Allowed track angle in XZ plane [xzMin, xzMax] in mrad. Tracks outside this range are ignored.")
    parser.add_argument('-yz', '--yz-range', nargs=2, type=float, default=[-1e12, 1e12], help="Allowed track angle in YZ plane [yzMin, yzMax] in mrad. Tracks outside this range are ignored.")
    parser.add_argument('-sc', '--scale', type=int, default=1, help="Scaling used for track reconstruction")

    args = parser.parse_args()


    get_nTracks_pipeline(args)
