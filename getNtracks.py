import argparse
import os

import ROOT

import pythonHelpers.bunch_struct
import pythonHelpers.general
import pythonHelpers.nTracks

this_dir = os.path.dirname(os.path.abspath(__file__))


def get_nTracks_pipeline(args):
    print("Starting nTracks extraction")

    # Assumes the analysis is compiled in build directory
    path_to_ntrackslib = os.path.join(this_dir, ".", "build", "nTracks", "libnTracks.so")
    path_to_ntrackslib = os.path.abspath(path_to_ntrackslib)
    ROOT.gSystem.Load(path_to_ntrackslib)
    print(hasattr(ROOT, "getNTracks"))

    if len(args.fout) == 0:
        print("No output file specified. Tracking efficiency results will be printed without saving.")
        outnames = args.fout
    else:
        outnames = " ".join(args.fout)
    outfile_root, outfile_csv = pythonHelpers.general.get_outfiles(outnames)
    print(outfile_root)
    print(outfile_csv)

    # run  = pythonHelpers.general.get_snd_run(args.input_files)
    # fill = pythonHelpers.general.get_lhc_fill(args.input_files)
    _, run, fill, acc_mode = pythonHelpers.general.load_run_info(args.input_files)

    vec = ROOT.getNTracks(
        args.input_files,
        *args.x_range, *args.y_range,
        *args.xz_range, *args.yz_range,
        *args.z_ref
    )
    bunches = ("IP1", "IP2", "B1Only", "B2noB1", "IP2B1B2")
    counts = {
        b: pythonHelpers.nTracks.get_track_counts(vec, b) for b in bunches
    }
    output = counts.copy()
    output.update({
        "Run": run, "Fill": fill, "scale": args.scale
    })
    for b in bunches:
        suffix = f"_{b}" if b != "IP1" else ""
        for tt, ntrks in output[b].items():
            output.update({f"nTracks{tt}{suffix}": ntrks})

    if args.bunch_correction:
        bunch_slots = pythonHelpers.bunch_struct.extract_bunch_struct(
            args.input_files,
        )

        if acc_mode==11: # Protons
            N_IP1_and_B1 = pythonHelpers.bunch_struct.get_bunch_subset_count(bunch_slots, include=("IP1", "B1"))
            N_B1Only = pythonHelpers.bunch_struct.get_bunch_subset_count(bunch_slots, include=("B1",), exclude=("IP1", "IP2", "B2"))
            N_IP1_and_B2 = pythonHelpers.bunch_struct.get_bunch_subset_count(bunch_slots, include=("IP1", "B2"))
            N_B2noB1 = pythonHelpers.bunch_struct.get_bunch_subset_count(bunch_slots, include=("B2",), exclude=("IP1", "IP2", "B1"))
            N_IP1_and_IP2 = pythonHelpers.bunch_struct.get_bunch_subset_count(bunch_slots, include=("IP1", "IP2"))
            N_IP2Only = pythonHelpers.bunch_struct.get_bunch_subset_count(bunch_slots, include=("IP2",), exclude=("IP1", "B1", "B2"))


            for tt, total_counts in counts["IP1"].items():
                if N_B1Only>0:
                    output["B1Only"][tt] *= (N_IP1_and_B1 / N_B1Only)
                if N_B2noB1>0:
                    output["B2noB1"][tt] *= (N_IP1_and_B2 / N_B2noB1)
                if N_IP2Only > 0:
                    output["IP2"][tt] *= (N_IP1_and_IP2 / N_IP2Only)

        elif acc_mode==12: # Heavy-Ions
            N_IP1_and_IP2_and_B1_and_B2 = pythonHelpers.bunch_struct.get_bunch_subset_count(bunch_slots, include=("IP1", "IP2", "B1", "B2"), exclude=())
            N_IP2_and_B1_and_B2 = pythonHelpers.bunch_struct.get_bunch_subset_count(bunch_slots, include=("IP2", "B1", "B2"), exclude=("IP1",))


            for tt, total_counts in counts:
                output["IP2B1B2"] *= (N_IP1_and_IP2_and_B1_and_B2 / N_IP2_and_B1_and_B2)

    print(output)

#    pythonHelpers.nTracks.write_output(outfile_root, run, fill, output, args.scale)
    pythonHelpers.nTracks.write_output(outfile_csv,  run, fill, output, args.scale)



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
    parser.add_argument('--bunch-correction', action='store_true', help="Whether to apply bunch structure correction to the track count.")

    args = parser.parse_args()


    get_nTracks_pipeline(args)
