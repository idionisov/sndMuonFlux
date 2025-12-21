import argparse
import os

import numpy as np
import ROOT

import pythonHelpers.general
import pythonHelpers.trkeff

this_dir = os.path.dirname(os.path.abspath(__file__))





def get_trkeff_pipeline(args):
    if args.mct:
        print("MC-Truth method is being used!")
    else:
        print("Tagging track method is being used!")

    # Assumes the analysis is compiled in build directory
    path_to_trkefflib = os.path.join(this_dir, ".", "build", "trkeff", "libtrkeffUtils.so")
    path_to_trkefflib = os.path.abspath(path_to_trkefflib)
    ROOT.gSystem.Load(path_to_trkefflib)

    if len(args.fout) == 0:
        print("No output file specified. Tracking efficiency results will be printed without saving.")
        outnames = args.fout
    else:
        outnames = " ".join(args.fout)
    outfile_root, outfile_csv = pythonHelpers.general.get_outfiles(outnames)

    ch, run, fill, acc_mode = pythonHelpers.general.load_run_info(args.input_files)
    if not pythonHelpers.general.is_mc(ch) or args.mct==False:
        vec = ROOT.computeTrackingEfficienciesPy(
            args.input_files,
            args.geofile,
            outfile_root,
            args.hist_params,
            *tuple(args.x_range), *tuple(args.y_range),
            *tuple(args.xz_range), *tuple(args.yz_range),
            *tuple(args.z_ref),
            args.veto_dist,
            args.us5_dist,
            args.sf_to_ds_dist,
            int(args.n_break)
        )
        effs = pythonHelpers.trkeff.get_effs_as_dict(vec)
    else:
        effs = pythonHelpers.trkeff.get_trkeff_mct(
            input_files = args.input_files,
            sigma = args.sigma,
            col_rate = args.col_rate,
            L_LHC = args.L_lhc,
            x_range = args.x_range,
            y_range = args.y_range,
            xz_range = args.xz_range,
            yz_range = args.yz_range,
            z_ref = args.z_ref
        )

    if outfile_csv:
        pythonHelpers.trkeff.save_trkeff_to_csv(outfile_csv, run, fill, effs)

    return effs



if __name__ == "__main__":
    default_hist_params = os.path.join(this_dir, "trkeff", "histParams.conf")


    parser = argparse.ArgumentParser(description="Script for computing the tracking efficiency.")

    parser.add_argument('-i', '--input-files', type=str, required=True, help="Regex pattern for input ROOT files with reconstructed tracks, e.g., '/path/to/files*.root'.")
    parser.add_argument('-g', '--geofile', type=str, required=True, help="Geofile to accurately compute distance of closest approach to MuFilter scintillator bars.")
    parser.add_argument('-o', '--fout', type=str, nargs="+", default=[], help="Optional output files (root by default to store histograms as well, but could be csv -- both formats simultaneously are supported). Csv files store only efficiencies while root files store root objets as well.")
    parser.add_argument('-z', '--z-ref', nargs="+", type=float, default=[430., 430., 450., 450.], help="Reference z-plane coordinates for each track type (types 1, 11, 3, 13). Provide 4 numbers: zRef1 zRef11 zRef3 zRef13.")
    parser.add_argument('-x', '--x-range', nargs=2, type=float, default=[-42., -10.], help="Fiducial x-coordinate range [xmin, xmax] in cm for tracks to be counted.")
    parser.add_argument('-y', '--y-range', nargs=2, type=float, default=[19., 48.], help="Fiducial y-coordinate range [ymin, ymax] in cm for tracks to be counted.")
    parser.add_argument('-xz', '--xz-range', nargs=2, type=float, default=[-1e12, 1e12], help="Allowed tagging track angle in XZ plane [xzMin, xzMax] in mrad. Tracks outside this range are ignored.")
    parser.add_argument('-yz', '--yz-range', nargs=2, type=float, default=[-1e12, 1e12], help="Allowed tagging track angle in YZ plane [yzMin, yzMax] in mrad. Tracks outside this range are ignored.")
    parser.add_argument('--veto-dist', type=float, default=3.0, help="Minimal distance to activated veto bar for DS tagging tracks.")
    parser.add_argument('--us5-dist', type=float, default=3.0, help="Minimal distance to activated us5 bar for SciFi tagging tracks.")
    parser.add_argument('--sf-to-ds-dist', type=float, default=3.0, help="Maximum distance between tagging and candidate tracks at reference plane for successful match.")
    parser.add_argument('--n-break', type=int, default=1e7, help="Breakpoint for the number of events processed.")
    parser.add_argument('--hist-params', type=str, default=default_hist_params, help="Histogram parameter config file.")
    parser.add_argument('-mct', '--mct', type=bool, default=False, help="Wether to use the Monte Carlo Truth method or not. Defaults to Tagging track method.")
    parser.add_argument('-x-sec', '--sigma', type=float, default=8e7, help="Cross section for inelastic hadron collisions for MC simulations.")
    parser.add_argument('-r', '--col-rate', type=float, default=100e6, help="Collision rate in Monte Carlo simulations.")
    parser.add_argument('--L-lhc', type=float, default=1, help="Luminosity to normalize to in nb.")

    args = parser.parse_args()


    effs = get_trkeff_pipeline(args)
    print(effs)
