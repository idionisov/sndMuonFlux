import argparse
import os
import time
import sys

# Add project root to path to allow importing from pythonHelpers
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import ROOT

import pythonHelpers.bunch_struct
import pythonHelpers.general
import pythonHelpers.lumi
import pythonHelpers.muon_flux
import pythonHelpers.nTracks
import pythonHelpers.trkeff

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

pythonHelpers.general.load_cpp_extension("libtrkeffUtils.so", "trkeff")
pythonHelpers.general.load_cpp_extension("libmuonFluxUtils.so", "muonFluxUtils")


def get_or_compute_track_counts(args, THIS_DIR: str):
    if any(n <= 0 for n in args.track_counts):
        print("\n" + "-"*50)
        print("Four valid track counts weren't provided. Starting track counter...")
        print("-"*50)

        pythonHelpers.general.load_cpp_extension("libnTracks.so", "nTracks")

        start_ts = args.t_range[0] if args.t_range[0] is not None else -1.0
        end_ts   = args.t_range[1] if args.t_range[1] is not None else -1.0

        std_vec = ROOT.getNTracks(
            args.input_files,
            *tuple(args.x_range), *tuple(args.y_range),
            *tuple(args.xz_range), *tuple(args.yz_range),
            *tuple(args.z_ref),
            start_ts,
            end_ts
        )

        return {
            "IP1":     pythonHelpers.nTracks.get_track_counts(std_vec, "IP1"),
            "IP2B1B2": pythonHelpers.nTracks.get_track_counts(std_vec, "IP2B1B2")
        }

    else:
        print("Using provided track counts")
        return {
            "IP1": {
                1:  args.track_counts[0],
                11: args.track_counts[1],
                3:  args.track_counts[2],
                13: args.track_counts[3]
            }
        }


def get_or_compute_efficiencies(args, outfile_root: str, THIS_DIR: str, is_mc: bool = False):
    if all(e > 0 for e in args.tracking_efficiencies) and \
       all(e > 0 for e in args.tracking_efficiency_errors):
        print("Using provided tracking efficiencies")
        effs = {
            1:  (args.tracking_efficiencies[0], args.tracking_efficiency_errors[0]),
            11: (args.tracking_efficiencies[1], args.tracking_efficiency_errors[1]),
            3:  (args.tracking_efficiencies[2], args.tracking_efficiency_errors[2]),
            13: (args.tracking_efficiencies[3], args.tracking_efficiency_errors[3])
        }
        return effs

    print("\n" + "-"*50)
    print("Evaluating tracking efficiencies...")
    print("-"*50)

    if not is_mc:
        effs_vec = ROOT.computeTrackingEfficiencies_TT(
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
        effs = pythonHelpers.trkeff.get_effs_as_dict(effs_vec)
    else:
        if args.py:
            effs = pythonHelpers.trkeff.get_trkeff_mct(
                input_files = args.input_files,
                sigma = args.sigma,
                col_rate = args.col_rate,
                L_LHC = args.L_lhc,
                x_range = args.x_range,
                y_range = args.y_range,
                z_ref = args.z_ref,
                xz_range = args.xz_range,
                yz_range = args.yz_range
            )
        else:
            effs_vec = ROOT.computeTrackingEfficiencies_MCT(
                args.input_files,
                args.sigma,
                args.col_rate,
                args.L_lhc,
                *tuple(args.x_range), *tuple(args.y_range),
                *tuple(args.z_ref),
                *tuple(args.xz_range), *tuple(args.yz_range)
            )
            effs = pythonHelpers.trkeff.get_effs_as_dict(effs_vec)

    return effs


def run_muon_flux_pipeline(args) -> dict:
    area = pythonHelpers.general.compute_area(args.x_range, args.y_range)

    # Output files
    if len(args.fout) == 0:
        print("No output file specified. Final results will be printed without saving.")
        outnames = args.fout
    else:
        outnames = " ".join(args.fout)
    outfile_root, outfile_csv = pythonHelpers.general.get_outfiles(outnames)

    # Load run, fill, accMode from TChain
    ch, run, fill, acc_mode = pythonHelpers.general.load_run_info(args.input_files)
    if ch.GetBranch("MCTrack"):
        is_mc = True
    else:
        is_mc = False


    if not is_mc:
        start_ts = args.t_range[0] if args.t_range else None
        end_ts   = args.t_range[1] if args.t_range else None
        lumi = pythonHelpers.lumi.get_lumi_eos(args.input_files, start_ts=start_ts, end_ts=end_ts)
        if acc_mode==12:
            lumi_err = 0.035*lumi
        elif acc_mode==11:
            lumi_err = 0.020*lumi
        else:
            lumi_err = 0


        # Track counts
        track_counts = get_or_compute_track_counts(args, THIS_DIR)

        # Tracking efficiencies
        effs = get_or_compute_efficiencies(args, outfile_root, THIS_DIR, is_mc)

        # Compute final muon flux
        muon_flux = pythonHelpers.muon_flux.compute_fluxes_per_track_type(
            track_counts, lumi, area, effs,
            lumi_err, args.scale,
            {1:1, 11:1, 3:1, 13:1} # No bunch correction
        )

        # Save outputs
        if outfile_root:
            pythonHelpers.trkeff.save_trkeff_to_root(outfile_root, run, fill, effs)
        
        pythonHelpers.muon_flux.write_output(
            outfile_root, run, fill, args.scale,
            track_counts, effs, muon_flux, "muonFlux"
        )
        pythonHelpers.muon_flux.write_output(
            outfile_csv, run, fill, args.scale,
            track_counts, effs, muon_flux, "muonFlux"
        )

        return {
            "run": run,
            "fill": fill,
            "acc_mode": acc_mode,
            "track_counts": track_counts,
            "effs": effs,
            "muon_flux": muon_flux
        }

    else:
        effs = get_or_compute_efficiencies(args, outfile_root, THIS_DIR, is_mc)

        muon_flux = pythonHelpers.muon_flux.compute_flux_mc_with_effs(
            input_files = args.input_files,
            sigma = args.sigma,
            col_rate = args.col_rate,
            L_LHC = args.L_lhc,
            x_range = args.x_range,
            y_range = args.y_range,
            z_ref = args.z_ref,
            xz_range = args.xz_range,
            yz_range = args.yz_range,
            effs = effs,
            python_eff = args.py
        )

        return {
            "run": np.nan,
            "fill": np.nan,
            "acc_mode": acc_mode,
            "effs": effs,
            "muon_flux": muon_flux
        }





if __name__ == "__main__":
    start = time.time()


    # Get project root (parent of scripts/)
    ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    default_hist_params = os.path.join(ROOT_DIR, "trkeff", "histParams.conf")

    parser = argparse.ArgumentParser(
        description="Script for computing the muon flux for an SND run.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )


    parser.add_argument('-i', '--input-files', type=str, required=True, help="Regex pattern for input ROOT files with reconstructed tracks, e.g., '/path/to/files*.root'.")
    parser.add_argument('-t', '--t-range', nargs=2, type=float, default=[None, None], help="(Optional) UTC timestamp range [start, end] in seconds for data selection.")
    parser.add_argument('-ntrks', '--track-counts', nargs=4, type=float, default=[0, 0, 0, 0], help="(Optional) Track counts for each track type (types 1, 11, 3, 13). Provide 4 numbers: nTrks1 nTrks11 nTrks3 nTrks13.")
    parser.add_argument('-effs', '--tracking-efficiencies', nargs=4, type=float, default=[0, 0, 0, 0], help="(Optional) Tracking efficiency for each track type (types 1, 11, 3, 13). Provide 4 numbers: eff1 eff11 eff3 eff13.")
    parser.add_argument('-effErr', '--tracking-efficiency-errors', nargs=4, type=float, default=[0, 0, 0, 0], help="(Optional) Tracking efficiency errors for each track type (types 1, 11, 3, 13). Provide 4 numbers: err1 err11 err3 err13.")
    parser.add_argument("-Lerr", "--lumi-relative-err", type=float, default=0, help="(Optional) Relative error on luminosity.")
    parser.add_argument('-g', '--geofile', type=str, default="", help="(Optional) Geofile to accurately compute distance of closest approach to MuFilter scintillator bars.")
    parser.add_argument('-o', '--fout', type=str, nargs="+", default=[], help="(Optional) Output files (root by default to store all objects, but could be csv -- both formats simultaneously are supported). Csv files store only efficiencies while root files store root objets as well.")
    parser.add_argument('-z', '--z-ref', nargs=4, type=float, default=[430., 430., 450., 450.], help="(Optional) Reference z-plane coordinates for each track type (types 1, 11, 3, 13). Provide 4 numbers: zRef1 zRef11 zRef3 zRef13.")
    parser.add_argument('-x', '--x-range', nargs=2, type=float, default=[-42., -10.], help="(Optional) Fiducial x-coordinate range [xmin, xmax] in cm for tracks to be counted.")
    parser.add_argument('-y', '--y-range', nargs=2, type=float, default=[19., 48.], help="(Optional) Fiducial y-coordinate range [ymin, ymax] in cm for tracks to be counted.")
    parser.add_argument('-xz', '--xz-range', nargs=2, type=float, default=[-1e12, 1e12], help="(Optional) Allowed tagging track angle in XZ plane [xzMin, xzMax] in mrad. Tracks outside this range are ignored.")
    parser.add_argument('-yz', '--yz-range', nargs=2, type=float, default=[-1e12, 1e12], help="(Optional) Allowed tagging track angle in YZ plane [yzMin, yzMax] in mrad. Tracks outside this range are ignored.")
    parser.add_argument('--tree-name', type=str, default="cbmsim", help="(Optional) Name of the TTree in input files")
    parser.add_argument('--veto-dist', type=float, default=3.0, help="(Optional) Minimal distance to activated veto bar for DS tagging tracks.")
    parser.add_argument('--us5-dist', type=float, default=3.0, help="(Optional) Minimal distance to activated us5 bar for SciFi tagging tracks.")
    parser.add_argument('--sf-to-ds-dist', type=float, default=3.0, help="(Optional) Maximum distance between tagging and candidate tracks at reference plane for successful match.")
    parser.add_argument('--n-break', type=int, default=1e7, help="(Optional) Breakpoint for the number of events processed (bunch structure and tracking efficiency).")
    parser.add_argument('-py', action='store_true', help="(Optional) Whether to use the python implementation, instead of the C++ one, of the tracking efficiency estimation.")
    parser.add_argument('--hist-params', type=str, default=default_hist_params, help="(Optional) Histogram parameter config file.")
    parser.add_argument('-sc', '--scale', type=int, default=1, help="(Optional) Scaling used for track reconstruction")
    parser.add_argument('-x-sec', '--sigma', type=float, default=8e7, help="(Optional) Cross section for inelastic hadron collisions for MC simulations.")
    parser.add_argument('-r', '--col-rate', type=float, default=100e6, help="(Optional) Collision rate in Monte Carlo simulations.")
    parser.add_argument('--L-lhc', type=float, default=1, help="(Optional) Luminosity to normalize to in nb.")

    args = parser.parse_args()
    results = run_muon_flux_pipeline(args)
    print(results['muon_flux'])




    elapsed = time.time() - start
    m, s = divmod(elapsed, 60)
    h, m = divmod(m, 60)
    print(f"Runtime: {int(h)}h {int(m)}m {int(s)}s  ({elapsed:.02f} seconds)")
