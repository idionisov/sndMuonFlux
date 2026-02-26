import csv
import os
from array import array

import numpy as np
import ROOT

import pythonHelpers.general
import pythonHelpers.trkeff


def compute_muon_flux(
    n_tracks: int, lumi: float, area: float,
    eff: float = 1.,
    eff_err: float = 0.,
    lumi_err: float = 0.,
    scale: int = 1,
    bunch_correction_factor: float = 1.,
):
    K = scale * bunch_correction_factor

    muon_flux = K * n_tracks / (lumi * area * eff)

    rel_N = 1.0 / np.sqrt(n_tracks) if n_tracks > 0 else 0.0
    rel_L = lumi_err / lumi
    rel_eff = eff_err / eff

    rel_total = np.sqrt(rel_N**2 + rel_L**2 + rel_eff**2)

    muon_flux_err = muon_flux * rel_total

    return muon_flux, muon_flux_err



def compute_fluxes_per_track_type(
    track_counts: dict,
    lumi: float,
    area: float,
    effs: dict,
    lumi_err: float,
    scale: int,
    correction: dict
) -> dict:
    muon_flux = {}
    for tt in track_counts["IP1"]:
        muon_flux[tt] = compute_muon_flux(
            track_counts["IP1"][tt],
            lumi, area,
            effs[tt][0], effs[tt][1],
            lumi_err, scale, correction[tt]
        )
    return muon_flux


def save_all_to_csv(
    filename: str, run: int, fill: int, scale: int,
    nTracks_dict: dict,
    effs_dict: dict,
    muon_flux_dict: dict,
    cols: list[str] = [
        "Run", "Fill", "scale",
        "nTracks1", "nTracks11", "nTracks3", "nTracks13",
        "trkeff1", "trkeff11", "trkeff3", "trkeff13",
        "trkeffErr1", "trkeffErr11", "trkeffErr3", "trkeffErr13",
        "flux1", "flux11", "flux3", "flux13",
        "fluxErr1", "fluxErr11", "fluxErr3", "fluxErr13"
    ]
):
    if not filename.endswith(".csv"):
        filename += ".csv"


    row = {
        "Run": run,
        "Fill": fill,
        "scale": scale,
        "nTracks1":      nTracks_dict.get(1,  0),
        "nTracks11":     nTracks_dict.get(11, 0),
        "nTracks3":      nTracks_dict.get(3,  0),
        "nTracks13":     nTracks_dict.get(13, 0),
        "trkeff1":       effs_dict.get(1,  0)[0],
        "trkeff11":      effs_dict.get(11, 0)[0],
        "trkeff3":       effs_dict.get(3,  0)[0],
        "trkeff13":      effs_dict.get(13, 0)[0],
        "trkeffErr1":    effs_dict.get(1,  0)[1],
        "trkeffErr11":   effs_dict.get(11, 0)[1],
        "trkeffErr3":    effs_dict.get(3,  0)[1],
        "trkeffErr13":   effs_dict.get(13, 0)[1],
        "flux1":         muon_flux_dict.get(1,  (0, 0))[0],
        "flux11":        muon_flux_dict.get(11, (0, 0))[0],
        "flux3":         muon_flux_dict.get(3,  (0, 0))[0],
        "flux13":        muon_flux_dict.get(13, (0, 0))[0],
        "fluxErr1":      muon_flux_dict.get(1,  (0, 0))[1],
        "fluxErr11":     muon_flux_dict.get(11, (0, 0))[1],
        "fluxErr3":      muon_flux_dict.get(3,  (0, 0))[1],
        "fluxErr13":     muon_flux_dict.get(13, (0, 0))[1]
    }

    file_exists = os.path.exists(filename)
    with open(filename, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=cols)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)



def save_muonFlux_to_root(
    filename: str,
    run: int,
    fill: int,
    muon_flux_dict: dict,
    cols: list[str] = ["Run", "Fill", "flux1", "flux11", "flux3", "flux13", "fluxErr1", "fluxErr11", "fluxErr3", "fluxErr13"],
    tree_name: str = "muonFlux"
):

    if not filename.endswith(".root"):
        filename = filename + ".root"

    def _entry_to_row(tree):
        return {
            "Run":       int(tree.GetLeaf(cols[0]).GetValue()),
            "Fill":      int(tree.GetLeaf(cols[1]).GetValue()),
            "flux1":     float(tree.GetLeaf(cols[2]).GetValue()),
            "flux11":    float(tree.GetLeaf(cols[3]).GetValue()),
            "flux3":     float(tree.GetLeaf(cols[4]).GetValue()),
            "flux13":    float(tree.GetLeaf(cols[5]).GetValue()),
            "fluxErr1":  float(tree.GetLeaf(cols[6]).GetValue()),
            "fluxErr11": float(tree.GetLeaf(cols[7]).GetValue()),
            "fluxErr3":  float(tree.GetLeaf(cols[8]).GetValue()),
            "fluxErr13": float(tree.GetLeaf(cols[9]).GetValue())
        }

    existing_rows = []
    if os.path.exists(filename):
        f = ROOT.TFile(filename, "READ")
        tree = f.Get(tree_name)
        if tree and isinstance(tree, ROOT.TTree):
            has_all_leaves = all(tree.GetLeaf(c) for c in cols)
            if has_all_leaves:
                for i in range(tree.GetEntries()):
                    tree.GetEntry(i)
                    existing_rows.append(_entry_to_row(tree))
        f.Close()

    new_row = {
        "Run":       int(run),
        "Fill":      int(fill),
        "flux1":     float(muon_flux_dict.get(1,  (0, 0))[0]),
        "flux11":    float(muon_flux_dict.get(11, (0, 0))[0]),
        "flux3":     float(muon_flux_dict.get(3,  (0, 0))[0]),
        "flux13":    float(muon_flux_dict.get(13, (0, 0))[0]),
        "fluxErr1":  float(muon_flux_dict.get(1,  (0, 0))[1]),
        "fluxErr11": float(muon_flux_dict.get(11, (0, 0))[1]),
        "fluxErr3":  float(muon_flux_dict.get(3,  (0, 0))[1]),
        "fluxErr13": float(muon_flux_dict.get(13, (0, 0))[1])
    }
    
    is_duplicate = any(r["Run"] == run for r in existing_rows)
    if not is_duplicate:
        existing_rows.append(new_row)

    outfile = ROOT.TFile(filename, "UPDATE")
    old_tree = outfile.Get(tree_name)
    if old_tree:
        outfile.Delete(f"{tree_name};*")
    tree = ROOT.TTree(tree_name, "Muon flux tree")
    ROOT.SetOwnership(tree, False)

    run_buf     = array('i', [0])
    fill_buf    = array('i', [0])
    mf1_buf     = array('f', [0])
    mf11_buf    = array('f', [0])
    mf3_buf     = array('f', [0])
    mf13_buf    = array('f', [0])
    mferr1_buf  = array('f', [0])
    mferr11_buf = array('f', [0])
    mferr3_buf  = array('f', [0])
    mferr13_buf = array('f', [0])

    tree.Branch(cols[0], run_buf,  f"{cols[0]}/I")
    tree.Branch(cols[1], fill_buf, f"{cols[1]}/I")
    tree.Branch(cols[2], mf1_buf,  f"{cols[2]}/F")
    tree.Branch(cols[3], mf11_buf, f"{cols[3]}/F")
    tree.Branch(cols[4], mf3_buf,  f"{cols[4]}/F")
    tree.Branch(cols[5], mf13_buf, f"{cols[5]}/F")
    tree.Branch(cols[6], mferr1_buf,  f"{cols[6]}/F")
    tree.Branch(cols[7], mferr11_buf, f"{cols[7]}/F")
    tree.Branch(cols[8], mferr3_buf,  f"{cols[8]}/F")
    tree.Branch(cols[9], mferr13_buf, f"{cols[9]}/F")

    for row in existing_rows:
        run_buf[0]     = int(row["Run"])
        fill_buf[0]    = int(row["Fill"])
        mf1_buf[0]     = float(row["flux1"])
        mf11_buf[0]    = float(row["flux11"])
        mf3_buf[0]     = float(row["flux3"])
        mf13_buf[0]    = float(row["flux13"])
        mferr1_buf[0]  = float(row["fluxErr1"])
        mferr11_buf[0] = float(row["fluxErr11"])
        mferr3_buf[0]  = float(row["fluxErr3"])
        mferr13_buf[0] = float(row["fluxErr13"])
        tree.Fill()

    tree.Write()
    outfile.Close()

def write_output(
    filename: str,
    run: int, fill: int, scale: int,
    track_counts: dict,
    effs: dict,
    muon_flux: dict,
    tree_name: str = ""
):
    if not filename:
        return

    if filename.endswith(".csv"):
        save_all_to_csv(
            filename, run, fill, scale,
            track_counts["IP1"], effs, muon_flux
        )

    else:
        if not filename.endswith(".root"):
            filename += ".root"

        save_muonFlux_to_root(
            filename, run, fill, muon_flux, tree_name=tree_name
        )



def compute_muon_flux_mc_asymm(Nr, dNr, L_LHC, A, eps, eps_err_up, eps_err_down):
    phi = Nr / (A * eps * L_LHC)

    phi_upper = (Nr + dNr) / (A * L_LHC * (eps - eps_err_down))
    phi_lower = (Nr - dNr) / (A * L_LHC * (eps + eps_err_up))

    sigma_phi_up = phi_upper - phi
    sigma_phi_down = phi - phi_lower

    phiStatVar = ( dNr / (A * eps * L_LHC) )**2
    phiEffVar = (((eps_err_up + eps_err_down)/2) * ( Nr  /  (A * L_LHC * eps**2) ))**2

    return {
        "results": (phi, sigma_phi_up, sigma_phi_down),
        "variances": {
            "stat": phiStatVar, "eff": phiEffVar
        }
    }


def compute_muon_flux_mc(Nr, dNr, L, A, eps, eps_err, Ks):
    phi = (Nr * Ks) / (A * eps * L)

    phiVarStat = ((dNr * Ks) / (A * eps * L))**2
    phiVarEff  = ((Ks * Nr * eps_err) / (A * eps**2 * L))**2

    phiErr = np.sqrt(phiVarStat + phiVarEff)

    return {
        "results": (phi, phiErr),
        "variances": {
            "stat": phiVarStat, "eff": phiVarEff
        }
    }

def compute_flux_mc_with_effs(
    input_files: str,
    sigma:       float        = 8e7,    # [nb]
    col_rate:    float        = 100e6,  # [s^-1]
    L_LHC:       float        = 1,      # [nb]
    x_range:     tuple[float, float] = (-42., -10.), # [cm]
    y_range:     tuple[float, float] = ( 19.,  48.), # [cm]
    z_ref:       tuple[float, float, float, float] = (430., 430., 450., 450.), # [cm]
    xz_range:    tuple[float, float] = (-1e12, 1e12), # [mrad]
    yz_range:    tuple[float, float] = (-1e12, 1e12), # [mrad]
    effs:        dict = {},
    python_eff:  bool = False
):
    L_MC = col_rate/sigma       # [nb^-1 s^-1]
    Ks = L_LHC / L_MC           # []

    trackTypes = (1, 11, 3, 13, "mcTrk")
    xy = {
        'min': {'x': x_range[0], 'y': y_range[0]},
        'max': {'x': x_range[1], 'y': y_range[1]}
    }

    z_ref_dict = {1: z_ref[0], 11: z_ref[1], 3: z_ref[2], 13: z_ref[3]}
    A = (xy['max']['x']-xy['min']['x']) * (xy['max']['y']-xy['min']['y'])

    ch = pythonHelpers.general.load_snd_TChain(input_files)
    n_entries = ch.GetEntries()

    eps = effs.copy()
    eps["mcTrk"] = (1.0, 0.0)

    Nrate               = {tt: 0.0 for tt in trackTypes}
    Nrate_err2          = {tt: 0.0 for tt in trackTypes}

    next_print = 0
    if python_eff:
        print("Using python implementation...")
        for i_entry, entry in enumerate(ch):
            progress = (i_entry * 100) // n_entries
            if progress >= next_print:
                print(f"{progress:.0f} %")
                next_print += 5

            if not entry.EventHeader.isIP1():
                continue

            w = 0.0
            for mctrack in entry.MCTrack:
                if (
                    mctrack.GetMotherId()==-1 and
                    abs(mctrack.GetPdgCode())==13
                ):
                    w = mctrack.GetWeight()

                    mcTrkZref = pythonHelpers.general.get_point_at_z(mctrack, 430)
                    if (
                        xy['min']['x'] <= mcTrkZref.X() <= xy['max']['x'] and
                        xy['min']['y'] <= mcTrkZref.Y() <= xy['max']['y']
                    ):
                        Nrate["mcTrk"] += w
                        Nrate_err2["mcTrk"] += w*w
                        break

            reco_track_types = set()
            for trk in entry.Reco_MuonTracks:
                if not (trk.getTrackFlag() and trk.getTrackMom().Z()):
                    continue

                tt = trk.getTrackType()
                ref = pythonHelpers.general.get_point_at_z(trk, z_ref_dict[tt])
                x = ref.X()
                y = ref.Y()

                if not (
                    xy["min"]["x"] <= x <= xy["max"]["x"] and
                    xy["min"]["y"] <= y <= xy["max"]["y"]
                ):
                    continue
                reco_track_types.add(tt)

            for tt in reco_track_types:
                Nrate[tt] += w
                Nrate_err2[tt] += w*w
    else:
        print("Using C++ implementation...")
        this_dir_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        lib_path = os.path.abspath(os.path.join(this_dir_path, "build", "muonFluxUtils", "libmuonFluxUtils.so"))
        ROOT.gSystem.Load(lib_path)

        cpp_res = ROOT.computeMCRates(
            input_files,
            x_range[0], x_range[1],
            y_range[0], y_range[1],
            z_ref_dict[1], z_ref_dict[11], z_ref_dict[3], z_ref_dict[13]
        )

        Nrate["mcTrk"] = cpp_res.nRate[0]
        Nrate_err2["mcTrk"] = cpp_res.nRateErr2[0]
        for tt in (1, 11, 3, 13):
            Nrate[tt] = cpp_res.nRate[tt]
            Nrate_err2[tt] = cpp_res.nRateErr2[tt]

    mf = {}
    for tt in trackTypes:
        Nr = Nrate[tt]
        dNr = np.sqrt(Nrate_err2[tt])

        D_sym  = compute_muon_flux_mc(
            Nr, dNr, L_LHC, A, eps[tt][0], eps[tt][1], Ks
        )
        phi, phiErr = D_sym["results"]
        mf[tt] = (phi, phiErr)

        print(f"\n > {tt}:\tΦ = {phi:.03f} ± {phiErr:.03f}    \033[1;39m[fb/cm²]\033[0m")

    return mf


def get_muon_flux_mc(
    input_files: str,
    sigma:       float        = 8e7,    # [nb]
    col_rate:    float        = 100e6,  # [s^-1]
    L_LHC:       float        = 1,      # [nb]
    x_range:     tuple[float, float] = (-42., -10.), # [cm]
    y_range:     tuple[float, float] = ( 19.,  48.), # [cm]
    z_ref:       tuple[float, float, float, float] = (430., 430., 450., 450.), # [cm]
    xz_range:    tuple[float, float] = (-1e12, 1e12), # [mrad]
    yz_range:    tuple[float, float] = (-1e12, 1e12), # [mrad]
    trk_eff:      tuple[float, float, float, float] = (0, 0, 0, 0),
    trkeff_err:  tuple[float, float, float, float] = (0, 0, 0, 0),
    python_eff:  bool = False
):
    L_MC = col_rate/sigma       # [nb^-1 s^-1]
    Ks = L_LHC / L_MC           # []

    trackTypes = (1, 11, 3, 13, "mcTrk")
    xy = {
        'min': {'x': x_range[0], 'y': y_range[0]},
        'max': {'x': x_range[1], 'y': y_range[1]}
    }

    z_ref = {1: z_ref[0], 11: z_ref[1], 3: z_ref[2], 13: z_ref[3]}
    A = (xy['max']['x']-xy['min']['x']) * (xy['max']['y']-xy['min']['y'])

    ch = pythonHelpers.general.load_snd_TChain(input_files)
    n_entries = ch.GetEntries()

    print(f"{ 'Entries':<10} : {n_entries:,}")
    print(f"{ 'Luminosity':<10} : {L_MC}")
    print(f"{ 'Area':<10} : {A}")
    print(f"{ 'Scaling':<10} : {Ks}")

    if (
        any(eff == 0 for eff in trk_eff) or
        any(err == 0 for err in trkeff_err)
    ):
        print("Tracking efficiencies were not provided and will be evaluated...")
        if python_eff:
            print("\tUsing the python implementation")
            eps = pythonHelpers.trkeff.get_trkeff_mct(
                input_files = input_files,
                sigma = sigma,
                col_rate = col_rate,
                L_LHC = L_LHC,
                x_range = x_range,
                y_range = y_range,
                z_ref = (z_ref[1], z_ref[11], z_ref[3], z_ref[13]),
                xz_range = xz_range,
                yz_range = yz_range,
            )
        else:
            print("\tUsing the C++ implementation")
            eps = ROOT.computeTrackingEfficiencies_MCT(
                input_files,
                sigma,
                col_rate,
                L_lhc,
                *tuple(x_range), *tuple(y_range),
                z_ref[1], z_ref[11], z_ref[3], z_ref[13],
                *tuple(xz_range), *tuple(yz_range)
            )
            eps = pythonHelpers.trkeff.get_effs_as_dict(eps)

    else:
        print("Adopting the provided efficiencies")
        eps = {
            1:  (trk_eff[0], trkeff_err[0]),
            11: (trk_eff[1], trkeff_err[1]),
            3:  (trk_eff[2], trkeff_err[2]),
            13: (trk_eff[3], trkeff_err[3])
        }
        print(eps)
    eps["mcTrk"] = (1, 0)

    eventMCmu = {}

    print("Starting evaluation of the track rate...")
    Nrate               = {tt: 0.0 for tt in trackTypes}
    Nrate_err2          = {tt: 0.0 for tt in trackTypes}
    Nrate["mcTrk"]      = 0.0
    Nrate_err2["mcTrk"] = 0.0

    next_print = 0
    if python_eff:
        print("\tStarting python implementation of track rate evaluation...")
        for i_entry, entry in enumerate(ch):
            progress = (i_entry * 100) // n_entries
            if progress >= next_print:
                print(f"{progress:.0f} %")
                next_print += 5

            if not entry.EventHeader.isIP1():
                continue

            for mctrack in entry.MCTrack:
                if (
                    mctrack.GetMotherId()==-1 and
                    abs(mctrack.GetPdgCode())==13
                ):
                    eventMCmu[i_entry] = mctrack.GetWeight()
                    w = eventMCmu.get(i_entry, 0.0)

                    mcTrkZref = pythonHelpers.general.get_point_at_z(mctrack, 430)
                    if (
                        xy['min']['x'] <= mcTrkZref.X() <= xy['max']['x'] and
                        xy['min']['y'] <= mcTrkZref.Y() <= xy['max']['y']
                    ):
                        Nrate["mcTrk"] += w
                        Nrate_err2["mcTrk"] += w*w
                        break

            reco_track_types = set()
            for trk in entry.Reco_MuonTracks:
                if not (trk.getTrackFlag() and trk.getTrackMom().Z()):
                    continue

                tt = trk.getTrackType()
                ref = pythonHelpers.general.get_point_at_z(trk, z_ref[tt])
                x = ref.X()
                y = ref.Y()

                if not (
                    xy["min"]["x"] <= x <= xy["max"]["x"] and
                    xy["min"]["y"] <= y <= xy["max"]["y"]
                ):
                    continue
                reco_track_types.add(tt)

            for tt in reco_track_types:
                Nrate[tt] += w
                Nrate_err2[tt] += w*w
    else:
        print("\tStarting C++ implementation of track rate evaluation...")
        this_dir_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        lib_path = os.path.abspath(os.path.join(this_dir_path, "build", "muonFluxUtils", "libmuonFluxUtils.so"))
        ROOT.gSystem.Load(lib_path)

        cpp_res = ROOT.computeMCRates(
            input_files,
            x_range[0], x_range[1],
            y_range[0], y_range[1],
            z_ref[1], z_ref[11], z_ref[3], z_ref[13]
        )

        Nrate["mcTrk"] = cpp_res.nRate[0]
        Nrate_err2["mcTrk"] = cpp_res.nRateErr2[0]
        for tt in (1, 11, 3, 13):
            Nrate[tt] = cpp_res.nRate[tt]
            Nrate_err2[tt] = cpp_res.nRateErr2[tt]

    mf, statVars, effVars = {}, {}, {}

    mf = {}
    statVars = {}
    effVars = {}

    for tt in trackTypes:
        Nr = Nrate[tt]
        dNr = np.sqrt(Nrate_err2[tt])


        D_sym  = compute_muon_flux_mc(
            Nr, dNr, L_LHC, A, eps[tt][0], eps[tt][1], Ks
        )
        phi, phiErr = D_sym["results"]
        statVars[tt] = D_sym["variances"]["stat"]
        effVars[tt] = D_sym["variances"]["eff"]

        mf[tt]     = (phi, phiErr)
        varPhiStat = statVars[tt]
        varPhiEff  = effVars[tt]

        # print(f"\n > {tt}:\tΦ = {phi:.03f} ± {phiErr:.03f}    \033[1;39m[fb/cm²]\033[0m")
        # # print(f"\n > {tt}:\tNr/A = {Nr/A}    \033[1;39m[cm^-2 s^-1]\033[0m")
        # # print(f"\n > {tt}:\tNr/(A.L_MC) = {Nr/(A*L_MC)}    \033[1;39m [fb cm^-2]\033[0m")
        # print(f"        \tError:  {phiErr*100/phi:.03f} %")
        # print(f"        \tStatistics: {varPhiStat} ({varPhiStat*100/(varPhiStat+varPhiEff):.02f} %)")
        # print(f"        \tEfficiency: {varPhiEff} ({varPhiEff*100/(varPhiStat+varPhiEff):.02f} %)")


    return mf
