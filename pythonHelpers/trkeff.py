import csv
import os

import numpy as np
import ROOT
from scipy.stats import beta

import pythonHelpers.general


def get_effs_as_dict(
    effs_std_vector: ROOT.std.vector(float)
) -> dict:
    return { # (eff, effErr)
        1:  (effs_std_vector.at(0), effs_std_vector.at(1)),
        11: (effs_std_vector.at(2), effs_std_vector.at(3)),
        3:  (effs_std_vector.at(4), effs_std_vector.at(5)),
        13: (effs_std_vector.at(6), effs_std_vector.at(7))
    }




def save_trkeff_to_csv(
    filename: str,
    run: int,
    fill: int,
    effs: dict,
    cols: list[str] = [
        "Run", "Fill", "trkeff1", "trkeff11", "trkeff3", "trkeff13",
        "trkeffErr1", "trkeffErr11", "trkeffErr3", "trkeffErr13"
    ]
):
    if not filename.endswith(".csv"):
        filename += ".csv"

    row = {
        "Run": run,
        "Fill": fill,
        "trkeff1":      effs.get(1,  (0,0))[0],
        "trkeff11":     effs.get(11, (0,0))[0],
        "trkeff3":      effs.get(3,  (0,0))[0],
        "trkeff13":     effs.get(13, (0,0))[0],
        "trkeffErr1":   effs.get(1,  (0,0))[1],
        "trkeffErr11":  effs.get(11, (0,0))[1],
        "trkeffErr3":   effs.get(3,  (0,0))[1],
        "trkeffErr13":  effs.get(13, (0,0))[1]
    }

    file_exists = os.path.exists(filename)
    with open(filename, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=cols)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)


def save_trkeff_to_root(
    filename: str,
    run: int,
    fill: int,
    effs: dict,
    tree_name: str = "trkeff"
):
    from array import array
    if not filename.endswith(".root"):
        filename += ".root"

    # Define columns to match CSV
    cols = [
        "Run", "Fill", "trkeff1", "trkeff11", "trkeff3", "trkeff13",
        "trkeffErr1", "trkeffErr11", "trkeffErr3", "trkeffErr13"
    ]

    def _entry_to_row(tree):
        return {
            "Run":        int(tree.GetLeaf(cols[0]).GetValue()),
            "Fill":       int(tree.GetLeaf(cols[1]).GetValue()),
            "trkeff1":    float(tree.GetLeaf(cols[2]).GetValue()),
            "trkeff11":   float(tree.GetLeaf(cols[3]).GetValue()),
            "trkeff3":    float(tree.GetLeaf(cols[4]).GetValue()),
            "trkeff13":   float(tree.GetLeaf(cols[5]).GetValue()),
            "trkeffErr1": float(tree.GetLeaf(cols[6]).GetValue()),
            "trkeffErr11":float(tree.GetLeaf(cols[7]).GetValue()),
            "trkeffErr3": float(tree.GetLeaf(cols[8]).GetValue()),
            "trkeffErr13":float(tree.GetLeaf(cols[9]).GetValue())
        }

    existing_rows = []
    if os.path.exists(filename):
        f = ROOT.TFile(filename, "READ")
        tree = f.Get(tree_name)
        if tree and isinstance(tree, ROOT.TTree):
            # Check if it has the leaves we expect for our standardization
            has_all_leaves = all(tree.GetLeaf(c) for c in cols)
            if has_all_leaves:
                for i in range(tree.GetEntries()):
                    tree.GetEntry(i)
                    existing_rows.append(_entry_to_row(tree))
            else:
                print(f"INFO: Tree '{tree_name}' exists but structure differs. Standardizing...")
        f.Close()

    # Add new data
    new_row = {
        "Run": run, "Fill": fill,
        "trkeff1":      effs.get(1,  (0,0))[0],
        "trkeff11":     effs.get(11, (0,0))[0],
        "trkeff3":      effs.get(3,  (0,0))[0],
        "trkeff13":     effs.get(13, (0,0))[0],
        "trkeffErr1":   effs.get(1,  (0,0))[1],
        "trkeffErr11":  effs.get(11, (0,0))[1],
        "trkeffErr3":   effs.get(3,  (0,0))[1],
        "trkeffErr13":  effs.get(13, (0,0))[1]
    }

    # Avoid duplicates if we are standardizing an existing file
    is_duplicate = any(r["Run"] == run for r in existing_rows)
    if not is_duplicate:
        existing_rows.append(new_row)

    f = ROOT.TFile(filename, "UPDATE")
    # If the tree already exists, we might want to overwrite it with the standardized version
    old_tree = f.Get(tree_name)
    if old_tree:
        f.Delete(f"{tree_name};*")

    tree = ROOT.TTree(tree_name, "Tracking efficiencies")

    # Buffers
    b_run = array('i', [0]); tree.Branch("Run", b_run, "Run/I")
    b_fill = array('i', [0]); tree.Branch("Fill", b_fill, "Fill/I")
    b_eff1 = array('f', [0]); tree.Branch("trkeff1", b_eff1, "trkeff1/F")
    b_eff11 = array('f', [0]); tree.Branch("trkeff11", b_eff11, "trkeff11/F")
    b_eff3 = array('f', [0]); tree.Branch("trkeff3", b_eff3, "trkeff3/F")
    b_eff13 = array('f', [0]); tree.Branch("trkeff13", b_eff13, "trkeff13/F")
    b_err1 = array('f', [0]); tree.Branch("trkeffErr1", b_err1, "trkeffErr1/F")
    b_err11 = array('f', [0]); tree.Branch("trkeffErr11", b_err11, "trkeffErr11/F")
    b_err3 = array('f', [0]); tree.Branch("trkeffErr3", b_err3, "trkeffErr3/F")
    b_err13 = array('f', [0]); tree.Branch("trkeffErr13", b_err13, "trkeffErr13/F")

    for row in existing_rows:
        b_run[0] = row["Run"]
        b_fill[0] = row["Fill"]
        b_eff1[0] = row["trkeff1"]
        b_eff11[0] = row["trkeff11"]
        b_eff3[0] = row["trkeff3"]
        b_eff13[0] = row["trkeff13"]
        b_err1[0] = row["trkeffErr1"]
        b_err11[0] = row["trkeffErr11"]
        b_err3[0] = row["trkeffErr3"]
        b_err13[0] = row["trkeffErr13"]
        tree.Fill()

    tree.Write()
    f.Close()



def there_is_a_muon(
    mc_event: ROOT.TChain
) -> bool:
    for mcTrack in mc_event.MCTrack:
        if abs(mcTrack.GetPdgCode()) == 13:
            return True
    return False


def mc_track_crossed_fiducial_area(
    entry: ROOT.TChain,
    z_ref: float = 430.,
    x_range: tuple = (-42., -10.),
    y_range: tuple = (18., 49.)
) -> bool:
    for mctrack in entry.MCTrack:
        if (
            mctrack.GetMotherId()==-1 and
            abs(mctrack.GetPdgCode()) == 13
        ):
            mcTrkZref = pythonHelpers.general.get_point_at_z(mctrack, z_ref)
            if (
                x_range[0] <= mcTrkZref.X() <= x_range[1] and
                y_range[0] <= mcTrkZref.Y() <= y_range[1]
            ):
                return True
    return False


def sf_track_is_reconstructible(
    mcEvent: ROOT.TChain
) -> dict:
    nMCPoints = {
        'h': {1:0, 2:0, 3:0, 4:0, 5:0},
        'v': {1:0, 2:0, 3:0, 4:0, 5:0}
    }

    nScifiPoints = {'h': 0, 'v': 0}

    for mcPoint in mcEvent.ScifiPoint:
        if not (
            abs(mcPoint.PdgCode()) == 13 and
            mcPoint.GetTrackID() == 0
        ): continue

        detID = mcPoint.GetDetectorID()

        # Checking for reconstructibility of SciFi tracks

        # Second digit:
        #   - type of the plane:
        #       - 0: Horizontal fiber plane
        #       - 1: Vertical fiber plane
        if int( (detID/100000)%2 ) == 0:
            nMCPoints['h'][int(detID/1e+6)]+=1
        elif int( (detID/100000)%2 ) == 1:
            nMCPoints['v'][int(detID/1e+6)]+=1

    for sfPlane in range(1, len(nMCPoints['v'])+1):
        if nMCPoints['v'][sfPlane]>0:
            nScifiPoints['v']+=1

        if nMCPoints['h'][sfPlane]>0:
            nScifiPoints['h']+=1

    if (nScifiPoints['h'] >= 3 and nScifiPoints['v'] >= 3):
        return True
    else:
        return False



def ds_track_is_reconstructible(
    mcEvent: ROOT.TChain
) -> dict:
    nMCPoints = {
        'h': {1:0, 2:0, 3:0, 4:0},
        'v': {1:0, 2:0, 3:0, 4:0}
    }
    nDSPoints = {'h': 0, 'v': 0}

    for mcPoint in mcEvent.MuFilterPoint:
        if (
            abs(mcPoint.PdgCode())==13 and
            mcPoint.GetTrackID()==0
        ):
            detID = mcPoint.GetDetectorID()
            if detID < 30000 or detID > 34999:
                continue

            station = int(str(detID)[1])+1
            barNum = int(str(detID)[-3:])

            if barNum > 59:
                nMCPoints['v'][station] += 1
            else:
                nMCPoints['h'][station] += 1


    for dsPlane in range(1, 5):
        if nMCPoints['v'][dsPlane] > 0: nDSPoints['v'] += 1
        if nMCPoints['h'][dsPlane] > 0: nDSPoints['h'] += 1

    if (nDSPoints['h'] >= 3 and nDSPoints['v'] >= 3):
        return True
    else:
        return False

def n_eff(weights):
    w = np.asarray(weights, dtype=float)
    s1 = w.sum()
    s2 = (w**2).sum()
    if s2 == 0:
        return 0.0
    return (s1*s1) / s2




def get_trkeff_mct(
    input_files: str,
    sigma:       float        = 8e7,
    col_rate:    float        = 100e6,
    L_LHC:       float        = 1,
    x_range:     tuple[float] = (-42., -10.),
    y_range:     tuple[float] = ( 19.,  48.),
    z_ref:       tuple[float] = (430., 430., 450., 450.),
    xz_range:    tuple[float] = (-1e12, 1e12),
    yz_range:    tuple[float] = (-1e12, 1e12)
):
    L_MC = col_rate / sigma
    Ks = L_LHC / L_MC
    trackTypes = (1, 11, 3, 13)

    xy = {
        'min': {'x': x_range[0], 'y': y_range[0]},
        'max': {'x': x_range[1], 'y': y_range[1]}
    }
    A = (xy['max']['x'] - xy['min']['x']) * (xy['max']['y'] - xy['min']['y'])
    xz_min, xz_max = xz_range
    yz_min, yz_max = yz_range

    ch = pythonHelpers.general.load_snd_TChain(input_files)
    n_entries = ch.GetEntries()
    print(f"{'Entries':<10} : {n_entries:,}")
    print(f"{'Luminosity':<10} : {L_MC}")
    print(f"{'Area':<10} : {A}")
    print(f"{'Scaling':<10} : {Ks}")

    weights = {tt: [] for tt in trackTypes}
    passed = {tt: [] for tt in trackTypes}

    next_print = 0
    for i_entry, entry in enumerate(ch):
        progress = (i_entry * 100) // n_entries
        if progress >= next_print:
            print(f"{progress:.0f} %")
            next_print += 5

        if not (
            entry.EventHeader.isIP1() and
            pythonHelpers.trkeff.there_is_a_muon(entry)
        ):
            continue

        _sf = pythonHelpers.trkeff.sf_track_is_reconstructible(entry)
        _ds = pythonHelpers.trkeff.ds_track_is_reconstructible(entry)
        muon_passed_within_A = {
            tt: mc_track_crossed_fiducial_area(entry,
                z_ref = z_ref[i_tt],
                x_range = x_range,
                y_range = y_range
            ) for i_tt, tt in enumerate((1, 11, 3, 13))
        }

        if _sf==False and _ds==False:
            continue

        reco = {
            1:  _sf and muon_passed_within_A[1],
            11: _sf and muon_passed_within_A[11],
            3:  _ds and muon_passed_within_A[3],
            13: _ds and muon_passed_within_A[13],
        }
        weight = entry.MCTrack[0].GetWeight() * Ks

        for tt in trackTypes:
            if not reco[tt]:
                continue

            eventPassed = False
            weights[tt].append(weight)
            for trk in entry.Reco_MuonTracks:
                if trk.getTrackType() != tt:
                    continue

                xz = trk.getAngleXZ()
                yz = trk.getAngleYZ()
                if not (
                    xz_min <= xz <= xz_max and
                    yz_min <= yz <= yz_max
                ):
                    continue

                eventPassed = True
                break
            passed[tt].append(eventPassed)


    for tt in trackTypes:
        weights[tt] = np.array(weights[tt], dtype=np.float64)
        passed[tt]  = np.array(passed[tt], dtype=bool)

    eff, effErrUp, effErrLow = {}, {}, {}
    for tt in (1, 11, 3, 13):

        w_pass = weights[tt][passed[tt]].sum()
        w_total = weights[tt].sum()

        # w2_pass = (weights[tt][passed[tt]] ** 2).sum()
        # w2_total = (weights[tt] ** 2).sum()

        p_hat = w_pass / w_total  if w_total>0 else 0.0
        n_eff = pythonHelpers.trkeff.n_eff(weights[tt])
        k_eff = p_hat * n_eff

        a = k_eff + 1.0
        b = n_eff - k_eff + 1.0
        low_q, high_q = 0.15865, 0.84135
        lo_b = beta.ppf(low_q, a, b)
        hi_b = beta.ppf(high_q, a, b)

        eff[tt], effErrUp[tt], effErrLow[tt] = p_hat, hi_b - p_hat, p_hat - lo_b

    return {tt: (eff[tt], (effErrUp[tt]+effErrLow[tt])/2) for tt in (1, 11, 3, 13)}
