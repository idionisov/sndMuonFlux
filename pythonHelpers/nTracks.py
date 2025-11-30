import csv
import os
from array import array

import ROOT


def get_track_counts(
    std_vector: ROOT.std.vector(ROOT.std.array(float, 4)),
    option: str = "IP1"
) -> dict:
    options = {
        "IP1": 0, "IP2": 1, "B1Only": 2, "B2noB1": 3, "IP2B1B2": 4
    }
    if option not in options:
        raise ValueError(f"Invalid option '{option}'. Valid options are: {', '.join(options.keys())}")

    counts = {
        1:  std_vector.at(options[option]).at(0),
        11: std_vector.at(options[option]).at(1),
        3:  std_vector.at(options[option]).at(2),
        13: std_vector.at(options[option]).at(3)
    }
    return counts



def save_nTracks_to_csv(
    filename: str,
    run: int,
    fill: int,
    counts: dict,
    scale: int = 1,
    cols: list[str] = ["Run", "Fill", "scale", "nTracks1", "nTracks11", "nTracks3", "nTracks13"]
):
    if not filename.endswith(".csv"):
        filename += ".csv"


    row = {
        "Run": run,
        "Fill": fill,
        "scale": scale,
        "nTracks1": counts.get(1, 0),
        "nTracks11": counts.get(11, 0),
        "nTracks3": counts.get(3, 0),
        "nTracks13": counts.get(13, 0)
    }

    file_exists = os.path.exists(filename)
    with open(filename, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=cols)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)




def save_nTracks_to_root(
    filename: str,
    run: int,
    fill: int,
    counts: dict,
    scale: int = 1,
    cols: list[str] = ["Run", "Fill", "scale", "nTracks1", "nTracks11", "nTracks3", "nTracks13"],
    tree_name: str = "nTracks"
):

    if not filename.endswith(".root"):
        filename = filename + ".root"

    def _entry_to_row(tree):
        return {
            "Run": int(getattr(tree, cols[0])),
            "Fill": int(getattr(tree, cols[1])),
            "scale": int(getattr(tree, cols[2])),
            "nTracks1": int(getattr(tree, cols[3])),
            "nTracks11": int(getattr(tree, cols[4])),
            "nTracks3": int(getattr(tree, cols[5])),
            "nTracks13": int(getattr(tree, cols[6]))
        }

    existing_rows = []

    if os.path.exists(filename):
        f = ROOT.TFile(filename, "READ")
        existing_tree = f.Get(tree_name)
        if existing_tree:
            n_entries = int(existing_tree.GetEntries())
            for i in range(n_entries):
                existing_tree.GetEntry(i)
                existing_rows.append(_entry_to_row(existing_tree))
        f.Close()

    new_row = {
        "Run": int(run),
        "Fill": int(fill),
        "scale": int(scale),
        "nTracks1": int(counts.get(1, 0)),
        "nTracks11": int(counts.get(11, 0)),
        "nTracks3": int(counts.get(3, 0)),
        "nTracks13": int(counts.get(13, 0))
    }
    existing_rows.append(new_row)

    outfile = ROOT.TFile(filename, "UPDATE")
    old_tree = outfile.Get(tree_name)
    if old_tree:
        outfile.Delete(f"{tree_name};*")
    tree = ROOT.TTree(tree_name, "Number of tracks tree")

    run_buf = array('i', [0])
    fill_buf = array('i', [0])
    scale_buf = array('i', [0])
    n1_buf = array('i', [0])
    n11_buf = array('i', [0])
    n3_buf = array('i', [0])
    n13_buf = array('i', [0])

    tree.Branch(cols[0], run_buf, f"{cols[0]}/I")
    tree.Branch(cols[1], fill_buf, f"{cols[1]}/I")
    tree.Branch(cols[2], scale_buf, f"{cols[2]}/I")
    tree.Branch(cols[3], n1_buf, f"{cols[3]}/I")
    tree.Branch(cols[4], n11_buf, f"{cols[4]}/I")
    tree.Branch(cols[5], n3_buf, f"{cols[5]}/I")
    tree.Branch(cols[6], n13_buf, f"{cols[6]}/I")

    for row in existing_rows:
        run_buf[0] = int(row["Run"])
        fill_buf[0] = int(row["Fill"])
        scale_buf[0] = int(row["scale"])
        n1_buf[0] = int(row["nTracks1"])
        n11_buf[0] = int(row["nTracks11"])
        n3_buf[0] = int(row["nTracks3"])
        n13_buf[0] = int(row["nTracks13"])
        tree.Fill()

    tree.Write()
    outfile.Close()



def write_output(
    filename: str,
    run: int,
    fill: int,
    counts: dict,
    scale: int = 1,
    cols: list[str] = ["Run", "Fill", "scale", "nTracks1", "nTracks11", "nTracks3", "nTracks13"],
    tree_name: str = "nTracks"
):
    if not filename:
        return

    if filename.endswith(".csv"):
        save_nTracks_to_csv(filename, run, fill, counts, scale, cols)
    else:
        if not filename.endswith(".root"):
            filename += ".root"
        save_nTracks_to_root(filename, run, fill, counts, scale, cols, tree_name)

    print(f"Track counts have been saved to {filename}")
