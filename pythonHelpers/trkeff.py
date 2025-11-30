import csv
import os

import ROOT


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
        "trkeff1":      effs.get(1,  0)[0],
        "trkeff11":     effs.get(11, 0)[0],
        "trkeff3":      effs.get(3,  0)[0],
        "trkeff13":     effs.get(13, 0)[0],
        "trkeffErr1":   effs.get(1,  0)[1],
        "trkeffErr11":  effs.get(11, 0)[1],
        "trkeffErr3":   effs.get(3,  0)[1],
        "trkeffErr13":  effs.get(13, 0)[1]
    }

    file_exists = os.path.exists(filename)
    with open(filename, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=cols)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)
