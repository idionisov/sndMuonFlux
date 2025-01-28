import ROOT
import pandas as pd
import numpy as np
import uproot
from sndUtils import sys, alg, system, algorithm, att
import matplotlib.ticker as mticker
from ddfRoot import getPandasFromUprootTH1
import matplotlib.pyplot as plt
plt.style.use("root")

mfout = "/eos/user/i/idioniso/mfout"
chi2ndf = {
    1:  [5,   10,  15,  20, 40,  60, 100],
    11: [8,   12,  16,  25, 50,  80, 120],
    3:  [2,   3,   4,   5,  7.5, 10, 15 ],
    13: [3.5, 4.5, 5.5, 7,  9,   12, 15 ]
}
dfEffChi2ndf = pd.read_csv(f"{mfout}/trkeff/trkeff-chi2ndfMax/trkeff-chi2ndf.csv")

chi2ndf = {
    tt: getPandasFromUprootTH1(
        uproot.open(f"/eos/user/i/idioniso/mfout/MuonFlux.root:TrackProperties/Run7080/h_chi2ndf_{tt}_7080")
    ) for tt in (1, 11, 3, 13)
}


fig, axs = plt.subplots(2, 2, figsize=(17,8))

for i_tt, tt in enumerate((1, 11, 3, 13)):
    _att = att(tt)
    ax_i, ax_j = divmod(i_tt, 2)

    l11 = axs[ax_i, ax_j].errorbar(
        x=dfEffChi2ndf[f'chi2ndf_{tt}'],
        y=dfEffChi2ndf[f'trkeff_{tt}'],
        yerr=dfEffChi2ndf[f'trkeffErr_{tt}'],
        marker="s", markersize=5, color="k", linestyle="--",
        label=f"{system(tt)} {algorithm(tt)}\nefficiency"
    )

    ax2 = axs[ax_i, ax_j].twinx()

    l2 = ax2.errorbar(
        x=chi2ndf[_att]["x"], y=chi2ndf[_att]["y"],
        xerr=chi2ndf[_att]["ex"], yerr=chi2ndf[_att]["ey"],
        label=f'{system(_att)} {algorithm(_att)}\n$\chi^2/ndf\;distribution$', color='blue',
        marker='o', markersize=5
    )

    ax2.tick_params(axis='y', labelcolor='blue', size=16)

    axs[ax_i, ax_j].set_title(f"{system(tt)} {algorithm(tt)} (run 7080)", fontsize=16)
    axs[ax_i, ax_j].set_ylabel(f"Tracking efficiency", fontsize=16)
    axs[ax_i, ax_j].set_xlabel(f"$(\chi^2/ndf)_{{max}}$ for tagging track ({sys(_att).upper()} {alg(tt).upper()})", fontsize=16)
    axs[ax_i, ax_j].tick_params(labelsize=16)
    axs[ax_i, ax_j].set_ylim(0, 1)
    axs[ax_i, ax_j].grid(axis="both")

    # Set the label color for the right axis (ax2) to blue
    ax2_label = f"{sys(_att).upper()} {alg(_att).upper()} " + r'$\chi^2 / ndf\;$'+'distribution'
    ax2.set_ylabel(ax2_label, color='navy')
    ax2.tick_params(labelsize=16)

    formatter = mticker.ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((3, 3))
    ax2.yaxis.set_major_formatter(formatter)
    ax2.yaxis.get_offset_text().set_size(16)

    h1, l1 = axs[ax_i, ax_j].get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()

    axs[ax_i, ax_j].legend(h1 + h2, l1 + l2, loc="center right", fontsize=16)

plt.tight_layout()
filename=f"{mfout}/trkeff_varChi2ndf.pdf"
plt.savefig(filename, format="pdf")
filename=f"{mfout}/trkeff_varChi2ndf.png"
plt.savefig(filename, format="png")
plt.show()
