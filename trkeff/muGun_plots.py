import ROOT
import numpy as np
from sndUtils import system, algorithm
import matplotlib.pyplot as plt
plt.style.use("root")

mfout = "/EOS/user/i/idioniso/mfout"
energies = [10, 12.5, 15, 17.5, 20, 25, 30, 35, 55, 100, 200, 300, 600, 800, 1010]

eff      = {a: {tt: np.zeros(len(energies), dtype=np.float64) for tt in (1, 11, 3, 13)} for a in ("rt", "tc")}
deff_up  = {a: {tt: np.zeros(len(energies), dtype=np.float64) for tt in (1, 11, 3, 13)} for a in ("rt", "tc")}
deff_low = {a: {tt: np.zeros(len(energies), dtype=np.float64) for tt in (1, 11, 3, 13)} for a in ("rt", "tc")}
E        = np.array(energies, dtype=np.float64)

for i_e, e in enumerate(energies):
    input_dir = f"{mfout}/trkeff_varE"
    input = {
        "rt": f"{input_dir}/trkeff_muGun.{e}GeV_rt.root",
        "tc": f"{input_dir}/trkeff_muGun.{e}GeV_tc.root"
    }

    trees = {}
    for a in ("rt", "tc"):
        trees[a] = {}

        for tt in (1, 11, 3, 13):
            trees[a][tt] = ROOT.TChain(f"eff_{tt}_muGun.{a}")
            trees[a][tt].Add(input[a])
            trees[a][tt].GetEntry(0)


            eff[a][tt][i_e]      = trees[a][tt].eff
            deff_up[a][tt][i_e]  = trees[a][tt].effErrUp
            deff_low[a][tt][i_e] = trees[a][tt].effErrLow



fig = plt.figure(figsize=(9, 8))

ax1 = fig.add_axes([0.1, 0.5, 0.9, 0.4], xticklabels=[])
ax2 = fig.add_axes([0.1, 0.1, 0.9, 0.4])

ax1.set_ylabel("Tracking Efficiency")
ax1.set_ylim(0.25, 1)
ax1.set_xlim(9.5, 1100)
ax1.set_xscale("log")
ax1.grid()

ax2.set_xlabel("Muon energy [GeV]")
ax2.set_ylabel("Tracking efficiency")
ax2.set_ylim(0.25, 1)
ax2.set_xlim(9.5, 1100)
ax2.set_xscale("log")
ax2.grid()


for i_tt, tt in enumerate((1, 11, 3, 13)):
    ax1.errorbar(
        E * (i_tt/75 + 1),
        eff["tc"][tt],
        yerr=(deff_low["tc"][tt], deff_up["tc"][tt]),
        label=f"{system(tt)} {algorithm(tt)}",
        marker="s", markersize=5
    )

    ax1.fill_between(E,
        eff["tc"][tt] - deff_low["tc"][tt],
        eff["tc"][tt] + deff_up["tc"][tt],
        alpha=0.375
    )




    ax2.errorbar(
        E * (i_tt/75 + 1),
        eff["rt"][tt],
        yerr=(deff_low["rt"][tt], deff_up["rt"][tt]),
        label=f"{system(tt)} {algorithm(tt)}",
        marker="s", markersize=5
    )

    ax2.fill_between(E,
        eff["rt"][tt] - deff_low["rt"][tt],
        eff["rt"][tt] + deff_up["rt"][tt],
        alpha=0.375
    )

ax1.text(0.1, 0.15, "TC method (Measured efficiency)", transform=ax1.transAxes)
ax2.text(0.1, 0.15, "RT method ('True' efficiency)", transform=ax2.transAxes)

ax1.legend()
ax2.legend()
plt.show()














fig, axs = plt.subplots(2, 2, figsize=(12, 6))

for i_tt, tt in enumerate([1, 11, 3, 13]):
    ax_i, ax_j = divmod(i_tt, 2)

    axs[ax_i][ax_j].set_title(f"{system(tt)} {algorithm(tt)}", fontsize=14)
    axs[ax_i][ax_j].set_ylabel("Tracking Efficiency", fontsize=12)
    axs[ax_i][ax_j].set_xlabel("Muon Energy [GeV]", fontsize=12)
    axs[ax_i][ax_j].grid(axis="y")
    axs[ax_i][ax_j].set_ylim(0.25, 1.)
    axs[ax_i][ax_j].set_xlim(9.6, 1100)
    axs[ax_i][ax_j].set_xscale("log")


    for em, em_name in zip(["tc", "rt"], ["TC Method", "RT Method"]):
        axs[ax_i][ax_j].errorbar(E,
            eff[em][tt],
            yerr=(deff_low[em][tt], deff_up[em][tt]),
            marker="s", markersize=5,
            label=em_name
        )

        axs[ax_i][ax_j].fill_between(E,
            eff[em][tt] - deff_low[em][tt],
            eff[em][tt] + deff_up[em][tt],
            alpha=0.5
        )

        axs[ax_i][ax_j].tick_params(labelsize=16)
        axs[ax_i][ax_j].legend(fontsize=12)
plt.tight_layout()

fname=f"{mfout}/pdf/trkeff_e-TCvsRT_pGun.pdf"
plt.savefig(fname, format="pdf")

plt.show()
print(fname)
