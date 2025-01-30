import ROOT
import pandas as pd
import numpy as np
import uproot
from sndUtils import sys, alg, system, algorithm, att
import matplotlib.ticker as mticker
from ddfRoot import getPandasFromUprootTH1, getNumpyFromUprootTH2
import matplotlib.pyplot as plt
plt.style.use("root")

mfout = "/eos/user/i/idioniso/mfout"



def plotTrkeffVarChi2ndf(
    inputTrkeff: str = "/eos/user/i/idioniso/mfout/trkeff/trkeff-chi2ndfMax/trkeff-chi2ndf.csv",
    inputPrpts: str = "/eos/user/i/idioniso/mfout/MuonFlux.root",
    chi2ndf: dict = {
        1:  [5,   10,  15,  20, 40,  60, 100],
        11: [8,   12,  16,  25, 50,  80, 120],
        3:  [2,   3,   4,   5,  7.5, 10, 15 ],
        13: [3.5, 4.5, 5.5, 7,  9,   12, 15 ]
    }
):

    dfEffChi2ndf = pd.read_csv(inputTrkeff)

    chi2ndf = {
        tt: getPandasFromUprootTH1(
            uproot.open(f"{inputPrpts}:TrackProperties/Run7080/h_chi2ndf_{tt}_7080")
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

    axs[0, 0].set_xlim(0, 20.5)

    plt.tight_layout()
    plt.savefig("/eos/user/i/idioniso/mfout/trkeff_varChi2ndf.pdf", format="pdf")
    plt.savefig("/eos/user/i/idioniso/mfout/trkeff_varChi2ndf.png", format="png")
    plt.show()






def plotTrkeffE(
    input: str = "/eos/user/i/idioniso/mfout/trkeff/trkeff-E/trkeff-E.csv",
):
    df = pd.read_csv(input)


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
            df["e"] * (i_tt/75 + 1),
            df[f"trkeff_{tt}_tc"],
            yerr=df[f"trkeffErr_{tt}_tc"],
            label=f"{system(tt)} {algorithm(tt)}",
            marker="s", markersize=5
        )

        ax1.fill_between(df["e"],
            df[f"trkeff_{tt}_tc"] - df[f"trkeffErr_{tt}_tc"]/2,
            df[f"trkeff_{tt}_tc"] + df[f"trkeffErr_{tt}_tc"]/2,
            alpha=0.375
        )




        ax2.errorbar(
            df["e"] * (i_tt/75 + 1),
            df[f"trkeff_{tt}_rt"],
            yerr=df[f"trkeffErr_{tt}_rt"],
            label=f"{system(tt)} {algorithm(tt)}",
            marker="s", markersize=5
        )

        ax2.fill_between(df["e"],
            df[f"trkeff_{tt}_rt"] - df[f"trkeffErr_{tt}_rt"]/2,
            df[f"trkeff_{tt}_rt"] + df[f"trkeffErr_{tt}_rt"]/2,
            alpha=0.375
        )

    ax1.text(0.1, 0.15, "TC method (Measured efficiency)", transform=ax1.transAxes)
    ax2.text(0.1, 0.15, "RT method ('True' efficiency)", transform=ax2.transAxes)

    ax1.legend()
    ax2.legend()

    plt.savefig("/eos/user/i/idioniso/mfout/trkeff_E_em.pdf", format="pdf")
    plt.savefig("/eos/user/i/idioniso/mfout/trkeff_E_em.png", format="png")

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
            axs[ax_i][ax_j].errorbar(df["e"],
                df[f"trkeff_{tt}_{em}"],
                yerr=df[f"trkeffErr_{tt}_{em}"],
                marker="s", markersize=5,
                label=em_name
            )

            axs[ax_i][ax_j].fill_between(df["e"],
                df[f"trkeff_{tt}_{em}"] - df[f"trkeffErr_{tt}_{em}"]/2,
                df[f"trkeff_{tt}_{em}"] + df[f"trkeffErr_{tt}_{em}"]/2,
                alpha=0.5
            )

            axs[ax_i][ax_j].tick_params(labelsize=16)
            axs[ax_i][ax_j].legend(fontsize=12)
    plt.tight_layout()

    plt.savefig("/eos/user/i/idioniso/mfout/trkeff_E_tt.pdf", format="pdf")
    plt.savefig("/eos/user/i/idioniso/mfout/trkeff_E_tt.png", format="png")

    plt.show()




def plotTrkeffXY(
    run: int = 7080,
    xmin: float = -75.,
    xmax: float =  5.,
    ymin: float = -5.,
    ymax: float =  75.
):
    eff, effErr, xy, xyEdgesX, xyEdgesY, X, Y = {}, {}, {}, {}, {}, {}, {}
    for tt in (1, 11, 3, 13):
        fin = uproot.open("/eos/user/i/idioniso/mfout/trkeffDataTC.root")
        eff[tt] = fin[f"Run{run}/eff_{tt}_data.tc/eff"].array()[0]
        effErr[tt] = fin[f"Run{run}/eff_{tt}_data.tc/effErr"].array()[0]

        xy[tt], xyEdgesX[tt], xyEdgesY[tt] = getNumpyFromUprootTH2(fin[f"Run{run}/eff_x.y_{tt}_{run}.tc"])
        X[tt], Y[tt] = np.meshgrid(xyEdgesX[tt], xyEdgesY[tt])



    fig, axs = plt.subplots(2, 2, figsize=(17, 8))

    for i_tt, tt in enumerate((1, 11, 3, 13)):
        ax_i, ax_j = divmod(i_tt, 2)

        pcm = axs[ax_i, ax_j].pcolormesh(X[tt], Y[tt], xy[tt], cmap='root', vmin=0, vmax=1)
        cbar = fig.colorbar(pcm, ax=axs[ax_i, ax_j], label='Tracking efficiency')
        cbar.set_label('Tracking efficiency', fontsize=18)
        cbar.ax.tick_params(labelsize=16)

        axs[ax_i, ax_j].set_title(f"{system(tt)} {algorithm(tt)} (run {run})", fontsize=18)
        axs[ax_i, ax_j].set_xlabel("$X_{ref}$ [cm]", fontsize=18)
        axs[ax_i, ax_j].set_ylabel("$Y_{ref}$ [cm]", fontsize=18)
        axs[ax_i, ax_j].grid(axis="both")
        axs[ax_i, ax_j].tick_params(axis='both', which='major', labelsize=16)

        axs[ax_i][ax_j].text(-65, 61, f"$\\varepsilon \\approx {eff[tt]:.03f} ± {effErr[tt]:.03f}$", fontsize=16, fontweight='bold')


    plt.tight_layout()

    fname=f"/eos/user/i/idioniso/mfout/trkeff_xy_run{run}_tc.pdf"
    plt.savefig(fname, format="pdf")
    fname=f"/eos/user/i/idioniso/mfout/trkeff_xy_run{run}_tc.png"
    plt.savefig(fname, format="png")

    plt.show()
