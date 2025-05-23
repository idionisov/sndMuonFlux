import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import uproot
import ROOT
from sndUtils import att, alg, sys, algorithm, system, nName
from ddfRoot import getAsPandas, getAsNumpy
from mplStyling import addMplLabel


gr = {"xz": {}, "yz": {}, "n": {}, "e": {}}
gr["n"] = {
    {a:
        getAsPandas(
            uproot.open(
                f"/eos/user/i/idioniso/mfout/trkeff/trkeff-E/trkeff_muGun_{a}.root:gr_eff_n_{tt}_muGun.{a}"
            )
        ) for a in ("rt", "tc")
    } for tt in (1, 11, 3, 13)
}

for ang in ("xz", "yz"):
    gr[ang] = {}
    for tt in (1, 11, 3, 13):
        gr[ang][tt] = getAsPandas(
            uproot.open(
                f"/eos/user/i/idioniso/mfout/MuonFlux.root:TrackingEfficiency/Run7080/gr_eff_yz_{tt}_7080.tc"
            )
        )



nDistr = {
    getAsPandas(
        uproot.open(
            f"/eos/user/i/idioniso/mfout/MuonFlux.root:TrackProperties/Run7080/h_n_{tt}_7080"
        )
    ) for tt in (1, 11, 3, 13)
}

for tt in nDistr:
    nDistr[tt]["y"] /= nDistr[tt]["y"].sum()
    nDistr[tt]["ey"] /= nDistr[tt]["y"].sum()



def plot_TrkeffVsN_TTvsMC int, pd.DataFrame, nDi pd.DataFrame):
    fig, ax1 = plt.subplots(figsize=(11, 5))
    ax2 = ax1.twinx()

    rt = gr["n"][tt]["rt"].query("y > 0 and eyh < 0.9 and eyl < 0.9")
    tc = gr["n"][tt]["tc"].query("y > 0 and eyh < 0.9 and eyl < 0.9")
    nD = nDistr[tt]

    step_plot, = ax2.step(
        nD['x']-0.5, nD['y'],
        where="mid",
        color="k",
        marker="",
        label=f"{nName(tt)} distribution"
    )
    ax2.set_ylim(0, 1.2*nD['y'].max())
    ax2.set_ylabel(f"{nName(tt)} Distribution", fontsize=16, fontname="serif", color="k")
    ax2.tick_params(axis='y', labelsize=16, labelcolor="k")

    rt_plot = ax1.errorbar(
        rt["x"]-0.5, rt["y"],
        yerr=(rt["eyl"], rt["eyh"]),
        label="'True' efficiency",
        marker="s",
        markersize=4,
        markerfacecolor="white",
        capsize=0,
        linestyle='dashed',
        linewidth=2
    )
    tc_plot = ax1.errorbar(
        tc["x"]-0.5, tc["y"],
        yerr=(tc["eyl"], tc["eyh"]),
        label="Estimated efficiency",
        marker="^",
        markersize=4,
        markerfacecolor="white",
        linestyle=":",
        capsize=0,
        linewidth=2
    )

    xmax = 30, 70 12, 25}
    ax1.set_xlim(np.min(rt["x"] - 1.05), xmax[tt]+0.55)

    ax1.set_xlabel(f"Number of {nName(tt)}", fontname="serif", fontsize=16)
    ax1.tick_params(axis='x', labelsize=16)

    ax1.set_ylim(0.3, 1.025)
    ax1.set_ylabel(f"Tracking efficiency", fontsize=16, fontname="serif")
    ax1.tick_params(axis='y', labelsize=16)

    handles = [rt_plot, tc_plot, step_plot]
    labels = [h.get_label() for h in handles]
    legPos = {1:"upper right", "upper right" "upper right", "lower right"}
    anchor = (1, 0.5), (1, 0.5) (1, 0.9), (1, 0.15)}
    ax1.legend(handles, labels, prop={'fami 'serif', 'si 15}, bbox_to_anchor=anchor[tt], loc=legPos[tt])

    addMplLabel(mainText="SND@LHC", extraText="Simulation", fontsize=16)
    plt.tight_layout()

    plt.savefig(f"/eos/user/i/idioniso/mfout/png/trkeff-{sys(tt)}.{alg(tt)}-muGun.N.png")
    plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/trkeff-{sys(tt)}.{alg(tt)}-muGun.N.pdf")




def plot_TrkeffVsE_TTvsMC int, pd.DataFrame):
    plt.figure(figsize=(11, 5))

    df = gr["e"][tt]

    plt.errorbar(
        df["e"], df[f"trkeff_{tt}_rt"],
        yerr=(df[f"trkeffErrLow_{tt}_rt"], df[f"trkeffErrUp_{tt}_rt"]),
        label="'True' efficiency",
        markersize=7,
        markerfacecolor="white",
        capsize=3
    )
    plt.errorbar(
        df["e"], df[f"trkeff_{tt}_tc"],
        yerr=(df[f"trkeffErrLow_{tt}_tc"], df[f"trkeffErrUp_{tt}_tc"]),
        label="Measured efficiency",
        markersize=7,
        markerfacecolor="white",
        capsize=3
    )

    plt.xscale("log")
    plt.xlabel("Muon Energy $E_{\mu}$ [GeV]", fontname="serif", fontsize=16)
    plt.xticks(fontsize=16, fontname="serif")

    plt.ylabel("Tracking efficiency", fontsize=16, fontname="serif")
    plt.yticks(fontsize=16, fontname="serif")
    plt.ylim(0.3, 1.025)

    plt.legend(prop={'fami 'serif', 'si 15})
    addMplLabel(mainText="SND@LHC", extraText="Simulation", fontsize=16)
    plt.tight_layout()

    plt.savefig(f"/eos/user/i/idioniso/mfout/png/trkeff-{sys(tt)}.{alg(tt)}-muGun.E.png")
    plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/trkeff-{sys(tt)}.{alg(tt)}-muGun.E.pdf")




def plot_TrkeffVsXZ(gr, alg: str = "ht", ang: str = "xz"):
    plt.figure(figsize=(11, 5))

    sftt, dstt = 1, 3
    if alg == "ht":
        sftt+=10
        dstt+=10

    sf = gr[ang][sftt]
    ds = gr[ang][dstt]

    sf = sf.query("y > 0 and eyl < 0.9 and eyh < 0.9")
    ds = ds.query("y > 0 and eyl < 0.9 and eyh < 0.9")

    plt.errorbar(
        sf["x"], sf["y"],
        xerr=(sf["exl"], sf["exh"]),
        yerr=(sf["eyl"], sf["eyh"]),
        label=f"{system(sftt)} {algorithm(sftt)}",
        marker="s",
        markersize=9,
        markerfacecolor="white",
        capsize=3,
        linestyle=""
    )
    plt.errorbar(
        ds["x"], ds["y"],
        xerr=(ds["exl"], ds["exh"]),
        yerr=(ds["eyl"], ds["eyh"]),
        label=f"{system(dstt)} {algorithm(dstt)}",
        marker="^",
        markersize=9,
        markerfacecolor="white",
        capsize=3,
        linestyle=""
    )

    plt.xticks(fontsize=16, fontname="serif")
    plt.xlabel(f"$\\theta_{{{ang.upper()}}}$ [$mrad$]", fontname="serif", fontsize=16)

    plt.yticks(fontsize=16, fontname="serif")
    plt.ylabel("Tracking efficiency", fontsize=16, fontname="serif")
    plt.ylim(0, 1.025)

    plt.legend(prop={'family': 'serif', 'size': 15})
    addMplLabel(mainText="SND@LHC", extraText="Run 7080", fontsize=16)

    plt.tight_layout()
    plt.savefig(f"/eos/user/i/idioniso/mfout/png/trkeff.{ang.upper()}-{alg}.png")
    plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/trkeff.{ang.upper()}-{alg}.pdf")
