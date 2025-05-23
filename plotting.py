import ROOT
from mplhep.utils import markers
import numpy as np
import matplotlib.pyplot as plt
from numpy.ma.core import log10
import pandas as pd
from pandas.core.construction import ma
from scipy.interpolate import Rbf
from scipy.stats import linregress
import uproot

from ddfRoot import getAsPandas, getAsNumpy
from sndUtils import sys, alg, system, algorithm, nName
from mplStyling import addMplLabel



df = pd.read_csv("/eos/user/i/idioniso/mfout/MuonFlux.csv")


def plot_FluxVsRuns(df, lumi: str = "eos"):
    fig, ax = plt.subplots(1, 2, sharey=True, figsize=(11, 4), gridspec_kw={'width_ratios': [1, 1]})

    df1 = df.query("Run <  10000")
    df2 = df.query("Run >=  10000")
    addMplLabel(mainText="SND@LHC", extraText="2023/2024 PbPb Data", fontsize=16, ax=ax[0])

    for tt in (1, 11, 3, 13):
        ax[0].errorbar(df1["Run"], df1[f"Flux{tt}_{lumi}"], yerr=df1[f"FluxErr{tt}_{lumi}"], label=f"{system(tt)} {algorithm(tt)}", markersize=7, markerfacecolor="white", capsize=3)
        ax[1].errorbar(df2["Run"], df2[f"Flux{tt}_{lumi}"], yerr=df2[f"FluxErr{tt}_{lumi}"], label=f"{system(tt)} {algorithm(tt)}", markersize=7, markerfacecolor="white", capsize=3)


    ax[0].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[0].tick_params(labelright=False, labelsize=16, labelfontfamily="serif")
    ax[1].tick_params(labelleft=False, labelsize=16, labelfontfamily="serif")


    d = .015  # Size of the diagonal line
    kwargs = dict(transform=ax[0].transAxes, color='k', clip_on=False)
    ax[0].plot((1-d, 1+d), (-d, +d), **kwargs)  # Top-right
    ax[0].plot((1-d, 1+d), (1-d, 1+d), **kwargs)  # Bottom-right
    ax[0].set_xlim(None, 7346)
    ax[0].set_ylim(0, 27e3)
    ax[0].set_xlabel("Run number", fontsize=16, fontname="serif")
    ax[0].set_ylabel("Muon Flux $\Phi_{\mu}$ [$nb/cm^{2}$]", fontsize=16, fontname="serif")

    kwargs.update(transform=ax[1].transAxes)
    ax[1].plot((-d, +d), (-d, +d), **kwargs)  # Top-left
    ax[1].plot((-d, +d), (1-d, 1+d), **kwargs)  # Bottom-left
    ax[1].set_xlim(10246, None)
    ax[1].set_ylim(0, 27e3)
    ax[1].set_xlabel("Run number", fontsize=16, fontname="serif")

    plt.legend(prop={'family': 'serif', 'size': 15})
    plt.tight_layout()

    plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/Flux.Run_{lumi}.pdf")
    plt.savefig(f"/eos/user/i/idioniso/mfout/png/Flux.Run_{lumi}.png")



def _plot_FluxVsRuns(lumi: str = "eos"):
    plt.figure(figsize=(11, 5))
    for tt in (1, 11, 3, 13):
        plt.errorbar(
            df["run"], df[f"Flux{tt}_{lumi}"],
            yerr=df[f"FluxErr{tt}_{lumi}"],
            label=f"{system(tt)} {algorithm(tt)}",
            lw=2
        )
    plt.plot(
        [df['run'].iloc[0], df['run'].iloc[-1]],
        [17775.731025833335, 17775.731025833335],
        color="k",
        linestyle="--",
        lw=2
    )
    plt.fill_between(
        [df['run'].iloc[0], df['run'].iloc[-1]],
        [17775.731025833335-826.0352532599704/2, 17775.731025833335-826.0352532599704/2],
        [17775.731025833335+826.0352532599704/2, 17775.731025833335+826.0352532599704/2],
        color="k",
        alpha=0.1
    )
    plt.xticks(fontsize=16, fontname="serif")
    plt.xlabel("Run number", fontname="serif", fontsize=16)
    plt.yticks(fontsize=16, fontname="serif")
    plt.ylabel("Muon Flux $\Phi_{\mu}$ [$nb/cm^{2}$]", fontsize=16, fontname="serif")
    plt.legend(prop={'family': 'serif', 'size': 15})
    addMplLabel(mainText="SND@LHC", extraText="2023 PbPb Data", fontsize=16)
    plt.tight_layout()
    plt.savefig(f"/eos/user/i/idioniso/mfout/png/Flux.Run_{lumi}.png")
    plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/Flux.Run_{lumi}.pdf")


def plot_FluxVsRuns_EosAndSt(tt: int = 11):
    fig, ax = plt.subplots(1, 2, sharey=True, figsize=(11, 4), gridspec_kw={'width_ratios': [1, 1]})

    df1 = df.query("run <  10000")
    df2 = df.query("run >=  10000")
    addMplLabel(mainText="SND@LHC", extraText="2023/2024 PbPb Data", fontsize=16, ax=ax[0])

    for lumi, Lumi in zip(["eos", "st"], ["eos", "LHC supertable"]):
        ax[0].errorbar(df1["run"], df1[f"Flux{tt}_{lumi}"], yerr=df1[f"FluxErr{tt}_{lumi}"], label=f"$\mathcal{{L}}_{{int}}$ from {Lumi}", markersize=7, markerfacecolor="white", capsize=3)
        ax[1].errorbar(df2["run"], df2[f"Flux{tt}_{lumi}"], yerr=df2[f"FluxErr{tt}_{lumi}"], label=f"$\mathcal{{L}}_{{int}}$ from {Lumi}", markersize=7, markerfacecolor="white", capsize=3)


    ax[0].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[0].tick_params(labelright=False, labelsize=16, labelfontfamily="serif")
    ax[1].tick_params(labelleft=False, labelsize=16, labelfontfamily="serif")


    d = .015  # Size of the diagonal line
    kwargs = dict(transform=ax[0].transAxes, color='k', clip_on=False)
    ax[0].plot((1-d, 1+d), (-d, +d), **kwargs)  # Top-right
    ax[0].plot((1-d, 1+d), (1-d, 1+d), **kwargs)  # Bottom-right
    ax[0].set_xlim(None, 7346)
    ax[0].set_ylim(0, 27e3)
    ax[0].set_xlabel("Run number", fontsize=16, fontname="serif")
    ax[0].set_ylabel("Muon Flux $\Phi_{\mu}$ [$nb/cm^{2}$]", fontsize=16, fontname="serif")

    kwargs.update(transform=ax[1].transAxes)
    ax[1].plot((-d, +d), (-d, +d), **kwargs)  # Top-left
    ax[1].plot((-d, +d), (1-d, 1+d), **kwargs)  # Bottom-left
    ax[1].set_xlim(10246, None)
    ax[1].set_ylim(0, 27e3)
    ax[1].set_xlabel("Run number", fontsize=16, fontname="serif")

    plt.legend(prop={'family': 'serif', 'size': 15})
    plt.tight_layout()

    plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/Flux{tt}.Run_LumiDiff.pdf")
    plt.savefig(f"/eos/user/i/idioniso/mfout/png/Flux{tt}.Run_LumiDiff.png")

def _plot_FluxVsRuns_EosAndSt(tt: int):
    plt.figure(figsize=(11, 5))
    plt.errorbar(
        df["run"], df[f"Flux{tt}_eos"],
        yerr=df[f"FluxErr{tt}_eos"],
        marker="s",
        markersize=5,
        markerfacecolor="white",
        label="$\mathcal{L}_{int}$ from eos",
        lw=2
    )
    plt.errorbar(
            df["run"], df[f"Flux{tt}_st"],
            yerr=df[f"FluxErr{tt}_st"],
            marker="^",
            markersize=5,
            markerfacecolor="white",
            label="$\mathcal{L}_{int}$ from supertable",
            linestyle="--",
            lw=2
        )
    plt.xticks(fontsize=16, fontname="serif")
    plt.xlabel("Run number", fontname="serif", fontsize=16)
    plt.yticks(fontsize=16, fontname="serif")
    plt.ylabel("Muon Flux $\Phi_{\mu}$ [$nb/cm^{2}$]", fontsize=16, fontname="serif")
    plt.legend(prop={'family': 'serif', 'size': 15})
    addMplLabel(mainText="SND@LHC", extraText="2023 PbPb Data", fontsize=16)
    plt.tight_layout()
    plt.savefig(f"/eos/user/i/idioniso/mfout/png/Flux{tt}.Run_LumiDiff.png")
    plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/Flux{tt}.Run_LumiDiff.pdf")


def plot_LumiVsT_Flux(tt: int, lumi: str = "eos"):
    plt.figure(figsize=(11, 5))

    plt.scatter(df[f"Lumi{lumi.capitalize()}"], df['Stable B. time [h]'], c=df[f"Flux{tt}_{lumi}"], cmap="jet", marker="s", s=45*np.ones(len(df)))
    plt.colorbar(label="Muon Flux $\Phi_{\mu}$ [$nb/cm^{2}$]")
    plt.xticks(fontsize=16, fontname="serif")
    plt.xlabel("Integrated Luminosiy $\mathcal{L}_{int}$ [$nb^{-1}$]", fontname="serif", fontsize=16)
    plt.yticks(fontsize=16, fontname="serif")
    plt.yscale("log")
    plt.ylabel("Stable beams duration [h]", fontsize=16, fontname="serif")
    addMplLabel(mainText="SND@LHC", extraText="2023 PbPb Data", fontsize=16)
    plt.tight_layout()
    plt.savefig(f"/eos/user/i/idioniso/mfout/png/Lumi.T.Flux{tt}.png")
    plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/Lumi.T.Flux{tt}.pdf")


def plot_LumiVsT_ColIP1(tt: int, lumi: str = "eos"):
    plt.figure(figsize=(11, 5))

    plt.scatter(df[f"Lumi{lumi.capitalize()}"], df['Stable B. time [h]'], c=df['# Coll. IP 1/5'], cmap="jet", marker="s", s=45*np.ones(len(df)))
    plt.colorbar(label="Muon Flux $\Phi_{\mu}$ [$nb/cm^{2}$]")
    plt.xticks(fontsize=16, fontname="serif")
    plt.xlabel("Integrated Luminosiy $\mathcal{L}_{int}$ [$nb^{-1}$]", fontname="serif", fontsize=16)
    plt.yticks(fontsize=16, fontname="serif")
    plt.yscale("log")
    plt.ylabel("Stable beams duration [h]", fontsize=16, fontname="serif")
    addMplLabel(mainText="SND@LHC", extraText="2023 PbPb Data", fontsize=16)
    plt.tight_layout()
    plt.savefig(f"/eos/user/i/idioniso/mfout/png/Lumi.T.Flux{tt}.png")
    plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/Lumi.T.Flux{tt}.pdf")



def plot_FluxOmega(tt: int, lumi: str = "eos"):
    plt.figure(figsize=(11, 6))
    omega = ((df['Intensity B1 [1e11 ppb]'] + df['Intensity B2 [1e11 ppb]']) * df['# Coll. IP 1/5'] * np.log(df['Stable B. time [h]']))/df[f"Lumi{lumi.capitalize()}"]

    plt.scatter(df[f"Lumi{lumi.capitalize()}"], omega, c=df[f"Flux{tt}_{lumi}"], cmap="jet", marker="s", s=45*np.ones(len(df)))
    plt.colorbar(label="Muon Flux $\Phi_{\mu}$ [$nb/cm^{2}$]")
    plt.xticks(fontsize=16, fontname="serif")
    plt.xlabel("Integrated Luminosiy $\mathcal{L}_{int}$ [$nb^{-1}$]", fontname="serif", fontsize=16)
    plt.yticks(fontsize=16, fontname="serif")
    plt.ylabel("$\Omega [1e11 ppb h nb]$", fontsize=16, fontname="serif")
    addMplLabel(mainText="SND@LHC", extraText="2023 PbPb Data", fontsize=16)
    plt.tight_layout()
    plt.savefig(f"/eos/user/i/idioniso/mfout/png/Lumi.Flux{tt}.Omega.png")
    plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/Lumi.Flux{tt}.Omega.pdf")



def plot_FluxZeta(tt: int, lumi: str = "eos"):
    plt.figure(figsize=(11, 6))
    omega = df['# Coll. IP 1/5'] * np.log(df['Stable B. time [h]']) / df[f"Flux{tt}_{lumi}"]

    plt.scatter(df[f"Lumi{lumi.capitalize()}"], omega, c=df[f"Flux{tt}_{lumi}"], cmap="jet", marker="s", s=45*np.ones(len(df)))
    plt.colorbar(label="Muon Flux $\Phi_{\mu}$ [$nb/cm^{2}$]")
    plt.xticks(fontsize=16, fontname="serif")
    plt.xlabel("Integrated Luminosiy $\mathcal{L}_{int}$ [$nb^{-1}$]", fontname="serif", fontsize=16)
    plt.yticks(fontsize=16, fontname="serif")
    plt.ylabel("$N_{IP1} \; \ln (T_{sb}) / \Phi_{\mu}$", fontsize=16, fontname="serif")
    addMplLabel(mainText="SND@LHC", extraText="2023 PbPb Data", fontsize=16)
    plt.tight_layout()
    plt.savefig(f"/eos/user/i/idioniso/mfout/png/Lumi.Flux{tt}.Omega.png")
    plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/Lumi.Flux{tt}.Omega.pdf")



def plot_FluxVsSBtime(df, lumi: str = "eos", year: int = 2023):

    if year==2023:
        df_new = df.query("Run < 10000")
    elif year==2024:
        df_new = df.query("Run >= 10000")
    else:
        return -9999

    df_new = df_new.sort_values(by='Stable B. time [h]')
    plt.figure(figsize=(12, 6))

    for tt in (1, 11, 3, 13):
        x = df_new['Stable B. time [h]']
        y = df_new[f"Flux{tt}_{lumi}"]
        ey = df_new[f"FluxErr{tt}_{lumi}"]
        plt.errorbar(x, y, yerr=ey, label=f"{system(tt)} {algorithm(tt)}", lw=2, markersize=7, markerfacecolor="white", capsize=3, alpha=0.8)

    plt.xlabel('Stable B. time (h)', fontsize=20)
    plt.ylabel("Muon flux ($1/nb$)", fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylim(0, 25000)

    plt.legend(prop={"size": 18})
    plt.tight_layout()
    plt.savefig(f"/eos/user/i/idioniso/mfout/png/Flux{lumi.capitalize()}.SBtime{year}.png")
    plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/Flux{lumi.capitalize()}.SBtime{year}.pdf")



def plot_AllRuns(tt: int, lumi: str = "eos"):
    if lumi.lower() in ["eos"]:
        lumiColumn = "LumiEosOld"
    elif lumi.lower() in ["st", "supertable"]:
        lumiColumn = 'ATLAS Int. Lumi [1/nb]'
    else:
        return None

    plt.figure(figsize=(12, 6))
    sc = plt.scatter(x=mf["Stable B. time [h]"], y=mf[f"Flux{tt}_{lumi}"], c=np.log(mf["Stable B. time [h]"])/mf[lumiColumn], s=np.sqrt(mf['# Coll. IP 1/5'])*4, cmap="viridis", marker="s", alpha=0.65)

    cbar = plt.colorbar(sc)
    cbar.set_label(
        "$\\frac{\\log(T_{sb})}{\\mathcal{L}_{int}}$",
        fontsize=20
    )
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel("Stable B. time $T_{sb}$ (h)", fontsize=20)
    plt.ylabel("Muon Flux $\Phi_{\mu}$ ($nb/cm^{2}$)", fontsize=20)
    addMplLabel(mainText="SND@LHC", extraText="2023/2024 PbPb Data", fontsize=18)

    plt.tight_layout()
    plt.show()


def plot_AllRuns2(tt: int, lumi: str = "eos"):
    if lumi.lower() in ["eos"]:
        lumiColumn = "LumiEosOld"
    elif lumi.lower() in ["st", "supertable"]:
        lumiColumn = 'ATLAS Int. Lumi [1/nb]'
    else:
        return None

    df_low = mf[mf["Run"] < 10000]
    df_high = mf[mf["Run"] >= 10000]

    plt.figure(figsize=(12, 6))

    for df, label, M in zip(
        [df_low, df_high],
        ["2023", "2024"],
        ["s", "o"]
    ):
        sc = plt.scatter(
            x=df[lumiColumn],
            y=df['# Coll. IP 1/5'] * np.log(df["Stable B. time [h]"]) / df[f"Flux{tt}_{lumi}"],
            c=np.log(df[f"Flux{tt}_{lumi}"]),
            s=np.log(df["Stable B. time [h]"]) * 100,
            cmap="viridis",
            marker=M,
            alpha=0.65,
            label=label
        )

    cbar = plt.colorbar(sc)
    cbar.set_label("Muon flux $\Phi_{\mu}$ ($nb/cm^{2}$)", fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel("$\mathcal{L}_{int}$ ($1/nb$)", fontsize=20)
    plt.ylabel("$N_{IP1} \; \ln(T_{sb})/\Phi_{\mu}$", fontsize=20)
    addMplLabel(mainText="SND@LHC", extraText="2023/2024 PbPb Data", fontsize=18)

    plt.legend(fontsize=18)
    plt.tight_layout()
    plt.savefig(f"/eos/user/i/idioniso/mfout/png/AllRuns2-Flux{tt}-{lumi.capitalize()}.png")
    plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/AllRuns2-Flux{tt}-{lumi.capitalize()}.pdf")





def _plot_AllRunsOmega(df0, I: int, tt: int, lumi: str = "eos", sb_cut: float = 0, saveFig: bool = True):
    df = df0.query(f"`Stable B. time [h]` > {sb_cut}")
    df = df.query("Flux11_eos < 50000")
    if lumi.lower() in ["eos"]:
        lumiColumn = "LumiEosOld"
    elif lumi.lower() in ["st", "supertable"]:
        lumiColumn = 'ATLAS Int. Lumi [1/nb]'
    else:
        return None

    df_low = df[df["Run"] < 10000]
    df_high = df[df["Run"] >= 10000]


    flux11eos23, fluxErr11eos23 = getMeanFlux(df_low, 11, "eos")
    flux11eos24, fluxErr11eos24 = getMeanFlux(df_high, 11, "eos")

    vmin = min(df_low[lumiColumn].min(), df_high[lumiColumn].min())
    vmax = max(df_low[lumiColumn].max(), df_high[lumiColumn].max())


    if saveFig:
        plt.figure(figsize=(12, 6))

    chi2ndf = {}
    for i, (df, year, label, M, color) in enumerate(
        zip(
            [df_low, df_high],
            [2023, 2024],
            [
                f"2023 ($< \Phi_{{\mu}} > = ({flux11eos23/1e3:.03f} \pm {fluxErr11eos23/1e3:.03f}) \\times 10^{{3}}$ $nb/cm^{{2}}$)",
                f"2024 ($< \Phi_{{\mu}} > = ({flux11eos24/1e3:.03f} \pm {fluxErr11eos24/1e3:.03f}) \\times 10^{{3}}$ $nb/cm^{{2}}$)"
            ],
            ["s", "o"],
            ["navy", "crimson"]
        )
    ):
        #x = df[f"Flux{tt}_{lumi}"]
        x = df[f"Flux{tt}_{lumi}"]
        y = df['# Coll. IP 1/5'] * np.log(df["Stable B. time [h]"]) / df[f"Flux{tt}_{lumi}"]
        size = df["Stable B. time [h]"] * 40
        color_value = df[lumiColumn]


        # Error bars
        xerr = df[f"FluxErr{tt}_{lumi}"]
        yerr = (df['# Coll. IP 1/5'] * df[f"FluxErr{tt}_{lumi}"] * np.log(df["Stable B. time [h]"])) / (df[f"Flux{tt}_{lumi}"] * df[f"Flux{tt}_{lumi}"])

        slope, intercept, _, _, _ = linregress(x, y)
        x_fit = np.linspace(x.min(), x.max(), 200)
        y_fit = slope * x_fit + intercept
        plt.plot(x_fit, y_fit, color=color, linestyle="--", linewidth=2, label=f"{year} fit")

        # Chi-squared (with error bars considered)
        y_pred = slope * x + intercept
        residuals = (y - y_pred) / yerr  # Normalize residuals by errors
        chi2 = np.sum(residuals**2)
        ndof = len(x) - 2
        chi2_red = chi2 / ndof if ndof > 0 else np.nan
        chi2ndf[year] = chi2_red

        if saveFig:

            plt.errorbar(
                x=x,
                y=y,
                xerr=xerr,
                yerr=yerr,
                fmt='none',
                ecolor=color,
                elinewidth=1.6,
                alpha=0.7,
                capsize=3
            )

            # Scatter points
            sc = plt.scatter(
                x=x,
                y=y,
                c=color_value,
                s=size,
                cmap="jet",
                marker=M,
                alpha=0.55,
                label=label,
                vmin=vmin,
                vmax=vmax
            )




            plt.text(
                0.5, 0.65 - 0.07 * i,
                f"$\\mathbf{{{year} \; fit:}}$ $\\chi^2/ndf$ = {chi2_red:.2f}",
                transform=plt.gca().transAxes,
                fontsize=16, color=color
            )



            plt.text(
                0.6, 0.35,
                f"Stable B. time > {sb_cut} h",
                transform=plt.gca().transAxes,
                fontsize=16, color="k", fontweight=750
            )

    if saveFig:
        # Colorbar
        cbar = plt.colorbar(sc)
        cbar.set_label("$\mathcal{L}_{int}$ ($nb^{-1}$)", fontsize=20)

        # Axes
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel("Muon flux $\Phi_{\mu}$ ($nb/cm^{2}$)", fontsize=20)
        plt.ylabel("$N_{IP1} \; \ln(T_{sb})$", fontsize=20)
        #plt.xlim(0, 0.1)
        #plt.yscale("log")
        #plt.ylim(-6e5, 1e5)
        addMplLabel(mainText="SND@LHC", extraText="2023/2024 PbPb Data", fontsize=18)

        plt.legend(fontsize=12, loc="lower right")
        plt.tight_layout()

        # Save
        plt.savefig(f"/eos/user/i/idioniso/mfout/png/T/{tt}/{I}_AllRuns2-Flux{tt}-{lumi.capitalize()}_T.{sb_cut}.png")
        plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/T/{tt}/{I}_AllRuns2-Flux{tt}-{lumi.capitalize()}_T.{sb_cut}.pdf")


    return chi2ndf[2023], chi2ndf[2024]





def _plot_AllRunsOmega_2(df0, I: int, tt: int, lumi: str = "eos", sb_cut: float = 0, saveFig: bool = True, lumiErr: float = 0.035):
    df = df0.query(f"`Stable B. time [h]` > {sb_cut}")
    df = df.query("Flux11_eos < 50000")
    if lumi.lower() in ["eos"]:
        lumiColumn = "LumiEosOld"
    elif lumi.lower() in ["st", "supertable"]:
        lumiColumn = 'ATLAS Int. Lumi [1/nb]'
    else:
        return None

    df_low = df[df["Run"] < 10000]
    df_high = df[df["Run"] >= 10000]


    flux11eos23, fluxErr11eos23 = getMeanFlux(df_low, 11, "eos")
    flux11eos24, fluxErr11eos24 = getMeanFlux(df_high, 11, "eos")

    vmin = min(df_low[f"Flux{tt}_{lumi}"].min(), df_high[f"Flux{tt}_{lumi}"].min())
    vmax = max(df_low[f"Flux{tt}_{lumi}"].max(), df_high[f"Flux{tt}_{lumi}"].max())


    if saveFig:
        plt.figure(figsize=(12, 6))

    chi2ndf = {}
    for i, (df, year, label, M, color) in enumerate(
        zip(
            [df_low, df_high],
            [2023, 2024],
            [
                f"2023 ($< \Phi_{{\mu}} > = ({flux11eos23/1e3:.03f} \pm {fluxErr11eos23/1e3:.03f}) \\times 10^{{3}}$ $nb/cm^{{2}}$)",
                f"2024 ($< \Phi_{{\mu}} > = ({flux11eos24/1e3:.03f} \pm {fluxErr11eos24/1e3:.03f}) \\times 10^{{3}}$ $nb/cm^{{2}}$)"
            ],
            ["s", "o"],
            ["navy", "crimson"]
        )
    ):
        x = np.log(df["Stable B. time [h]"]) * df['# Coll. IP 1/5']
        y = df[lumiColumn]
        size = df["Stable B. time [h]"] * 40
        color_value = df[f"Flux{tt}_{lumi}"]


        # Error bars
        yerr = lumiErr * df[lumiColumn]

        slope, intercept, _, _, _ = linregress(x, y)
        x_fit = np.linspace(x.min(), x.max(), 200)
        y_fit = slope * x_fit + intercept
        plt.plot(x_fit, y_fit, color=color, linestyle="--", linewidth=2, label=f"{year} fit")

        # Chi-squared (with error bars considered)
        y_pred = slope * x + intercept
        residuals = (y - y_pred) / yerr  # Normalize residuals by errors
        chi2 = np.sum(residuals**2)
        ndof = len(x) - 2
        chi2_red = chi2 / ndof if ndof > 0 else np.nan
        chi2ndf[year] = chi2_red

        if saveFig:

            plt.errorbar(
                x=x,
                y=y,
                yerr=yerr,
                fmt='none',
                ecolor=color,
                elinewidth=1.6,
                alpha=0.7,
                capsize=3
            )

            # Scatter points
            sc = plt.scatter(
                x=x,
                y=y,
                c=color_value,
                s=size,
                cmap="jet",
                marker=M,
                alpha=0.55,
                label=label,
                vmin=vmin,
                vmax=vmax
            )




            plt.text(
                0.05, 0.6 - 0.07 * i,
                f"$\\mathbf{{{year} \; fit:}}$ $\\chi^2/ndf$ = {chi2_red:.2f}",
                transform=plt.gca().transAxes,
                fontsize=16, color=color
            )



            plt.text(
                0.05, 0.35,
                f"Stable B. time > {sb_cut} h",
                transform=plt.gca().transAxes,
                fontsize=16, color="k", fontweight=750
            )

    if saveFig:
        # Colorbar
        cbar = plt.colorbar(sc)
        cbar.set_label("Muon flux $\Phi_{\mu}$ ($nb/cm^{2}$)", fontsize=20)

        # Axes
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.ylabel("$\mathcal{L}_{int}$ ($nb^{-1}$)", fontsize=20)
        plt.xlabel("$N_{IP1} \; \ln(T_{sb})$", fontsize=20)
        #plt.xlim(-1800, 2600)
        plt.ylim(-0.01, 0.1)
        addMplLabel(mainText="SND@LHC", extraText="2023/2024 PbPb Data", fontsize=18)

        plt.legend(fontsize=12, loc="upper left")
        plt.tight_layout()

        # Save
        plt.savefig(f"/eos/user/i/idioniso/mfout/png/T/{tt}/{I}_AllRuns2-Flux{tt}-{lumi.capitalize()}_T.{sb_cut}.png")
        plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/T/{tt}/{I}_AllRuns2-Flux{tt}-{lumi.capitalize()}_T.{sb_cut}.pdf")


    return chi2ndf[2023], chi2ndf[2024]



def _plot_AllRunsOmega_3(df0, I: int, tt: int, lumi: str = "eos", sb_cut: float = 0, saveFig: bool = True, lumiErr: float = 0.035):
    df = df0.query(f"`Stable B. time [h]` > {sb_cut}")
    df = df.query("Flux11_eos < 50000")
    if lumi.lower() in ["eos"]:
        lumiColumn = "LumiEosOld"
    elif lumi.lower() in ["st", "supertable"]:
        lumiColumn = 'ATLAS Int. Lumi [1/nb]'
    else:
        return None

    df_low = df[df["Run"] < 10000]
    df_high = df[df["Run"] >= 10000]


    flux11eos23, fluxErr11eos23 = getMeanFlux(df_low, 11, "eos")
    flux11eos24, fluxErr11eos24 = getMeanFlux(df_high, 11, "eos")

    vmin = min(df_low[f"Flux{tt}_{lumi}"].min(), df_high[f"Flux{tt}_{lumi}"].min())
    vmax = max(df_low[f"Flux{tt}_{lumi}"].max(), df_high[f"Flux{tt}_{lumi}"].max())


    if saveFig:
        plt.figure(figsize=(12, 6))

    chi2ndf = {}
    for i, (df, year, label, M, color) in enumerate(
        zip(
            [df_low, df_high],
            [2023, 2024],
            [
                f"2023 ($< \Phi_{{\mu}} > = ({flux11eos23/1e3:.03f} \pm {fluxErr11eos23/1e3:.03f}) \\times 10^{{3}}$ $nb/cm^{{2}}$)",
                f"2024 ($< \Phi_{{\mu}} > = ({flux11eos24/1e3:.03f} \pm {fluxErr11eos24/1e3:.03f}) \\times 10^{{3}}$ $nb/cm^{{2}}$)"
            ],
            ["s", "o"],
            ["navy", "crimson"]
        )
    ):
        x = df[lumiColumn]
        y = df['# Coll. IP 1/5'] * np.log(df["Stable B. time [h]"]) / df[f"Flux{tt}_{lumi}"]
        size = df["Stable B. time [h]"] * 40
        color_value = df[f"Flux{tt}_{lumi}"]


        # Error bars
        xerr = lumiErr * df[lumiColumn]
        yerr = abs(
            (df['# Coll. IP 1/5'] * np.log(df["Stable B. time [h]"]) * df[f"FluxErr{tt}_{lumi}"]) / df[f"Flux{tt}_{lumi}"]**2
        )

        slope, intercept, _, _, _ = linregress(x, y)
        x_fit = np.linspace(x.min(), x.max(), 200)
        y_fit = slope * x_fit + intercept
        plt.plot(x_fit, y_fit, color=color, linestyle="--", linewidth=2, label=f"{year} fit")

        # Chi-squared (with error bars considered)
        y_pred = slope * x + intercept
        residuals = (y - y_pred) / yerr  # Normalize residuals by errors
        chi2 = np.sum(residuals**2)
        ndof = len(x) - 2
        chi2_red = chi2 / ndof if ndof > 0 else np.nan
        chi2ndf[year] = chi2_red

        if saveFig:

            plt.errorbar(
                x=x,
                y=y,
                xerr=xerr,
                yerr=yerr,
                fmt='none',
                ecolor=color,
                elinewidth=1.6,
                alpha=0.7,
                capsize=3
            )

            # Scatter points
            sc = plt.scatter(
                x=x,
                y=y,
                c=color_value,
                s=size,
                cmap="jet",
                marker=M,
                alpha=0.55,
                label=label,
                vmin=vmin,
                vmax=vmax
            )




            plt.text(
                0.05, 0.6 - 0.07 * i,
                f"$\\mathbf{{{year} \; fit:}}$ $\\chi^2/ndf$ = {chi2_red:.2f}",
                transform=plt.gca().transAxes,
                fontsize=16, color=color
            )



            plt.text(
                0.05, 0.35,
                f"Stable B. time > {sb_cut} h",
                transform=plt.gca().transAxes,
                fontsize=16, color="k", fontweight=750
            )

    if saveFig:
        # Colorbar
        cbar = plt.colorbar(sc)
        cbar.set_label("Muon flux $\Phi_{\mu}$ ($nb/cm^{2}$)", fontsize=20)

        # Axes
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.ylabel("$N_{IP1} \; \ln(T_{sb}) / \Phi_{\mu}$", fontsize=20)
        plt.xlabel("$\mathcal{L}_{int}$ ($nb^{-1}$)", fontsize=20)
        # plt.xlim(-1800, 2600)
        # plt.ylim(-0.01, 0.1)
        addMplLabel(mainText="SND@LHC", extraText="2023/2024 PbPb Data", fontsize=18)

        plt.legend(fontsize=12, loc="upper left")
        plt.tight_layout()

        # Save
        plt.savefig(f"/eos/user/i/idioniso/mfout/png/T/{tt}/{I}_AllRuns2-Flux{tt}-{lumi.capitalize()}_T.{sb_cut}.png")
        plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/T/{tt}/{I}_AllRuns2-Flux{tt}-{lumi.capitalize()}_T.{sb_cut}.pdf")


    return chi2ndf[2023], chi2ndf[2024]



def getLumiArr(run: int):
    return uproot.open(f"/eos/experiment/sndlhc/atlas_lumi/fill_{getFill(run):06d}.root:atlas_lumi").arrays(library="pd")


def getLumiDdf(run, t0):
    fill = getFill(run)
    df = getLumiArr(run)
    df = df.loc[
        df["unix_timestamp"] >= t0
    ]
    delivered_inst_lumi = []
    delivered_unix_timestamp = []
    delivered_run_number = []
    delivered_fill_number = []
    for i in df.index:
        delivered_inst_lumi.append(df.at[i, "var"])
        delivered_unix_timestamp.append(df.at[i, "unix_timestamp"])
    recorded_mask = np.array(True)
    delivered_inst_lumi = np.array(delivered_inst_lumi)
    delivered_unix_timestamp = np.array(delivered_unix_timestamp)
    delivered_deltas = delivered_unix_timestamp[1:] - delivered_unix_timestamp[:-1]
    delivered_mask = delivered_deltas < 600
    delivered_run = np.logical_and(delivered_unix_timestamp[1:] > fill, delivered_mask)


    return np.cumsum(
        np.multiply(
            delivered_deltas[delivered_run], delivered_inst_lumi[1:][delivered_run]
        )
    )[-1]/1e3



def _plot_AllRunsOmega_4(df0,
    I: int, tt: int, lumi: str = "eos", sb_cut: float = 0,
    saveFig: bool = True,
    lumiErr: float = 0.035,
    legendPos: str = "lower right",
    cutTextPos: tuple = (0.05, 0.4),
    fitTextPos: tuple = (0.6, 0.55),
    xlim: tuple = (),
    ylim: tuple = (),
    printRuns: dict = {2023: False, 2024: False},
    excludedRuns = [6975, 7144, 10181, 10204, 10208],
    plot2023: bool = True,
    plot2024: bool = True
):
    zeta = {}

    if lumi.lower() in ["eos"]:
        lumiColumn = "LumiEosOld"
    elif lumi.lower() in ["st", "supertable"]:
        lumiColumn = 'ATLAS Int. Lumi [1/nb]'
    else:
        return None


    df = df0.loc[df0["Stable B. time [h]"] > sb_cut]
    for R in excludedRuns:
        df = df.loc[df["Run"] != R]
    #df = df.query("Flux11_eos < 50000")


    epsH1 = df['Emittance B1 H [um]']
    epsV1 = df['Emittance B1 V [um]']
    epsH2 = df['Emittance B2 H [um]']
    epsV2 = df['Emittance B2 V [um]']
    Nb1   = df['Intensity B1 [1e11 ppb]']
    Nb2   = df['Intensity B2 [1e11 ppb]']
    nb    = df['# Coll. IP 1/5']
    f     = np.ones(len(df)) * 1e9/50
    beta  = df['beta* IP 1 [m]']
    Tsb   = df["Stable B. time [h]"]
    df["LumiDdf"]  = df[lumiColumn]
    L_peak = df['ATLAS Peak Lumi [Hz/ub]']





    #lnT   = np.log(df["Stable B. time [h]"])
    zeta = (f*Nb1*Nb2*nb) / (4 * np.pi * beta * np.sqrt((epsH1 + epsH2)*(epsV1 + epsV2)))
    x = zeta

    for el in (np.log(Tsb) * zeta - x):
        print(el)

    df["Omega"] = x
    # df = df.query("Omega < 2e8")



    df_low = df[df["Run"] < 10000]
    df_high = df[df["Run"] >= 10000]


    flux11eos23, fluxErr11eos23 = getMeanFlux(df_low, 11, "eos")
    flux11eos24, fluxErr11eos24 = getMeanFlux(df_high, 11, "eos")

    vmin = min(df_low[f"Flux{tt}_{lumi}"].min(), df_high[f"Flux{tt}_{lumi}"].min())
    vmax = max(df_low[f"Flux{tt}_{lumi}"].max(), df_high[f"Flux{tt}_{lumi}"].max())


    if saveFig:
        plt.figure(figsize=(12, 6))

    chi2ndf = {}
    for i, (df, year, label, M, color) in enumerate(
        zip(
            [df_low, df_high],
            [2023, 2024],
            [
                f"2023 ($< \Phi_{{\mu}} > = ({flux11eos23/1e3:.03f} \pm {fluxErr11eos23/1e3:.03f}) \\times 10^{{3}}$ $nb/cm^{{2}}$)",
                f"2024 ($< \Phi_{{\mu}} > = ({flux11eos24/1e3:.03f} \pm {fluxErr11eos24/1e3:.03f}) \\times 10^{{3}}$ $nb/cm^{{2}}$)"
            ],
            ["s", "o"],
            ["navy", "crimson"]
        )
    ):
        if plot2023==False and year==2023:
            continue
        if plot2024==False and year==2024:
            continue

        epsH1 = df['Emittance B1 H [um]']
        epsV1 = df['Emittance B1 V [um]']
        epsH2 = df['Emittance B2 H [um]']
        epsV2 = df['Emittance B2 V [um]']
        Nb1   = df['Intensity B1 [1e11 ppb]']
        Nb2   = df['Intensity B2 [1e11 ppb]']
        nb    = df['# Coll. IP 1/5']
        f     = np.ones(len(df)) * 1e9/50
        beta  = df['beta* IP 1 [m]']
        lnT   = np.log(df["Stable B. time [h]"])
        runs  = df["Run"]

        x = df["Omega"]

        y = df['ATLAS Peak Lumi [Hz/ub]']
        size = df["Stable B. time [h]"] * 40
        color_value = df[f"Flux{tt}_{lumi}"]


        # Error bars
        yerr = lumiErr * y

        slope, intercept, _, _, _ = linregress(x, y)
        x_fit = np.linspace(x.min(), x.max(), 200)
        y_fit = slope * x_fit + intercept
        plt.plot(x_fit, y_fit, color=color, linestyle="--", linewidth=2, label=f"{year} fit")

        # Chi-squared (with error bars considered)
        y_pred = slope * x + intercept
        residuals = (y - y_pred) / yerr  # Normalize residuals by errors
        chi2 = np.sum(residuals**2)
        ndof = len(x) - 2
        chi2_red = chi2 / ndof if ndof > 0 else np.nan
        chi2ndf[year] = chi2_red

        if saveFig:

            plt.errorbar(
                x=x,
                y=y,
                yerr=yerr,
                fmt='none',
                ecolor=color,
                elinewidth=1.6,
                alpha=0.7,
                capsize=3
            )

            # Scatter points
            sc = plt.scatter(
                x=x,
                y=y,
                c=color_value,
                s=size,
                cmap="jet",
                marker=M,
                alpha=0.55,
                label=label,
                vmin=vmin,
                vmax=vmax
            )

            if printRuns[year]:
                for j, run in enumerate(runs):
                    # Manually add an offset to the x and y coordinates
                    x_offset = 0.0001  # Adjust as needed for spacing
                    y_offset = 0.0001  # Adjust as needed for spacing

                    # Add text annotation
                    plt.text(
                        x.iloc[j] + x_offset,  # Adjust x position
                        y.iloc[j] + y_offset,  # Adjust y position
                        str(run),  # Annotate with the run number
                        fontsize=9,  # Adjust the font size
                        color=color,  # Color of the annotation
                        alpha=0.7,  # Transparency of the annotation
                        ha='center',  # Horizontal alignment
                        va='center',  # Vertical alignment
                    )



            plt.text(
                fitTextPos[0], fitTextPos[1] - 0.07 * i,
                f"$\\mathbf{{{year} \; fit:}}$ $\\chi^2/ndf$ = {chi2_red:.2f}",
                transform=plt.gca().transAxes,
                fontsize=16, color=color
            )



            plt.text(
                cutTextPos[0], cutTextPos[1],
                f"Stable B. time > {sb_cut} h",
                transform=plt.gca().transAxes,
                fontsize=16, color="k", fontweight=750
            )

    if saveFig:
        # Colorbar
        cbar = plt.colorbar(sc)
        cbar.set_label("Muon flux $\Phi_{\mu}$ ($nb/cm^{2}$)", fontsize=20)

        # Axes
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.ylabel("$\mathcal{L}_{int}$ ($nb^{-1}$)", fontsize=20)
        plt.xlabel("$\ln(T_{sb}) \\times \\frac{f \; N_{B1} \; N_{B2} \; n_{b}}{4 \pi \\beta^{*} \sqrt{(\\varepsilon_1^H + \\varepsilon_2^H)(\\varepsilon_1^V + \\varepsilon_2^V)}}$", fontsize=20)
        if xlim:
            plt.xlim(xlim[0], xlim[1])
        if ylim:
            plt.ylim(ylim[0], ylim[1])
        addMplLabel(mainText="SND@LHC", extraText="2023/2024 PbPb Data", fontsize=18)

        plt.legend(fontsize=12, loc=legendPos)
        plt.tight_layout()

        # Save
        plt.savefig(f"/eos/user/i/idioniso/mfout/png/T/{tt}/{I}_AllRuns2-Flux{tt}-{lumi.capitalize()}_T.{sb_cut}.png")
        plt.savefig(f"/eos/user/i/idioniso/mfout/pdf/T/{tt}/{I}_AllRuns2-Flux{tt}-{lumi.capitalize()}_T.{sb_cut}.pdf")

    plt.show()
    return chi2ndf[2023], chi2ndf[2024]





def _plot_zetaVsL(df0, year: int = 2023):
    if year == 2023:
        df = df0.loc[df0["Run"] < 10000]
    elif year == 2024:
        df = df0.loc[df0["Run"] >= 10000]
    else:
        raise ValueError("Year isn't 2023 or 2024!")

    x = df['ATLAS Peak Lumi [Hz/ub]'] * np.log(df["Stable B. time [h]"])
    y = df["LumiEosOld"]
    z = df["Stable B. time [h]"]
    ey = 0.035 * y

    plt.figure(figsize=(12, 6))
    plt.errorbar(
        x=x,
        y=y,
        yerr=ey,
        fmt='none',
        ecolor="k",
        elinewidth=1.6,
        alpha=0.9,
        capsize=3
    )
    sc = plt.scatter(x, y, c=z, cmap="jet")


    cbar = plt.colorbar(sc)
    cbar.set_label("Stable B. time (h)", fontsize=20)

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylabel("$\mathcal{L}_{int}$ ($nb^{-1}$)", fontsize=18)
    plt.xlabel("$\mathcal{L}_{peak} \; \ln(T_{sb})$", fontsize=18)
    addMplLabel(mainText="SND@LHC", extraText="2023 PbPb Data", fontsize=18)

    plt.legend(fontsize=12)
    plt.tight_layout()


    plt.show()





def plot_NeventsVsTsb(df, year = 2023, Tsb_cut: float = 0, N_cut: int = 0):
    if year==2023:
        df0 = df.loc[df["Run"] <  10000]
    elif year==2024:
        df0 = df.loc[df["Run"] >= 10000]
    elif year not in ["both", "all"]:
        raise ValueError("Year isn't 2023 or 2024!")
    else:
        df0 = df.copy()

    df1 = df0.loc[
        (df0["Stable B. time [h]"] > Tsb_cut) &
        (df0["nEvents"] > N_cut)
    ]
    Flux0 = df0["Flux11_eos"].mean()
    FluxErr0 = df0["Flux11_eos"].std() + df0["FluxErr11_eos"].mean()

    Flux1 = df1["Flux11_eos"].mean()
    FluxErr1 = df1["Flux11_eos"].std() + df0["FluxErr11_eos"].mean()

    plt.figure(figsize=(8, 7))

    if year==2022 or year==2024:
        plt.scatter(y=df0["Stable B. time [h]"], x=df0["nEvents"], c=df0["Flux11_eos"], cmap="jet", marker="s", s=100, label=f"{year}")

        plt.colorbar(label="Muon flux ($nb/cm^{2}$)")

        plt.axvspan(0, 1e7, color='gray', alpha=0.15)
        plt.axhspan(0, 1.75, color='gray', alpha=0.15)

        plt.text(
            0.05, 0.85,
            rf"$\mathrm{{No\:cut:\;\quad\langle \Phi_\mu \rangle = {Flux0/1e3:.03f} \pm {FluxErr0/1e3:.03f}}}$",
            transform=plt.gca().transAxes,
            fontsize=14,
            family='monospace'
        )

        plt.text(
            0.05, 0.8,
            rf"$\mathrm{{With\:cut:\;\langle \Phi_\mu \rangle = {Flux1/1e3:.03f} \pm {FluxErr1/1e3:.03f}}}$",
            transform=plt.gca().transAxes,
            fontsize=14,
            family='monospace'
        )


    else:
        df0_2023 = df0.loc[df0["Run"] < 10000]
        df0_2024 = df0.loc[df0["Run"] >= 10000]

        df1_2023 = df0_2023.loc[
            (df0_2023["Stable B. time [h]"] > Tsb_cut) &
            (df0_2023["nEvents"] > N_cut)
        ]
        df1_2024 = df0_2024.loc[
            (df0_2024["Stable B. time [h]"] > Tsb_cut) &
            (df0_2024["nEvents"] > N_cut)
        ]

        vmin = df0["Flux11_eos"].min() * 0.95
        vmax = df0["Flux11_eos"].max() * 1.05

        plt.scatter(y=df0_2023["Stable B. time [h]"], x=df0_2023["nEvents"], c=df0_2023["Flux11_eos"], cmap="jet", marker="s", s=100, label=f"2023", vmin=vmin, vmax=vmax)
        plt.scatter(y=df0_2024["Stable B. time [h]"], x=df0_2024["nEvents"], c=df0_2024["Flux11_eos"], cmap="jet", marker="x", s=100, label=f"2024", vmin=vmin, vmax=vmax)

        plt.colorbar(label="Muon flux ($nb/cm^{2}$)")

        plt.axvspan(0, 1e7, color='gray', alpha=0.15)
        plt.axhspan(0, 1.75, color='gray', alpha=0.15)


        Flux0_2023 = df0_2023["Flux11_eos"].mean()
        FluxErr0_2023 = df0_2023["Flux11_eos"].std() + df0_2023["FluxErr11_eos"].mean()

        Flux1_2023 = df1_2023["Flux11_eos"].mean()
        FluxErr1_2023 = df1_2023["Flux11_eos"].std() + df0_2023["FluxErr11_eos"].mean()

        Flux0_2024 = df0_2024["Flux11_eos"].mean()
        FluxErr0_2024 = df0_2024["Flux11_eos"].std() + df0_2024["FluxErr11_eos"].mean()

        Flux1_2024 = df1_2024["Flux11_eos"].mean()
        FluxErr1_2024 = df1_2024["Flux11_eos"].std() + df0_2024["FluxErr11_eos"].mean()

        plt.text(
            0.05, 0.85,
            rf"$\mathrm{{(2023) No\:cut:\;\quad\langle \Phi_\mu \rangle = {Flux0_2023/1e3:.03f} \pm {FluxErr0_2023/1e3:.03f}}}$",
            transform=plt.gca().transAxes,
            fontsize=14,
            family='monospace'
        )

        plt.text(
            0.05, 0.8,
            rf"$\mathrm{{(2023) With\:cut:\;\langle \Phi_\mu \rangle = {Flux1_2023/1e3:.03f} \pm {FluxErr1_2023/1e3:.03f}}}$",
            transform=plt.gca().transAxes,
            fontsize=14,
            family='monospace'
        )

        plt.text(
            0.05, 0.75,
            rf"$\mathrm{{(2024) No\:cut:\;\quad\langle \Phi_\mu \rangle = {Flux0_2024/1e3:.03f} \pm {FluxErr0_2024/1e3:.03f}}}$",
            transform=plt.gca().transAxes,
            fontsize=14,
            family='monospace'
        )

        plt.text(
            0.05, 0.7,
            rf"$\mathrm{{(2024) With\:cut:\;\langle \Phi_\mu \rangle = {Flux1_2024/1e3:.03f} \pm {FluxErr1_2024/1e3:.03f}}}$",
            transform=plt.gca().transAxes,
            fontsize=14,
            family='monospace'
        )

    plt.ylabel("Stable B. duration (h)", fontsize=20)

    plt.xlabel("SND@LHC events")
    plt.ylim(0, 11.5)
    plt.xlim(0, 5.2e7)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(prop={'family': 'serif', 'size': 20}, loc="lower right")
    plt.tight_layout()
    plt.show()
