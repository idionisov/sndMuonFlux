import uproot
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def getFormattedData(
    inputFile: str = "/eos/user/i/idioniso/mfout/muonFlux.csv"
):
    dfMuonFlux = pd.read_csv(f"/eos/user/i/idioniso/mfout/muonFlux.csv", index_col="run")
    runs = dfMuonFlux.index.tolist()

    trackTypes = {
        1:  "Scifi simple\ntracking",
        11: "Scifi Hough\ntransform",
        3:  "DS simple\ntracking",
        13: "DS Hough\ntransform"
    }

    muonFlux = []
    for run in runs:
        for tt, TT in trackTypes.items():
            if run==-1:
                muonFlux.append((TT, f"Monte Carlo (PbPb)", dfMuonFlux.at[run, f"flux_{tt}"], dfMuonFlux.at[run, f"fluxErr_{tt}"]))
            else:
                muonFlux.append((TT, f"Run {run}", dfMuonFlux.at[run, f"flux_{tt}"], dfMuonFlux.at[run, f"fluxErr_{tt}"]))

    columns=['Track type', 'Dataset', 'Muon flux', 'Errors']
    return pd.DataFrame(muonFlux, columns=columns)


def plotMuonFlux(
    inputFile: str = "/eos/user/i/idioniso/mfout/muonFlux.csv",
    style: str = "ddf-pub"
):
    plt.style.use(style)

    columns=['Track type', 'Dataset', 'Muon flux', 'Errors']
    data = getFormattedData(inputFile)
    formats = ("png", "pdf")
    plt.figure(figsize=(10, 6))

    ax1 = sns.barplot(
        data=data,
        x=columns[0],
        hue=columns[1],
        y=columns[2],
        errorbar=None,
        gap=0.2
    )

    categories = data.iloc[:, 0].unique()
    types 	   = data.iloc[:, 1].unique()
    ntypes	   = len(types)
    c1 = -0.00138889*(ntypes**3) + 0.02482143*(ntypes**2) - 0.16128968*ntypes + 0.03285714
    c2 = -0.00296296*(ntypes**3) + 0.05303571*(ntypes**2) - 0.33578704*ntypes + 0.88238095

    for i_category, category in enumerate(categories):
        for i_type, type in enumerate(types):
            y    = data[(data.iloc[:, 0] == category) & (data.iloc[:, 1] == type)].iloc[:, 2].values[0]
            yerr = data[(data.iloc[:, 0] == category) & (data.iloc[:, 1] == type)].iloc[:, 3].values[0]
            ax1.errorbar(
                i_category + c1 + c2*i_type,
                y, yerr=yerr, fmt='none', c='black', capsize=5
            )


    plt.xticks(rotation=0, ha='center')
    plt.ylabel('Muon flux [$nb \; cm^{-2}$]')
    plt.xlabel('')
    plt.grid(axis="y")
    plt.ylim(1e4, 2.05e4)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, +1.135), ncol=5, fontsize=10)
    plt.tight_layout()

    for fmt in formats:
        fname = f"/eos/user/i/idioniso/mfout/muonFluxPlot.{fmt}"
        plt.savefig(fname, format=fmt)
        print(fname)
    plt.show()


def main():
    plotMuonFlux(f"/eos/user/i/idioniso/mfout/muonFlux.csv", style="ddf-pub")

if __name__=="__main__":
    main()
