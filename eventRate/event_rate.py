import argparse
import ROOT
from sndUtils import SndData


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run', type=int, default=7080)
    args = parser.parse_args()

    data=SndData(Run=args.run)
