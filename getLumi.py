import pythonHelpers.lumi

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input_files", type=str, help="Input files for luminosity extraction.")
    args = parser.parse_args()

    lumi = pythonHelpers.lumi.get_lumi_eos(args.input_files)
    print(lumi)
