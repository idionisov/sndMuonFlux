import os, csv

xy_eff_range = {
    'min': {'x': -42., 'y': 19.},
    'max': {'x': -10., 'y': 48.}
}

def eosIsMounted() -> bool:
    eos = "/eos/user/i/idioniso"

    if os.path.exists(eos) and os.path.isdir(eos):
        if not os.listdir(eos):
            return False
        else:
            return True
    else:
        return False


def getEosDir(remoteEos: bool) -> str:
    if not remoteEos:
        return "/EOS/user/i/idioniso"
    else:
        if eosIsMounted():
            return "/eos/user/i/idioniso"
        else:
            print("Remote eos is requested, but it is not mounted!")
            print("Falling back to local /EOS/user/i/idioniso")
            return "/EOS/user/i/idioniso"


def getInputDir(input: str, remoteEos: bool) -> str:
    if input:
        return input
    return f"{getEosDir(remoteEos)}/1_Data/Tracks"


def getOutput(output: str, remoteEos: bool) -> str:
    if output:
        return output
    return f"{getEosDir(remoteEos)}/mfout/nTracks.csv"


def saveToCsv(header, data, fout: str):
    foutExists = os.path.isfile(fout)

    with open(fout, mode='a', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)

        if not foutExists or os.stat(fout).st_size == 0:
            writer.writerow(header)

        writer.writerows(data)
    print(f"Output: {fout}")
