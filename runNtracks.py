import os
import pandas

df = pd.read_csv("/eos/user/i/idioniso/tmp/pp23.csv")
runs = df["Run"]
scales = df["scale"]

for run, scale in zip(runs, scales):
    os.system(
        f'python getNtracks.py -i /eos/user/i/idioniso/1_Data/Tracks/run_{run:06d}/\*.root --scale {scale} -o ./trks-test.csv'
    )

