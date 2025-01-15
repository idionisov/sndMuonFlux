import ROOT
from sndUtils import SndData, DdfTrack





input_dir = "/eos/user/i/idioniso/1_Data/Tracks/"
data = SndData(Run=7080, InputDir=input_dir)

for i_event, event in enumerate(data.Tree):
    if not event.Reco_MuonTracks:
        continue

    # Create a DdfTrack from the first track
    trk = DdfTrack(event.Reco_MuonTracks[0], event)
    trk.Print()

    break
