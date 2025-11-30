#pragma once
#include <chrono>
#include <map>
#include <cstddef>
#include "TEfficiency.h"
#include "sndRecoTrack.h"
#include "TH2.h"

extern const std::map<int, int> oppositeTrackType;
using EffPair = std::pair<double,double>;     // (efficiency, error)
using EffResults = std::array<EffPair, 4>;

bool isNearVetoBar(sndRecoTrack *track, TClonesArray *mfHits, double distance = 3.0);
bool isNearUS5Bar(sndRecoTrack *track, TClonesArray *mfHits, double distance = 3.0);


EffPair computeTrackingEfficiency(
    const TEfficiency& teff,
    double xmin = -42.0,
    double xmax = -10.0,
    double ymin = 19.0,
    double ymax = 48.0
);


EffResults createAndSaveTEffs(
    int runNum,
    const std::array<int,4>& trackTypes,
    double xmin, double xmax,
    double ymin, double ymax,
    TFile* outFile = nullptr,
    TTree* effTree = nullptr
);

std::array<int, 5> getNsfHits(TClonesArray* sfHits);
std::array<int, 4> getNdsHits(TClonesArray* mfHits);

struct TrkeffConfig {
    TString inputStr;
    TString geoFile;
    TString outFileName;
    TString histParamsFile;
    double xmin, xmax, ymin, ymax;
    double xzMin, xzMax, yzMin, yzMax;
    std::array<double,4> zRef;
    double vetoBarDistance;
    double us5BarDistance;
    double scifiToDSTrackDistance;
    long nBreak;
};
