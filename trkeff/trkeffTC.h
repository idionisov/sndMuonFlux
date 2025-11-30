#include "trkeffUtils.h"

EffResults computeTrackingEfficiencies(const TrkeffConfig& cfg);
std::vector<double> computeTrackingEfficienciesPy(
    TString inputStr,
    TString geoFile,
    TString outFileName,
    TString histParamsFile,
    double xmin, double xmax, double ymin, double ymax,
    double xzMin, double xzMax, double yzMin, double yzMax,
    double zRef1, double zRef11, double zRef3, double zRef13,
    double vetoBarDistance,
    double us5BarDistance,
    double scifiToDSTrackDistance,
    long nBreak
);
