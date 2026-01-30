#include "trkeffUtils.h"

// EffResults computeTrackingEfficiencies(const TrkeffConfig& cfg);
std::vector<float> computeTrackingEfficienciesMCT(
    std::string input_files,
    double sigma,
    double col_rate,
    double L_LHC,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    double zRef1,
    double zRef11,
    double zRef3,
    double zRef13,
    double xzMin,
    double xzMax,
    double yzMin,
    double yzMax
);
