#ifndef TRKEFFMCT_H
#define TRKEFFMCT_H

#include "TString.h"
#include <vector>

std::vector<double> computeTrackingEfficiencies_MCT(
    TString input_files,
    double sigma    =  8e7,
    double col_rate =  100e6,
    double L_LHC    =  1.0,
    double xmin     = -42.,
    double xmax     = -10.,
    double ymin     =  19.,
    double ymax     =  48.,
    double zRef1    =  430.,
    double zRef11   =  430.,
    double zRef3    =  450.,
    double zRef13   =  450.,
    double xzMin    = -1e12,
    double xzMax    =  1e12,
    double yzMin    = -1e12,
    double yzMax    =  1e12
);

#endif
