#include <vector>
#include <array>
#include "TString.h"

std::vector<std::array<int, 4>> getNTracks(
    TString inputStr,
    double xmin, double xmax, double ymin, double ymax,
    double xzMin, double xzMax, double yzMin, double yzMax,
    double zRef1, double zRef11, double zRef3, double zRef13
);
