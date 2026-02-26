#include <iostream>
#include <vector>
#include "TString.h"
#include "trkeffTT.h"
#include "muonFluxUtils.h"

int main(int argc, char** argv)
{
    if (argc < 21) {
        std::cerr << "Too few arguments " << argc << "/20" << std::endl;
        return -1;
    }

    TString inputStr            = argv[1];
    TString geoFile             = argv[2];
    TString outFileName         = argv[3];
    double xmin                = atof(argv[4]);
    double xmax                = atof(argv[5]);
    double ymin                = atof(argv[6]);
    double ymax                = atof(argv[7]);
    double xzMin               = atof(argv[8]);
    double xzMax               = atof(argv[9]);
    double yzMin               = atof(argv[10]);
    double yzMax               = atof(argv[11]);
    double zRef1               = atof(argv[12]);
    double zRef11              = atof(argv[13]);
    double zRef3               = atof(argv[14]);
    double zRef13              = atof(argv[15]);
    double vetoBarDistance     = atof(argv[16]);
    double us5BarDistance      = atof(argv[17]);
    double scifiToDSTrackDistance = atof(argv[18]);
    long nBreak              = atol(argv[19]);
    TString histParamsFile      = argv[20];

    auto effResults = computeTrackingEfficiencies_TT(
        inputStr, geoFile, outFileName, histParamsFile,
        xmin, xmax, ymin, ymax, xzMin, xzMax, yzMin, yzMax,
        zRef1, zRef11, zRef3, zRef13, vetoBarDistance,
        us5BarDistance, scifiToDSTrackDistance, nBreak
    );

    unsigned i = 0;
    while (i < effResults.size()){
        double eff    = effResults.at(i);
        double effErr = effResults.at(i+1);

        std::cout << "Track type " << trackTypes.at(i/2) << ":	"
                  << eff << " +/- " << effErr
                  << std::endl;

        i+=2;
    }

    return 0;
}
