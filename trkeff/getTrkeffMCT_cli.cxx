#include <iostream>
#include <vector>
#include "TString.h"
#include "trkeffMCT.h"

int main(int argc, char** argv) {
    if (argc < 17) {
        std::cerr << "Usage: " << argv[0] << " input_files sigma col_rate L_LHC xmin xmax ymin ymax zRef1 zRef11 zRef3 zRef13 xzMin xzMax yzMin yzMax" << std::endl;
        return 1;
    }

    TString input_files = argv[1];
    double sigma = std::atof(argv[2]);
    double col_rate = std::atof(argv[3]);
    double L_LHC = std::atof(argv[4]);
    double xmin = std::atof(argv[5]);
    double xmax = std::atof(argv[6]);
    double ymin = std::atof(argv[7]);
    double ymax = std::atof(argv[8]);
    double zRef1 = std::atof(argv[9]);
    double zRef11 = std::atof(argv[10]);
    double zRef3 = std::atof(argv[11]);
    double zRef13 = std::atof(argv[12]);
    double xzMin = std::atof(argv[13]);
    double xzMax = std::atof(argv[14]);
    double yzMin = std::atof(argv[15]);
    double yzMax = std::atof(argv[16]);

    auto results = computeTrackingEfficiencies_MCT(
        input_files, sigma, col_rate, L_LHC, xmin, xmax, ymin, ymax,
        zRef1, zRef11, zRef3, zRef13, xzMin, xzMax, yzMin, yzMax
    );

    for (size_t i = 0; i < results.size(); i += 2) {
        std::cout << "Track Type Result: " << results[i] << " +/- " << results[i+1] << std::endl;
    }

    return 0;
}
