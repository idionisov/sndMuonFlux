#include <cmath>
// #include <chrono>

#include "TChain.h"
#include "TString.h"
#include <TError.h>
#include "muonFluxUtils.h"
#include "sndRecoTrack.h"

std::array<int, 4> trackTypes = {1, 11, 3, 13};
std::map<int, int> indexToTrackType = {{0, 1}, {1, 11}, {2, 3}, {3, 13}};
std::map<int, int> trackTypeToIndex = {{1, 0}, {11, 1}, {3, 2}, {13, 3}};


bool hasBranch(
    TTree* tree,
    const TString& branchName
)
{
    if (!tree) {
        std::cerr << "Error: TTree pointer is null!" << std::endl;
        return false;
    }
    return tree->GetBranch(branchName) != nullptr;
}



TChain* loadChainWithFallback(const TString& inputStr)
{
    std::vector<TString> trees;
    trees.push_back("cbmsim");
    trees.push_back("rawConv");

    for (const auto& t : trees) {
        TChain* ch = new TChain(t);
        long nAdded = ch->Add(inputStr);
        if (nAdded == 0) {
            delete ch;
            continue;
        }

        long nEntries = ch->GetEntries();
        if (nEntries > 0) {
            std::cout << "Successfully loaded tree '" << t
                      << "' with " << nEntries << " entries." << std::endl;
            return ch;
        }

        std::cout << "No entries in tree '" << t << "'; trying next." << std::endl;
        delete ch;
    }

    return nullptr;
}


void getSecAsHMS(double seconds, int &h, int &m, int &s)
{
    h = static_cast<int>(seconds / 3600);
    seconds -= h * 3600;
    m = static_cast<int>(seconds / 60);
    seconds -= m * 60;
    s = std::round(seconds);
}


bool trackIsConverged(sndRecoTrack *track)
{
    bool isConverged = track->getTrackFlag();
    bool hasNonZeroMomZ = track->getTrackMom().Z() != 0;

    return isConverged && hasNonZeroMomZ;
}



bool trackIsWithinArea(
    sndRecoTrack *track, double zRef, double xmin, double xmax, double ymin, double ymax
)
{
    TVector3 refPoint = track->getPointAtZ(zRef);
    return refPoint.X() >= xmin && refPoint.X() <= xmax &&
           refPoint.Y() >= ymin && refPoint.Y() <= ymax;
}

bool trackIsWithinAngleRange(
    sndRecoTrack *track, double xzMin, double xzMax, double yzMmin, double yzMax
)
{
    double xz = track->getAngleXZ()*1e3; // mrad
    double yz = track->getAngleYZ()*1e3; // mrad
    return xz >= xzMin && xz <= xzMax &&
           yz >= yzMmin && yz <= yzMax;
}


// void printStatusWithTime(
//     unsigned int i, unsigned int iMax,
//     const std::chrono::steady_clock::time_point &start_time
// )
// {
//     if (iMax == 0) throw std::invalid_argument("iMax cannot be zero!");
//     if (i > iMax) throw std::invalid_argument("i exceeds iMax!");

//     double elapsed = std::chrono::duration<double>(
//         std::chrono::steady_clock::now() - start_time
//     ).count();
//     double percent = (static_cast<double>(i) / iMax) * 100.0;

//     int h, m, s;
//     getSecAsHMS(elapsed, h, m, s);

//     std::cout << "\r >> [" << percent << "%] "
//             << i+1 << "/" << iMax
//             << " " << h << ":" << m << ":" << s
//             << std::flush;
//     if (i+1 == iMax) std::cout << std::endl;
// }

// size_t printStatus(
//     unsigned int i, unsigned int iMax,
//     const std::chrono::steady_clock::time_point &start_time,
//     size_t count
// )
// {
//     double elapsed = std::chrono::duration<double>(
//         std::chrono::steady_clock::now() - start_time
//     ).count();
//     size_t elapsed_int = static_cast<size_t>(std::floor(elapsed));

//     if (elapsed_int > count) {
//         printStatusWithTime(i, iMax, start_time);
//         count += 1;
//     } else if (i+1 == iMax) {
//         printStatusWithTime(i, iMax, start_time);
//     }
//     return count;
// }
