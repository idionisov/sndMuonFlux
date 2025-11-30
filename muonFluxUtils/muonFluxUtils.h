#include <iostream>
#include <iostream>
#include <cmath>
#include <stdexcept>
// #include <chrono>
#include <map>
#include <cstddef>

#include "TChain.h"
#include "TString.h"
#include "sndRecoTrack.h"

// To easily obtain track type from track type index and the other way around
extern std::array<int, 4> trackTypes;
extern std::map<int, int> indexToTrackType;
extern std::map<int, int> trackTypeToIndex;

bool hasBranch(TTree* tree, const TString& branchName);
TChain* loadChainWithFallback(const TString& inputStr);
void getSecAsHMS(double seconds, int &h, int &m, int &s);
bool trackIsConverged(sndRecoTrack *track);
bool trackIsWithinArea(sndRecoTrack *track, double zRef, double xmin, double xmax, double ymin, double ymax);
bool trackIsWithinAngleRange(sndRecoTrack *track, double xzMin, double xzMax, double yzMmin, double yzMax);

// void printStatusWithTime(
//     unsigned int i, unsigned int iMax,
//     const std::chrono::steady_clock::time_point &start_time
// );

// size_t printStatus(
//     unsigned int i, unsigned int iMax,
//     const std::chrono::steady_clock::time_point &start_time,
//     size_t count = 0
// );
