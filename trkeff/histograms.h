#pragma once

#include <array>
#include <map>
#include <string>
#include <utility>
#include "TObjArray.h"
#include "TEnv.h"

using HistPair = std::pair<TObjArray*, TObjArray*>;
using HistCreator = HistPair(*)(int runNum, const std::array<int,4>& trackTypes, const TEnv& hParams);


extern std::map<std::string, HistPair> histGroups;

HistPair createXYHists(int run, const std::array<int, 4>& trackTypes, const TEnv& hParams);
HistPair createXHists(int run, const std::array<int, 4>& trackTypes, const TEnv& hParams);
HistPair createYHists(int run, const std::array<int, 4>& trackTypes, const TEnv& hParams);
HistPair createNhitsAllHists(int run, const std::array<int, 4>& trackTypes, const TEnv& hParams);
HistPair createNhitsPlaneHists(int run, const std::array<int, 4>& trackTypes, const TEnv& hParams);
HistPair createXvsNhitsAllHists(int run, const std::array<int, 4>& trackTypes, const TEnv& hParams);
HistPair createXvsNhitsPlaneHists(int run, const std::array<int, 4>& trackTypes, const TEnv& hParams);
HistPair createYvsNhitsAllHists(int run, const std::array<int, 4>& trackTypes, const TEnv& hParams);
HistPair createYvsNhitsPlaneHists(int run, const std::array<int, 4>& trackTypes, const TEnv& hParams);
HistPair createChi2ndfHists(int run, const std::array<int, 4>& trackTypes, const TEnv& hParams);

extern std::map<std::string, HistCreator> histFactory;

int getHGroupDim(const std::string& s);
TH1D* getHist1D(const std::string& key, int idx, bool passed);
TH2D* getHist2D(const std::string& key, int idx, bool passed);
