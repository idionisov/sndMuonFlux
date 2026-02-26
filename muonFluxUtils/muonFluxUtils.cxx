#include <cmath>
// #include <chrono>

#include "TChain.h"
#include "TString.h"
#include "SNDLHCEventHeader.h"
#include "ShipMCTrack.h"
#include "TClonesArray.h"
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



//TChain* loadChainWithFallback(const TString& inputStr)
//{
//    std::vector<TString> trees;
//    trees.push_back("cbmsim");
//    trees.push_back("rawConv");
//
//    for (const auto& t : trees) {
//        TChain* ch = new TChain(t);
//        long nAdded = ch->Add(inputStr);
//        if (nAdded == 0) {
//            delete ch;
//            continue;
//        }
//
//        long nEntries = ch->GetEntries();
//        if (nEntries > 0) {
//            std::cout << "Successfully loaded tree '" << t
//                      << "' with " << nEntries << " entries." << std::endl;
//            return ch;
//        }
//
//        std::cout << "No entries in tree '" << t << "'; trying next." << std::endl;
//        delete ch;
//    }
//
//    return nullptr;
//}
TChain* loadChainWithFallback(const TString& inputStr)
{
    TChain dummyChain("dummy");
    long nFiles = dummyChain.Add(inputStr);

    if (nFiles == 0) return nullptr;

    const char* firstFileName = dummyChain.GetListOfFiles()->At(0)->GetTitle();
    TFile* f = TFile::Open(firstFileName);
    if (!f || f->IsZombie()) {
        std::cout << "Error: Could not open representative file " << firstFileName << std::endl;
        if (f) delete f;
        return nullptr;
    }

    TString foundTree = "";
    std::vector<TString> trees = {"cbmsim", "rawConv"};

    for (const auto& t : trees) {
        if (f->GetListOfKeys()->FindObject(t)) {
            foundTree = t;
            break;
        }
    }
    delete f;

    if (foundTree.Length() > 0) {
        TChain* ch = new TChain(foundTree);
        ch->Add(inputStr);

        std::cout << "Successfully loaded tree '" << foundTree << "'" << std::endl;
        return ch;
    }

    std::cout << "No valid tree found in files." << std::endl;
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



MCRateResult computeMCRates(
    TString inputStr,
    double xmin, double xmax,
    double ymin, double ymax,
    double zRef1, double zRef11, double zRef3, double zRef13
) {
    MCRateResult res;
    std::vector<int> trackTypes = {1, 11, 3, 13, 0}; // 0 will represent "mcTrk"
    for (int tt : trackTypes) {
        res.nRate[tt] = 0.0;
        res.nRateErr2[tt] = 0.0;
    }

    TChain* ch = loadChainWithFallback(inputStr);
    if (!ch) return res;

    SNDLHCEventHeader* header = new SNDLHCEventHeader();
    TClonesArray* mcTracks = new TClonesArray("ShipMCTrack");
    TClonesArray* recoTracks = new TClonesArray("sndRecoTrack");

    ch->SetBranchAddress("EventHeader.", &header);
    ch->SetBranchAddress("MCTrack", &mcTracks);
    ch->SetBranchAddress("Reco_MuonTracks", &recoTracks);

    std::map<int, double> zRefs = {{1, zRef1}, {11, zRef11}, {3, zRef3}, {13, zRef13}};

    long nEntries = ch->GetEntries();
    for (long i = 0; i < nEntries; ++i) {
        ch->GetEntry(i);

        if (!header->isIP1()) continue;

        double w = 0.0;
        bool hasPrimaryMuon = false;

        for (int j = 0; j < mcTracks->GetEntriesFast(); ++j) {
            ShipMCTrack* mct = (ShipMCTrack*)mcTracks->At(j);
            if (mct->GetMotherId() == -1 && std::abs(mct->GetPdgCode()) == 13) {
                w = mct->GetWeight();
                hasPrimaryMuon = true;

                // mcTrk calculation (fiducial at z=430)
                double pz = mct->GetPz();
                if (pz != 0) {
                    double t = (430.0 - mct->GetStartZ()) / pz;
                    double x = mct->GetStartX() + mct->GetPx() * t;
                    double y = mct->GetStartY() + mct->GetPy() * t;

                    if (x >= xmin && x <= xmax && y >= ymin && y <= ymax) {
                        res.nRate[0] += w;
                        res.nRateErr2[0] += w * w;
                    }
                }
                break;
            }
        }

        if (!hasPrimaryMuon) continue;

        std::set<int> seenTypes;
        for (int j = 0; j < recoTracks->GetEntriesFast(); ++j) {
            sndRecoTrack* trk = (sndRecoTrack*)recoTracks->At(j);
            if (!trk->getTrackFlag() || trk->getTrackMom().Z() == 0) continue;

            int tt = trk->getTrackType();
            if (zRefs.find(tt) == zRefs.end()) continue;

            TVector3 pos = trk->getPointAtZ(zRefs[tt]);
            if (pos.X() >= xmin && pos.X() <= xmax && pos.Y() >= ymin && pos.Y() <= ymax) {
                seenTypes.insert(tt);
            }
        }

        for (int tt : seenTypes) {
            res.nRate[tt] += w;
            res.nRateErr2[tt] += w * w;
        }
    }

    // delete ch;
    // delete eh;
    // delete mcTracks;
    // delete recoTracks;

    return res;
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
