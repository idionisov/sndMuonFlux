#include <cmath>
#include <iostream>
#include <chrono>

#include "TObject.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"
#include "TClonesArray.h"
#include "TString.h"
#include "SNDLHCEventHeader.h"
#include "sndRecoTrack.h"

#include "muonFluxUtils.h"
// #include "utilsNTracks.h"


std::vector<std::array<int, 4>> getNTracks(
    TString inputStr,
    double xmin = -42.0,
    double xmax = -10.0,
    double ymin =  19.0,
    double ymax =  48.0,
    double xzMin = -1e12,
    double xzMax =  1e12,
    double yzMin = -1e12,
    double yzMax =  1e12,
    double zRef1  = 430.0,
    double zRef11 = 430.0,
    double zRef3  = 450.0,
    double zRef13 = 450.0
) {
    std::array<double, 4> zRef = {zRef1, zRef11, zRef3, zRef13};
    std::vector<std::array<int, 4>> result;

    // track types {sf st, sf ht, ds st, ds ht} or equivalently {1, 11, 3, 13}
    std::array<int, 4> trackCounts = {0, 0, 0, 0};

    // For bunch structure corrections
    std::array<int, 4> trackCountsIP2     = {0, 0, 0, 0}; // Solely for IP2
    std::array<int, 4> trackCountsB1Only  = {0, 0, 0, 0}; // B1Only
    std::array<int, 4> trackCountsB2noB1  = {0, 0, 0, 0}; // B2noB1
    std::array<int, 4> trackCountsIP2B1B2 = {0, 0, 0, 0}; // IP2 and B1 and B2 but not IP1

    // Loading input Files
    TChain *ch = loadChainWithFallback(inputStr);
    if (!ch) {
        std::cerr << "No valid tree found in input files. Exiting." << std::endl;
        return result;
    }
    long nEntries = ch->GetEntries();


    // Linking event header and tracks
    TClonesArray * tracks = new TClonesArray("sndRecoTrack");
    ch->SetBranchAddress("Reco_MuonTracks", &tracks);

    SNDLHCEventHeader * eventHeader = new SNDLHCEventHeader();
    ch->SetBranchAddress("EventHeader.", &eventHeader);

    ch->GetEntry(0);
    int run = eventHeader->GetRunId();


    // Event loop
    std::cout << "Starting event loop..." << std::endl;
    auto startTime = std::chrono::steady_clock::now();
    int lastStatusPercentage = -1;
    for (unsigned i_entry = 0; i_entry < nEntries; ++i_entry) {

        // Print status every 5%
        int currentStatusPercentage = static_cast<int>(i_entry*100.0 / nEntries);
        if (currentStatusPercentage % 5 == 0 && currentStatusPercentage != lastStatusPercentage){
            std::cout << currentStatusPercentage << " %" << std::endl;
            lastStatusPercentage = currentStatusPercentage;
        }

        ch->GetEntry(i_entry);
        if (!eventHeader->isIP1()) continue;

        // Process tracks
        for (int i = 0; i < tracks->GetEntries(); ++i) {
            sndRecoTrack * track = (sndRecoTrack *) tracks->At(i);
            int tt = track->getTrackType();
            int i_tt = trackTypeToIndex.at(tt);

            if (!trackIsConverged(track)) continue;
            if (!trackIsWithinArea(track, zRef.at(i_tt), xmin, xmax, ymin, ymax)) continue;


            trackCounts[i_tt]++;
            if (
                eventHeader->isIP2() && !eventHeader->isB1() && !eventHeader->isB2()
            ) { trackCountsIP2[i_tt]++; }
            if (eventHeader->isB1Only()) trackCountsB1Only[i_tt]++;
            if (eventHeader->isB2noB1()) trackCountsB2noB1[i_tt]++;
            if (
                eventHeader->isIP2() && eventHeader->isB1() && eventHeader->isB2()
            ) { trackCountsIP2B1B2[i_tt]++; }
        }
    }
    std::cout << "Event loop completed." << std::endl;

    result.push_back(trackCounts);
    result.push_back(trackCountsIP2);
    result.push_back(trackCountsB1Only);
    result.push_back(trackCountsB2noB1);
    result.push_back(trackCountsIP2B1B2);

    // Print elapsed Time
    auto endTime = std::chrono::steady_clock::now();
    double elapsedSec = std::chrono::duration<double>(endTime - startTime).count();
    int h, m, s;
    getSecAsHMS(elapsedSec, h, m, s);
    std::cout << "Elapsed time: "
              << h << "h " << m << "m " << s << "s ("
              << elapsedSec << " seconds)" << std::endl;

    // 0 -> IP1
    // 1 -> IP2
    // 2 -> B1
    // 3 -> B2
    // 4 -> IP2B1B2
    return result;
}


int main(int argc, char **argv)
{
    if (argc < 14) {
        std::cerr << "Too few arguments!" << std::endl;
        exit(-1);
    }
    std::cout << "Starting track counter..." << std::endl;

    TString inputStr(argv[1]);
    double xmin = atof(argv[2]);
    double xmax = atof(argv[3]);
    double ymin = atof(argv[4]);
    double ymax = atof(argv[5]);
    double xzMin = atof(argv[6]);
    double xzMax = atof(argv[7]);
    double yzMin = atof(argv[8]);
    double yzMax = atof(argv[9]);
    double zRef1  = atof(argv[10]);
    double zRef11 = atof(argv[11]);
    double zRef3  = atof(argv[12]);
    double zRef13 = atof(argv[13]);

    // 0 -> IP1
    // 1 -> IP2
    // 2 -> B1
    // 3 -> B2
    // 4 -> IP2B1B2
    std::vector<std::array<int, 4>> allTrackCounts = getNTracks(
        inputStr,
        xmin, xmax, ymin, ymax,
        xzMin, xzMax, yzMin, yzMax,
        zRef1, zRef11, zRef3, zRef13
    );

    // Print results
    for (unsigned i_tt = 0; i_tt < 4; ++i_tt) {
        int tt = indexToTrackType.at(i_tt);

        std::cout << "Track type " << tt << ":\t" << allTrackCounts.at(0).at(i_tt) << "\t"
            << "(IP2: " << allTrackCounts.at(1).at(i_tt)
            << ", B1Only: " << allTrackCounts.at(2).at(i_tt)
            << ", B2noB1: " << allTrackCounts.at(3).at(i_tt)
            << ", IP2B1B2: " << allTrackCounts.at(4).at(i_tt)
            << ")" << std::endl;
    }



    return 0;
}
