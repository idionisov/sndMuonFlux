#include <climits>
#include <cmath>
#include <map>

#include "TMath.h"
#include "TEfficiency.h"
#include "TClonesArray.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TVector3.h"

#include "MuFilterHit.h"
#include "sndScifiHit.h"
#include "sndRecoTrack.h"
#include "muonFluxUtils.h"
#include "trkeffUtils.h"
#include "histograms.h"

// To easily obtain opposing track types
extern const std::map<int, int> oppositeTrackType = {
    // SF ST <-> DS ST & SF HT <-> DS HT
    {1, 3}, {11, 13}, {3, 1}, {13, 11}
};


bool isNearVetoBar(sndRecoTrack *track, TClonesArray *mfHits, double distance){
    bool result = false;
    for (int i = 0; i < mfHits->GetEntries(); ++i){
        auto mfHit = dynamic_cast<MuFilterHit*>(mfHits->At(i));
        if (
            mfHit->GetSystem()==1 && track->getDoca(mfHit)<=distance
        ) {
            result = true;
        }
    }
    return result;
}



bool isNearUS5Bar(sndRecoTrack *track, TClonesArray *mfHits, double distance){
    bool result = false;
    for (int i = 0; i < mfHits->GetEntries(); ++i){
        auto mfHit = dynamic_cast<MuFilterHit*>(mfHits->At(i));
        if (
            mfHit->GetSystem()==2 &&
            mfHit->GetPlane()==4 &&
            track->getDoca(mfHit)<=distance
        ) {
            result = true;
        }
    }
    return result;
}




std::pair<double, double> computeTrackingEfficiency(
    const TEfficiency& teff,
    double xmin, double xmax,
    double ymin, double ymax
){
    const TH2* hPassed = dynamic_cast<const TH2*>(teff.GetPassedHistogram());
    const TH2* hTotal  = dynamic_cast<const TH2*>(teff.GetTotalHistogram());


    int binx_min = hTotal->GetXaxis()->FindBin(xmin);
    int binx_max = hTotal->GetXaxis()->FindBin(xmax);
    int biny_min = hTotal->GetYaxis()->FindBin(ymin);
    int biny_max = hTotal->GetYaxis()->FindBin(ymax);

    double sum = 0.0;
    double sumsq = 0.0;
    int count = 0;

    for (int ix = binx_min; ix <= binx_max; ix++) {
        for (int iy = biny_min; iy <= biny_max; iy++) {

            int bin = hTotal->GetBin(ix, iy);

            if (!teff.GetTotalHistogram()->GetBinContent(bin))
                continue;

            double eff = teff.GetEfficiency(bin);
            double err = teff.GetEfficiencyErrorLow(bin);

            sum   += eff;
            sumsq += err*err;
            count++;
        }
    }

    if (count == 0)
        return {0.0, 0.0};

    double mean_eff = sum / count;
    double stdev_eff = std::sqrt(sumsq) / count;

    return {mean_eff, stdev_eff};
}



EffResults createAndSaveTEffs(
    int runNum,
    const std::array<int,4>& trackTypes,
    double xmin, double xmax,
    double ymin, double ymax,
    TFile* outFile,
    TTree* effTree
)
{
    // ----------------------------
    //  Create or get output directory
    // ----------------------------
    EffResults results;
    TDirectory* runDir = nullptr;
    if (outFile) {
        outFile->cd();
        runDir = dynamic_cast<TDirectory*>(outFile->Get(Form("Run%d", runNum)));
        if (!runDir)
            runDir = outFile->mkdir(Form("Run%d", runNum));
        }

    for (const auto& [groupName, histPair] : histGroups) {
        TObjArray* histsPassed = histPair.first;
        TObjArray* histsTotal  = histPair.second;

        // Write histograms if output directory exists
        if (runDir) {
            runDir->cd();
            TDirectory* passedDir = runDir->GetDirectory("Passed");
            if (!passedDir) passedDir = runDir->mkdir("Passed");
            passedDir->cd();
            histsPassed->Write();

            TDirectory* totalDir = runDir->GetDirectory("Total");
            if (!totalDir) totalDir = runDir->mkdir("Total");
            totalDir->cd();
            histsTotal->Write();
        }

        int dim = getHGroupDim(groupName);
        for (int i = 0; i < static_cast<int>(trackTypes.size()); ++i) {
            int trackType = trackTypes[i];

            TH1* hPassed1 = dynamic_cast<TH1*>(histsPassed->At(i));
            TH1* hTotal1  = dynamic_cast<TH1*>(histsTotal->At(i));
            if (!hPassed1 || !hTotal1) {
                std::cerr << "Missing hist for group " << groupName << std::endl;
                continue;
            }

            TEfficiency* teff = nullptr;
            TString teffName;
            TString teffTitle = TString::Format("Efficiency %d, run %d", trackType, runNum);

            if (dim == 1) {
                teffName = TString::Format("teff.1D_%s_%d_%d", groupName.c_str(), runNum, trackType);
                TH1D* hP = dynamic_cast<TH1D*>(hPassed1);
                TH1D* hT = dynamic_cast<TH1D*>(hTotal1);
                if (!hP || !hT) {
                    std::cerr << "Expected 1D histograms for " << groupName << std::endl;
                    continue;
                }
                teff = new TEfficiency(*hP, *hT);

            } else if (dim == 2) {
                teffName = TString::Format("teff.2D_%s_%d_%d", groupName.c_str(), runNum, trackType);
                TH2D* hP = dynamic_cast<TH2D*>(hPassed1);
                TH2D* hT = dynamic_cast<TH2D*>(hTotal1);
                if (!hP || !hT) {
                    std::cerr << "Expected 2D histograms for " << groupName << std::endl;
                    continue;
                }
                teff = new TEfficiency(*hP, *hT);
            } else {
                continue;
            }

            teff->SetName(teffName);
            teff->SetTitle(teffTitle);
            teff->SetStatisticOption(runNum == 0 ? TEfficiency::kBBayesian : TEfficiency::kFCP);
            teff->SetConfidenceLevel(0.683);

            if (runDir) {
                runDir->cd();
                teff->Write();
            }

            if (groupName == "x.y") {
                results.at(i) = computeTrackingEfficiency(*teff, xmin, xmax, ymin, ymax);
            }

            delete teff;
        }
    }

    // ----------------------------
    //  Fill TTree if provided
    // ----------------------------
    if (effTree) {
        int branchRunNum;
        int branchTrackType;
        double branchEff;
        double branchErr;

        effTree->Branch("runNum", &branchRunNum, "runNum/I");
        effTree->Branch("trackType", &branchTrackType, "trackType/I");
        effTree->Branch("eff", &branchEff, "eff/D");
        effTree->Branch("effErr", &branchErr, "effErr/D");

        for (size_t i = 0; i < results.size(); ++i) {
            branchRunNum = runNum;
            branchTrackType = trackTypes[i];
            branchEff = results.at(i).first;
            branchErr = results.at(i).second;
            effTree->Fill();
        }
    }

    return results;
}



std::array<int, 5> getNsfHits(TClonesArray* sfHits){
    std::array<int, 5> hitCounts = {0, 0, 0, 0, 0};
    for (int i = 0; i < sfHits->GetEntries(); i++) {
        auto sfHit = static_cast<sndScifiHit*>(sfHits->At(i));
        hitCounts.at(sfHit->GetStation()-1)++;
    }
    return hitCounts;
}

std::array<int, 4> getNdsHits(TClonesArray* mfHits){
    std::array<int, 4> hitCounts = {0, 0, 0, 0};
    for (int i = 0; i < mfHits->GetEntries(); i++) {
        auto dsHit = static_cast<MuFilterHit*>(mfHits->At(i));
        if (dsHit->GetSystem() != 3) continue;
        hitCounts.at(dsHit->GetPlane())++;
    }
    return hitCounts;
}


bool thereIsAMuon(TClonesArray* shipMCTracks) {
    if (!shipMCTracks) return false;

    for (unsigned i = 0; i < shipMCTracks->GetEntries(); ++i) {
        ShipMCTrack* mcTrack = (ShipMCTrack*) shipMCTracks->At(i);
        if (std::abs(mcTrack->GetPdgCode()) == 13) {
            return true;
        }
    }
    return false;
}

bool shipMCTracksCrossedFiducialArea(
    TClonesArray* shipMCTracks,
    double z_ref = 430.0,
    double xmin = -42.0,
    double xmax = -10.0,
    double ymin =  18.0,
    double ymax =  49.0
) {
    if (!shipMCTracks) return false;

    for (unsigned i = 0; i < shipMCTracks->GetEntries(); ++i) {
        ShipMCTrack* mctrack = (ShipMCTrack*) shipMCTracks->At(i);

        if (mctrack->GetMotherId() == -1 && std::abs(mctrack->GetPdgCode()) == 13) {

            TVector3 mcTrkZref = mctrack.GetPointAtZ(z_ref);

            if (
                mcTrkZref.X() >= xmin &&
                mcTrkZref.X() <= xmax &&
                mcTrkZref.Y() >= ymin &&
                mcTrkZref.Y() <= ymax
            ) {
                return true;
            }
        }
    }
    return false;
}


bool sfTrackIsReconstructible(
    TClonesArray* sfHits,
    TClonesArray* sndRecoTracks,
    TClonesArray* scifiPoints
) {
    if (!scifiPoints) return false;

    // Use maps to count hits per station (keys 1-5)
    std::map<int, int> nMCPoints_h;
    std::map<int, int> nMCPoints_v;

    // Initialize maps
    for(unsigned i=1; i<=5; ++i) { nMCPoints_h[i]=0; nMCPoints_v[i]=0; }

    for (int i = 0; i < scifiPoints->GetEntries(); ++i) {
        ScifiPoint* mcPoint = (ScifiPoint*) scifiPoints->At(i);

        if ( !(std::abs(mcPoint->GetPdgCode()) == 13 && mcPoint->GetTrackID() == 0) ) {
            continue;
        }

        int detID = mcPoint->GetDetectorID();

        // Station is the digit at 1,000,000 position
        int station = detID / 1000000;

        // Plane type logic: (detID / 100000) % 2
        // 0 = Horizontal, 1 = Vertical
        int planeType = (detID / 100000) % 2;

        if (planeType == 0) {
            nMCPoints_h[station]++;
        } else if (planeType == 1) {
            nMCPoints_v[station]++;
        }
    }

    int nScifiPoints_h = 0;
    int nScifiPoints_v = 0;

    // Count stations with hits
    for (int sfPlane = 1; sfPlane <= 5; ++sfPlane) {
        if (nMCPoints_v[sfPlane] > 0) nScifiPoints_v++;
        if (nMCPoints_h[sfPlane] > 0) nScifiPoints_h++;
    }

    return (nScifiPoints_h >= 3 && nScifiPoints_v >= 3);
}

bool dsTrackIsReconstructible(
    TClonesArray* mfHits,
    TClonesArray* sndRecoTracks,
    TClonesArray* muFilterPoints
) {
    if (!muFilterPoints) return false;

    std::map<int, int> nMCPoints_h;
    std::map<int, int> nMCPoints_v;

    // Initialize for stations 1 to 4
    for(unsigned i=1; i<=4; ++i) { nMCPoints_h[i]=0; nMCPoints_v[i]=0; }

    for (unsigned i = 0; i < muFilterPoints->GetEntries(); ++i) {
        MuFilterPoint* mcPoint = (MuFilterPoint*) muFilterPoints->At(i);

        if (std::abs(mcPoint->GetPdgCode()) == 13 && mcPoint->GetTrackID() == 0) {
            int detID = mcPoint->GetDetectorID();

            // Check ID range for Downstream
            if (detID < 30000 || detID > 34999) continue;

            // Extract Station and BarNum
            // Using standard string conversion
            std::string detStr = std::to_string(detID);

            // station = int(str(detID)[1])+1
            // detStr[1] is the char at index 1. Subtract '0' to convert char to int.
            // barNum = int(str(detID)[-3:])
            int station = (detStr[1] - '0') + 1;
            int barNum = std::stoi(detStr.substr(detStr.length() - 3));

            if (barNum > 59) {
                nMCPoints_v[station]++;
            } else {
                nMCPoints_h[station]++;
            }
        }
    }

    int nDSPoints_h = 0;
    int nDSPoints_v = 0;

    for (int dsPlane = 1; dsPlane <= 4; ++dsPlane) {
        if (nMCPoints_v[dsPlane] > 0) nDSPoints_v++;
        if (nMCPoints_h[dsPlane] > 0) nDSPoints_h++;
    }

    return (nDSPoints_h >= 3 && nDSPoints_v >= 3);
}
