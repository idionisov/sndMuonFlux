#include <climits>
#include <cmath>
#include <map>
#include <numeric>
#include <vector>

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
    double ymin, double ymax,
    bool verbose // for debugging
){
    const TH2* hPassed = dynamic_cast<const TH2*>(teff.GetPassedHistogram());
    const TH2* hTotal  = dynamic_cast<const TH2*>(teff.GetTotalHistogram());

    if (!hPassed || !hTotal) return {0.0, 0.0};

    int binx_min = hTotal->GetXaxis()->FindBin(xmin);
    int binx_max = hTotal->GetXaxis()->FindBin(xmax);
    int biny_min = hTotal->GetYaxis()->FindBin(ymin);
    int biny_max = hTotal->GetYaxis()->FindBin(ymax);

    struct BinData { double eff; double w; };
    std::vector<BinData> bins;

    double sum_w = 0.0;     // Total tracks
    double sum_w_eff = 0.0; // Total passed tracks
    double sum_w_sq = 0.0;  // Sum of squared weights (needed for variance correction)

    for (int ix = binx_min; ix <= binx_max; ix++) {
        for (int iy = biny_min; iy <= biny_max; iy++) {
            int bin = hTotal->GetBin(ix, iy); // Global bin index
            double t = hTotal->GetBinContent(bin);
            double p = hPassed->GetBinContent(bin);

            if (t > 0) {
                // Bypass TEfficiency errors and calculate raw bin efficiency
                double eff = p / t; 
                
                // Track-Weighting: Weight is exactly the number of tracks
                double w = t; 

                sum_w += w;
                sum_w_sq += (w * w);
                sum_w_eff += (w * eff);

                bins.push_back({eff, w});
            }
        }
    }

    if (sum_w <= 0.0 || bins.empty()) return {0.0, 0.0};

    // 1. The Track-Weighted Mean (Mathematically identical to Total Passed / Total Tracks)
    double mean_eff = sum_w_eff / sum_w;

    // 2. Global Statistical Error on the Mean (Standard binomial statistics)
    double global_stat_var = (mean_eff * (1.0 - mean_eff)) / sum_w;

    // 3. Track-Weighted Spatial Spread (Population Variance)
    double spatial_var = 0.0;
    if (bins.size() > 1) {
        double weighted_sq_diff_sum = 0.0;
        for (const auto& b : bins) {
            weighted_sq_diff_sum += b.w * std::pow(b.eff - mean_eff, 2);
        }
        
        // Unbiased weighted sample variance formula (reliability weights)
        double correction = sum_w / ((sum_w * sum_w) - sum_w_sq);
        spatial_var = correction * weighted_sq_diff_sum;
    }

    double stat_err_mean = std::sqrt(global_stat_var);
    double syst_err_spatial = std::sqrt(spatial_var);
    
    // Unified error for standard error propagation
    double total_unified_variance = global_stat_var + spatial_var;
    double total_unified_err = std::sqrt(total_unified_variance);

    if (verbose) {
        std::cout << "[Eff]: " << mean_eff << std::endl;
        std::cout << "  -> Total combined error: +/- " << total_unified_err << std::endl;
        std::cout << "  -> Error on the mean: +/- " << stat_err_mean << " (stat)" << std::endl;
        std::cout << "  -> Spatial spread across xy: +/- " << syst_err_spatial << " (syst)" << std::endl;
        std::cout << "[Bins used]: " << bins.size() << std::endl;
        std::cout << "[Total tracks in ROI]: " << sum_w << std::endl;
    }

    return {mean_eff, total_unified_err}; 
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
                results.at(i) = computeTrackingEfficiency(*teff, xmin, xmax, ymin, ymax, true);
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

        if (!effTree->GetBranch("runNum")) effTree->Branch("runNum", &branchRunNum, "runNum/I");
        if (!effTree->GetBranch("trackType")) effTree->Branch("trackType", &branchTrackType, "trackType/I");
        if (!effTree->GetBranch("eff")) effTree->Branch("eff", &branchEff, "eff/D");
        if (!effTree->GetBranch("effErr")) effTree->Branch("effErr", &branchErr, "effErr/D");

        // Set addresses for existing branches to ensure Fill() uses local variables
        effTree->SetBranchAddress("runNum", &branchRunNum);
        effTree->SetBranchAddress("trackType", &branchTrackType);
        effTree->SetBranchAddress("eff", &branchEff);
        effTree->SetBranchAddress("effErr", &branchErr);

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
