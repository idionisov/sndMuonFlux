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
    double ymin, double ymax
){
    const TH2* hPassed = dynamic_cast<const TH2*>(teff.GetPassedHistogram());
    const TH2* hTotal  = dynamic_cast<const TH2*>(teff.GetTotalHistogram());

    if (!hPassed || !hTotal) return {0.0, 0.0};

    int binx_min = hTotal->GetXaxis()->FindBin(xmin);
    int binx_max = hTotal->GetXaxis()->FindBin(xmax);
    int biny_min = hTotal->GetYaxis()->FindBin(ymin);
    int biny_max = hTotal->GetYaxis()->FindBin(ymax);

    double total_passed = 0.0;
    double total_total = 0.0;

    // Store bin data to calculate variance later
    struct BinData { double p; double t; };
    std::vector<BinData> bins;

    // 1. Pool data
    for (int ix = binx_min; ix <= binx_max; ix++) {
        for (int iy = biny_min; iy <= biny_max; iy++) {
            int bin = hTotal->GetBin(ix, iy);
            double t = hTotal->GetBinContent(bin);
            double p = hPassed->GetBinContent(bin);

            if (t > 0) {
                total_total += t;
                total_passed += p;
                bins.push_back({p, t});
            }
        }
    }

    if (total_total <= 0.0) return {0.0, 0.0};

    // 2. Calculate Global Efficiency and Global Statistical Error
    double mean_eff = total_passed / total_total;
    // Using standard binomial error for the global stat error
    double global_stat_err_sq = (mean_eff * (1.0 - mean_eff)) / total_total;

    // If there is only 1 bin, there is no spatial variance.
    if (bins.size() <= 1) {
        return {mean_eff, std::sqrt(global_stat_err_sq)};
    }

    // 3. Calculate Spatial (Systematic) Variance
    double weighted_sq_diff_sum = 0.0;
    double expected_stat_var_sum = 0.0;

    for (const auto& b : bins) {
        double eff_i = b.p / b.t;
        double weight = b.t / total_total; // Weight bins by their track count

        // Observed variance of this bin from the mean
        weighted_sq_diff_sum += weight * (eff_i - mean_eff) * (eff_i - mean_eff);

        // Expected statistical variance of this specific bin
        // (If eff_i is 1 or 0, this is 0, which is handled correctly)
        expected_stat_var_sum += weight * (eff_i * (1.0 - eff_i) / b.t);
    }

    // True spatial variance = Observed variance - Expected Statistical noise
    double spatial_syst_err_sq = weighted_sq_diff_sum - expected_stat_var_sum;

    // If spatial_syst_err_sq < 0, it means all observed spread is purely
    // explained by statistical fluctuations. We floor it at 0.
    if (spatial_syst_err_sq < 0.0) {
        spatial_syst_err_sq = 0.0;
    }

    // 4. Combine in Quadrature
    double total_err = std::sqrt(global_stat_err_sq + spatial_syst_err_sq);

    return {mean_eff, total_err};
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
