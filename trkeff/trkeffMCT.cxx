#include <iostream>

#include "TString.h"
#include "TMath.h"
#include "TPython.h"
#include "TChain.h"
#include "TTree.h"
#include "TObject.h"
#include "Rtypes.h"
#include "TROOT.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TDirectory.h"
#include "TEnv.h"

#include "MuFilter.h"
#include "MuFilterHit.h"
#include "sndScifiHit.h"
#include "SNDLHCEventHeader.h"
#include "sndRecoTrack.h"
#include "ShipMCTrack.h"
#include "muonFluxUtils.h"
#include "trkeffTC.h"

std::vector<float> computeTrackingEfficienciesMCT(
    std::string inputStr,
    double sigma    =  8e7,
    double colRate  =  100e6,
    double L_LHC    =  1.0,
    double xmin     = -42.,
    double xmax     = -10.,
    double ymin     =  19.,
    double ymax     =  48.,
    double zRef1    =  430.,
    double zRef11   =  430.,
    double zRef3    =  450.,
    double zRef13   =  450.,
    double xzMin    = -1e12,
    double xzMax    =  1e12,
    double yzMin    = -1e12,
    double yzMax    =  1e12
) {
    std::vector<double> results = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> zRef = {zRef1, zRef11, zRef3, zRef13};

    // 1. Constants & Setup
    double L_MC = col_rate / sigma;
    double Ks = L_LHC / L_MC;

    // Calculate Area (A)
    double A = (xmax-xmin) * (ymax-ymin);

    // Load Chain
    TChain *ch = loadChainWithFallback(inputStr);
    if (!ch) {
        std::cerr << "No valid tree found!" << std::endl;
        return results;
    }
    long nEntries = ch->GetEntries();

    // std::cout << "Entries    : " << nEntries << std::endl;
    // std::cout << "Luminosity : " << L_MC << std::endl;
    // std::cout << "Area       : " << A << std::endl;
    // std::cout << "Scaling    : " << Ks << std::endl;

    SNDEventHeader* eventHeader = nullptr;
    TClonesArray* mcTracks = nullptr;
    TClonesArray* recoTracks = nullptr;
    TClonesArray* scifiPoints = nullptr;
    TClonesArray* muFilterPoints = nullptr;

    ch->SetBranchAddress("EventHeader", &eventHeader);
    ch->SetBranchAddress("MCTrack", &mcTracks);
    ch->SetBranchAddress("Reco_MuonTracks", &recoTracks);
    ch->SetBranchAddress("ScifiPoint", &scifiPoints);
    ch->SetBranchAddress("MuFilterPoint", &muFilterPoints);

    std::map<int, std::vector<double>> w_total_vec;
    std::map<int, std::vector<double>> w_pass_vec;

    for(int tt : trackTypes) {
        w_total_vec[tt] = std::vector<double>();
        w_pass_vec[tt]  = std::vector<double>();
    }

    int next_print = 0;
    for (Long64_t i = 0; i < nEntries; ++i) {
        ch->GetEntry(i);

        // Progress printing
        int progress = (i * 100) / nEntries;
        if (progress >= next_print) {
            std::cout << progress << " %" << std::endl;
            next_print += 5;
        }

        // Global Event Cuts
        // Note: isIP1 is likely a function of SNDEventHeader
        if (!eventHeader->isIP1()) {
            continue;
        }

        bool _sf = sfTrackIsReconstructible(sfHits, recoTracks, scifiPoints);
        bool _ds = DsTrackIsReconstructible(mcTracks);

        if (!_sf && !_ds) continue;

        std::map<int, bool> reco_map;
        reco_map[1]  = _sf;
        reco_map[11] = _sf;
        reco_map[3]  = _ds;
        reco_map[13] = _ds;

        ShipMCTrack* mcTrk0 = (ShipMCTrack*)mcTracks->At(0);
        weight = mcTrk0->GetWeight() * Ks;

        for (int tt : trackTypes) {
            if (!reco_map[tt]) continue;

            w_total_vec[tt].push_back(weight);

            bool eventPassed = false;

            for (int j = 0; j < recoTracks->GetEntriesFast(); ++j) {
                sndRecoTrack* trk = (sndRecoTrack*)recoTracks->At(j);

                if (trk->getTrackType() != tt) continue;

                double xz = trk->getAngleXZ();
                double yz = trk->getAngleYZ();

                if (
                    xz >= xzMin &&
                    xz <= xzMax &&
                    yz >= yzMin &&
                    yz <= yzMax
                ) {
                    eventPassed = true;
                    break;
                }
            }

            // - "weights" list contains ALL weights.
            // - "passed" list contains booleans.
            if (eventPassed) {
                w_pass_vec[tt].push_back(weight);
            }
        }
    }

    std::map<int, std::vector<double>> results;

    // Beta quantile constants
    double low_q = 0.15865;
    double high_q = 0.84135;

    for (int tt : trackTypes) {
        double sum_w_total = 0;
        double sum_w2_total = 0;
        double sum_w_pass = 0;

        for (double w : w_total_vec[tt]) {
            sum_w_total += w;
            sum_w2_total += w * w;
        }
        for (double w : w_pass_vec[tt]) {
            sum_w_pass += w;
        }

        // Calculation: (Sum W)^2 / Sum (W^2)
        double Neff = (sum_w2_total > 0) ? std::pow(sum_w_total, 2) / sum_w2_total : 0.0;

        // Point estimate
        double p_hat = (sum_w_total > 0) ? sum_w_pass / sum_w_total : 0.0;

        double k_eff = p_hat * Neff;

        // Bayesian Beta Parameters
        double a = k_eff + 1.0;
        double b = Neff - k_eff + 1.0;

        // Quantiles using ROOT::Math
        double lo_b = ROOT::Math::beta_quantile(low_q, a, b);
        double hi_b = ROOT::Math::beta_quantile(high_q, a, b);

        double errUp = hi_b - p_hat;
        double errLow = p_hat - lo_b;

        results[tt] = {p_hat, (errUp + errLow) / 2.0};

        std::cout << "TrackType " << tt << ": "
                  << "Eff=" << p_hat << ", "
                  << "Neff=" << Neff << ", "
                  << "Err=" << results[tt].error << std::endl;
    }

    return results;
}
