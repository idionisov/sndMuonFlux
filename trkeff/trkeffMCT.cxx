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

std::vector<float> get_trkeff_mct(
    std::string input_files,
    double sigma    =  8e7,
    double col_rate =  100e6,
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
    // 1. Constants & Setup
    double L_MC = col_rate / sigma;
    double Ks = L_LHC / L_MC;
    std::vector<int> trackTypes = {1, 11, 3, 13};

    // Calculate Area (A)
    double A = (x_range.second - x_range.first) * (y_range.second - y_range.first);

    // Load Chain
    TChain* ch = new TChain("cbmsim");
    ch->Add(input_files.c_str());
    Long64_t n_entries = ch->GetEntries();

    std::cout << "Entries    : " << n_entries << std::endl;
    std::cout << "Luminosity : " << L_MC << std::endl;
    std::cout << "Area       : " << A << std::endl;
    std::cout << "Scaling    : " << Ks << std::endl;

    // 2. Branch Setup
    SNDEventHeader* eventHeader = nullptr;
    TClonesArray* mcTracks = nullptr;
    TClonesArray* recoTracks = nullptr;

    ch->SetBranchAddress("EventHeader", &eventHeader);
    ch->SetBranchAddress("MCTrack", &mcTracks);
    ch->SetBranchAddress("Reco_MuonTracks", &recoTracks);

    // 3. Storage for weights and pass/fail
    // We store weights separately for passed events and total events to compute N_eff later
    std::map<int, std::vector<double>> w_total_vec;
    std::map<int, std::vector<double>> w_pass_vec;

    // Initialize vectors
    for(int tt : trackTypes) {
        w_total_vec[tt] = std::vector<double>();
        w_pass_vec[tt]  = std::vector<double>(); // Only stores weights if passed is true
    }

    // 4. Event Loop
    int next_print = 0;
    for (Long64_t i = 0; i < n_entries; ++i) {
        ch->GetEntry(i);

        // Progress printing
        int progress = (i * 100) / n_entries;
        if (progress >= next_print) {
            std::cout << progress << " %" << std::endl;
            next_print += 5;
        }

        // Global Event Cuts
        // Note: isIP1 is likely a function of SNDEventHeader
        if (!eventHeader->isIP1() || !ThereIsAMuon(eventHeader, mcTracks)) {
            continue;
        }

        bool _sf = SfTrackIsReconstructible(mcTracks);
        bool _ds = DsTrackIsReconstructible(mcTracks);

        if (!_sf && !_ds) continue;

        // Map track type to reconstructibility boolean
        std::map<int, bool> reco_map;
        reco_map[1] = _sf;
        reco_map[11] = _sf;
        reco_map[3] = _ds;
        reco_map[13] = _ds;

        // Get MC Weight
        double weight = 1.0;
        if (mcTracks->GetEntriesFast() > 0) {
            ShipMCTrack* mc = (ShipMCTrack*)mcTracks->At(0);
            weight = mc->GetWeight() * Ks;
        }

        // Loop over track types
        for (int tt : trackTypes) {
            if (!reco_map[tt]) continue;

            // Store weight in "total" denominator for this track type
            w_total_vec[tt].push_back(weight);

            bool eventPassed = false;

            // Check Reco Tracks
            for (int j = 0; j < recoTracks->GetEntriesFast(); ++j) {
                sndRecoTrack* trk = (sndRecoTrack*)recoTracks->At(j);

                if (trk->getTrackType() != tt) continue;

                double xz = trk->getAngleXZ();
                double yz = trk->getAngleYZ();

                if (xz >= xz_range.first && xz <= xz_range.second &&
                    yz >= yz_range.first && yz <= yz_range.second) {
                    eventPassed = true;
                    break;
                }
            }

            // Store weight in "passed" if success (or 0 if fail? No, typical N_eff logic uses vectors)
            // To match your python logic:
            // - "weights" list contains ALL weights.
            // - "passed" list contains booleans.
            // I will implement N_eff calculation exactly as typically done in Python.
            if (eventPassed) {
                w_pass_vec[tt].push_back(weight);
            }
        }
    }

    // 5. Compute Statistics (Bayesian / Beta Distribution)
    std::map<int, EfficiencyResult> results;

    // Beta quantile constants
    double low_q = 0.15865;
    double high_q = 0.84135;

    for (int tt : trackTypes) {
        // Calculate Sums
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

        // N_eff calculation: (Sum W)^2 / Sum (W^2)
        // Note: Using total weights for N_eff is standard for efficiencies
        double n_eff = (sum_w2_total > 0) ? std::pow(sum_w_total, 2) / sum_w2_total : 0.0;

        // Point estimate
        double p_hat = (sum_w_total > 0) ? sum_w_pass / sum_w_total : 0.0;

        // k_eff
        double k_eff = p_hat * n_eff;

        // Bayesian Beta Parameters
        double a = k_eff + 1.0;
        double b = n_eff - k_eff + 1.0;

        // Quantiles using ROOT::Math
        double lo_b = ROOT::Math::beta_quantile(low_q, a, b);
        double hi_b = ROOT::Math::beta_quantile(high_q, a, b);

        double errUp = hi_b - p_hat;
        double errLow = p_hat - lo_b;

        results[tt] = {p_hat, (errUp + errLow) / 2.0};

        std::cout << "TrackType " << tt << ": "
                  << "Eff=" << p_hat << ", "
                  << "N_eff=" << n_eff << ", "
                  << "Err=" << results[tt].error << std::endl;
    }

    return results;
}
