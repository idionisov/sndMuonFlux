#include <iostream>
#include <cmath>
#include <vector>
#include <map>

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
#include "TClonesArray.h"
#include "Math/QuantFunc.h"

#include "MuFilter.h"
#include "MuFilterHit.h"
#include "MuFilterPoint.h"
#include "sndScifiHit.h"
#include "ScifiPoint.h"
#include "SNDLHCEventHeader.h"
#include "sndRecoTrack.h"
#include "ShipMCTrack.h"
#include "muonFluxUtils.h"
#include "trkeffMCT.h"

bool ThereIsAMuon(SNDLHCEventHeader* header, TClonesArray* mcTracks) {
    for (int i = 0; i < mcTracks->GetEntriesFast(); ++i) {
        ShipMCTrack* trk = (ShipMCTrack*)mcTracks->At(i);
        if (trk->GetPdgCode() == 13) return true;
    }
    return false;
}

bool SfTrackIsReconstructible(TClonesArray* scifiPoints) {
    std::map<int, int> nPointsV;
    std::map<int, int> nPointsH;

    for (int i = 0; i < scifiPoints->GetEntriesFast(); ++i) {
        ScifiPoint* p = (ScifiPoint*)scifiPoints->At(i);
        if (p->PdgCode() != 13 || p->GetTrackID() != 0) continue;

        int detID = p->GetDetectorID();
        int station = detID / 1000000;
        int planeType = (detID / 100000) % 2; // 0: H, 1: V

        if (planeType == 0) nPointsH[station]++;
        else nPointsV[station]++;
    }

    int nV = 0, nH = 0;
    for (int s = 1; s <= 5; ++s) {
        if (nPointsV[s] > 0) nV++;
        if (nPointsH[s] > 0) nH++;
    }
    return (nV >= 3 && nH >= 3);
}

bool DsTrackIsReconstructible(TClonesArray* muFilterPoints) {
    std::map<int, int> nPointsV;
    std::map<int, int> nPointsH;

    for (int i = 0; i < muFilterPoints->GetEntriesFast(); ++i) {
        MuFilterPoint* p = (MuFilterPoint*)muFilterPoints->At(i);
        if (p->PdgCode() != 13 || p->GetTrackID() != 0) continue;

        int detID = p->GetDetectorID();
        if (detID < 30000 || detID > 34999) continue;

        int station = (detID / 1000) % 10 + 1;
        int barNum = detID % 1000;

        if (barNum > 59) nPointsV[station]++;
        else nPointsH[station]++;
    }

    int nV = 0, nH = 0;
    for (int s = 1; s <= 4; ++s) {
        if (nPointsV[s] > 0) nV++;
        if (nPointsH[s] > 0) nH++;
    }
    return (nV >= 3 && nH >= 3);
}

bool McTrackCrossedFiducialArea(TClonesArray* mcTracks, double zRef, double xmin, double xmax, double ymin, double ymax) {
    for (int i = 0; i < mcTracks->GetEntriesFast(); ++i) {
        ShipMCTrack* trk = (ShipMCTrack*)mcTracks->At(i);
        if (trk->GetMotherId() == -1 && trk->GetPdgCode() == 13) {
            double px = trk->GetPx();
            double py = trk->GetPy();
            double pz = trk->GetPz();

            double x0 = trk->GetStartX();
            double y0 = trk->GetStartY();
            double z0 = trk->GetStartZ();

            if (pz == 0) continue;
            double t = (zRef - z0) / pz;
            double x = x0 + px * t;
            double y = y0 + py * t;

            if (x >= xmin && x <= xmax && y >= ymin && y <= ymax) {
                return true;
            }
        }
    }
    return false;
}

std::vector<double> computeTrackingEfficiencies_MCT(
    TString input_files,
    double sigma,
    double col_rate,
    double L_LHC,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    double zRef1,
    double zRef11,
    double zRef3,
    double zRef13,
    double xzMin,
    double xzMax,
    double yzMin,
    double yzMax
) {
    // 1. Constants & Setup
    double L_MC = col_rate / sigma;
    double Ks = L_LHC / L_MC;
    std::vector<int> trackTypes = {1, 11, 3, 13};

    // Calculate Area (A)
    double A = (xmax - xmin) * (ymax - ymin);

    // Load Chain
    TChain* ch = new TChain("cbmsim");
    ch->Add(input_files.Data());
    Long64_t n_entries = ch->GetEntries();

    std::cout << "Entries    : " << n_entries << std::endl;
    std::cout << "Luminosity : " << L_MC << std::endl;
    std::cout << "Area       : " << A << std::endl;
    std::cout << "Scaling    : " << Ks << std::endl;

    // 2. Branch Setup
    SNDLHCEventHeader* eventHeader = nullptr;
    TClonesArray* mcTracks = nullptr;
    TClonesArray* recoTracks = nullptr;
    TClonesArray* scifiPoints = nullptr;
    TClonesArray* muFilterPoints = nullptr;

    ch->SetBranchAddress("EventHeader.", &eventHeader);
    ch->SetBranchAddress("MCTrack", &mcTracks);
    ch->SetBranchAddress("Reco_MuonTracks", &recoTracks);
    ch->SetBranchAddress("ScifiPoint", &scifiPoints);
    ch->SetBranchAddress("MuFilterPoint", &muFilterPoints);

    // 3. Storage for weights and pass/fail
    std::map<int, std::vector<double>> w_total_vec;
    std::map<int, std::vector<double>> w_pass_vec;

    for(int tt : trackTypes) {
        w_total_vec[tt] = std::vector<double>();
        w_pass_vec[tt]  = std::vector<double>();
    }

    // 4. Event Loop
    int next_print = 0;
    for (Long64_t i = 0; i < n_entries; ++i) {
        ch->GetEntry(i);

        int progress = (i * 100) / n_entries;
        if (progress >= next_print) {
            std::cout << progress << " %" << std::endl;
            next_print += 5;
        }

        if (!eventHeader->isIP1() || !ThereIsAMuon(eventHeader, mcTracks)) {
            continue;
        }

        bool _sf = SfTrackIsReconstructible(scifiPoints);
        bool _ds = DsTrackIsReconstructible(muFilterPoints);

        if (!_sf && !_ds) continue;

        std::map<int, bool> reco_map;
        reco_map[1]  = _sf && McTrackCrossedFiducialArea(mcTracks, zRef1,  xmin, xmax, ymin, ymax);
        reco_map[11] = _sf && McTrackCrossedFiducialArea(mcTracks, zRef11, xmin, xmax, ymin, ymax);
        reco_map[3]  = _ds && McTrackCrossedFiducialArea(mcTracks, zRef3,  xmin, xmax, ymin, ymax);
        reco_map[13] = _ds && McTrackCrossedFiducialArea(mcTracks, zRef13, xmin, xmax, ymin, ymax);

        double weight = 1.0;
        if (mcTracks->GetEntriesFast() > 0) {
            ShipMCTrack* mc = (ShipMCTrack*)mcTracks->At(0);
            weight = mc->GetWeight() * Ks;
        }

        for (int tt : trackTypes) {
            if (!reco_map[tt]) continue;

            w_total_vec[tt].push_back(weight);

            bool eventPassed = false;
            for (int j = 0; j < recoTracks->GetEntriesFast(); ++j) {
                sndRecoTrack* trk = (sndRecoTrack*)recoTracks->At(j);
                if (trk->getTrackType() != tt) continue;

                double xz = trk->getAngleXZ() * 1e3;
                double yz = trk->getAngleYZ() * 1e3;

                if (xz >= xzMin && xz <= xzMax &&
                    yz >= yzMin && yz <= yzMax) {
                    eventPassed = true;
                    break;
                }
            }

            if (eventPassed) {
                w_pass_vec[tt].push_back(weight);
            }
        }
    }

    // 5. Compute Statistics
    std::map<int, std::vector<double>> results;
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

        double n_eff = (sum_w2_total > 0) ? std::pow(sum_w_total, 2) / sum_w2_total : 0.0;
        double p_hat = (sum_w_total > 0) ? sum_w_pass / sum_w_total : 0.0;
        double k_eff = p_hat * n_eff;

        double a = k_eff + 1.0;
        double b = n_eff - k_eff + 1.0;

        double lo_b = ROOT::Math::beta_quantile(low_q, a, b);
        double hi_b = ROOT::Math::beta_quantile(high_q, a, b);

        double errUp = hi_b - p_hat;
        double errLow = p_hat - lo_b;

        results[tt] = {p_hat, (errUp + errLow) / 2.0};

        // std::cout << "TrackType " << tt << ": "
        //           << "Eff=" << p_hat << ", "
        //           << "N_eff=" << n_eff << ", "
        //           << "Err=" << results[tt][1] << std::endl;
    }

    std::vector<double> final_results;
    for (int tt : {1, 11, 3, 13}) {
        if (results.count(tt)) {
            final_results.push_back(results[tt][0]);
            final_results.push_back(results[tt][1]);
        } else {
            final_results.push_back(0.0);
            final_results.push_back(0.0);
        }
    }

    return final_results;
}
