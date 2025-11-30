#include <iostream>
#include <algorithm>

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
#include "trkeffUtils.h"
#include "histograms.h"
#include "trkeffTC.h"



EffResults computeTrackingEfficiencies(const TrkeffConfig& cfg) {
    EffResults results;
    for (auto &p : results) p = {0.0, 0.0}; // Initialize to zeros
    auto startTime = std::chrono::steady_clock::now();

    // ----------------------------
    //  Load histogram parameters
    // ----------------------------
    TEnv hParams(cfg.histParamsFile);
    if (hParams.GetValue("xy.xAxis.bins", -1) == -1) {
        std::cerr << "Failed to read xy.xAxis.bins from " << cfg.histParamsFile.Data() << std::endl;
        exit(-1);
    } else {
        std::cout << "TEnv loaded successfully!" << std::endl;
    }

    // ----------------------------
    //  Load input chain
    // ----------------------------
    TChain *ch = loadChainWithFallback(cfg.inputStr);
    if (!ch) {
        std::cerr << "No valid tree found!" << std::endl;
        return results;
    }
    long nEntries = ch->GetEntries();
    long nBreak = std::min(cfg.nBreak, nEntries);
    bool isMC = hasBranch(ch, "MCTrack");

    // ----------------------------
    //  Set branch addresses
    // ----------------------------
    TClonesArray *tracks = new TClonesArray("sndRecoTrack");
    ch->SetBranchAddress("Reco_MuonTracks", &tracks);

    TClonesArray *sfHits = new TClonesArray("sndScifiHit");
    ch->SetBranchAddress("Digi_ScifiHits", &sfHits);

    TClonesArray *mfHits = new TClonesArray("MuFilterHit");
    ch->SetBranchAddress("Digi_MuFilterHits", &mfHits);

    SNDLHCEventHeader *eventHeader = new SNDLHCEventHeader();
    ch->SetBranchAddress("EventHeader.", &eventHeader);

    TClonesArray *mcTracks = nullptr;
    if (isMC) {
        mcTracks = new TClonesArray("ShipMCTrack");
        ch->SetBranchAddress("MCTrack", &mcTracks);
    }

    // ----------------------------
    //  Prepare histograms
    // ----------------------------
    ch->GetEntry(0);
    int runNum = eventHeader->GetFillNumber() == 0 ?
        0 : eventHeader->GetRunId();

    for (auto& [name, creator] : histFactory)
        histGroups[name] = creator(runNum, trackTypes, hParams);


    // ----------------------------
    //  Geometry
    // ----------------------------
    TPython::Exec("from SndlhcGeo import GeoInterface");
    TPython::Exec(Form("geo = GeoInterface('%s')", cfg.geoFile.Data()));

    MuFilter* MufiDet = (MuFilter*) gROOT->GetListOfGlobals()->FindObject("MuFilter");
    MufiDet->InitEvent(eventHeader);

    int sfst = 0;
    int dsst = 0;
    // ----------------------------
    //  Event loop
    // ----------------------------
    std::cout << "Starting event loop..." << std::endl;
    int lastStatusPercentage = -1;
    double weight = 1.0;
    for (unsigned long i_entry = 0; i_entry < nEntries; ++i_entry) {
        if (i_entry==nBreak) break;

        // Print status every 5%
        int currentStatusPercentage = i_entry*100 / nBreak;
        if (currentStatusPercentage % 5 == 0 && currentStatusPercentage != lastStatusPercentage){
            std::cout << currentStatusPercentage << " %" << std::endl;
            lastStatusPercentage = currentStatusPercentage;
        }
        ch->GetEntry(i_entry);
        if (!eventHeader->isIP1()) continue;

        if (isMC) {
            weight = dynamic_cast<ShipMCTrack*>(mcTracks->At(0))->GetWeight();
        }

        for (const auto& [candTrackType, i_candTrackType] : trackTypeToIndex) {

            for (unsigned i_tagTrack = 0; i_tagTrack < tracks->GetEntries(); ++i_tagTrack) {
                sndRecoTrack* tagTrack = dynamic_cast<sndRecoTrack*>(tracks->At(i_tagTrack));

                int tagTrackType = tagTrack->getTrackType();
                int i_tagTrackType = trackTypeToIndex.at(tagTrackType);

                if (oppositeTrackType.at(candTrackType) != tagTrackType)
                    continue;

                double z = cfg.zRef[i_candTrackType];
                TVector3 tagRefPoint = tagTrack->getPointAtZ(z);
                double xTag = tagRefPoint.X();
                double yTag = tagRefPoint.Y();

                double xz = tagTrack->getAngleXZ();
                double yz = tagTrack->getAngleYZ();

                std::array<int, 5> sfHitCounts = getNsfHits(sfHits);
                std::array<int, 4> dsHitCounts = getNdsHits(mfHits);
                int nSfHitsAll = std::accumulate(sfHitCounts.begin(), sfHitCounts.end(), 0);
                int nDsHitsAll = std::accumulate(dsHitCounts.begin(), dsHitCounts.end(), 0);
                int nSfHitsMax = *std::max_element(sfHitCounts.begin(), sfHitCounts.end());
                int nDsHitsMax = *std::max_element(dsHitCounts.begin(), dsHitCounts.end());

                if (!trackIsConverged(tagTrack)) continue;
                if (!trackIsWithinAngleRange(tagTrack, cfg.xzMin, cfg.xzMax, cfg.yzMin, cfg.yzMax)) continue;



                if (tagTrackType==3 || tagTrackType==13){
                    if (!isNearVetoBar(tagTrack, mfHits, cfg.vetoBarDistance)) continue;
                } else if (tagTrackType==1 || tagTrackType==11) {
                    if (!isNearUS5Bar(tagTrack, mfHits, cfg.us5BarDistance)) continue;
                }
                else continue;

                getHist2D("x.y", i_candTrackType, false)->Fill(xTag, yTag, weight);
                getHist1D("x",   i_candTrackType, false)->Fill(xTag, weight);
                getHist1D("y",   i_candTrackType, false)->Fill(yTag, weight);

                double chi2ndf = tagTrack->getChi2Ndf();
                if (trackIsWithinArea(tagTrack, cfg.zRef.at(i_candTrackType), cfg.xmin, cfg.xmax, cfg.ymin, cfg.ymax)) {


                    getHist1D("chi2ndf",   i_candTrackType, false)->Fill(chi2ndf, weight);

                    if (candTrackType==1 or candTrackType==11){
                        getHist1D("NhitsAll",   i_candTrackType, false)->Fill(nSfHitsAll, weight);
                        getHist1D("NhitsPlane",   i_candTrackType, false)->Fill(nSfHitsMax, weight);
                        getHist2D("x.NhitsAll",   i_candTrackType, false)->Fill(xTag, nSfHitsAll, weight);
                        getHist2D("x.NhitsPlane",   i_candTrackType, false)->Fill(xTag, nSfHitsMax, weight);
                        getHist2D("y.NhitsAll",   i_candTrackType, false)->Fill(yTag, nSfHitsAll, weight);
                        getHist2D("y.NhitsPlane",   i_candTrackType, false)->Fill(yTag, nSfHitsMax, weight);
                    } else {
                        getHist1D("NhitsAll",   i_candTrackType, false)->Fill(nDsHitsAll, weight);
                        getHist1D("NhitsPlane",   i_candTrackType, false)->Fill(nDsHitsMax, weight);
                        getHist2D("x.NhitsAll",   i_candTrackType, false)->Fill(xTag, nDsHitsAll, weight);
                        getHist2D("x.NhitsPlane",   i_candTrackType, false)->Fill(xTag, nDsHitsMax, weight);
                        getHist2D("y.NhitsAll",   i_candTrackType, false)->Fill(yTag, nDsHitsAll, weight);
                        getHist2D("y.NhitsPlane",   i_candTrackType, false)->Fill(yTag, nDsHitsMax, weight);
                    }
                }


                for (unsigned i_candTrack = 0; i_candTrack < tracks->GetEntries(); ++i_candTrack) {
                    auto candTrack = dynamic_cast<sndRecoTrack*>(tracks->At(i_candTrack));

                    if (
                        candTrack->getTrackType() != candTrackType ||
                        !trackIsConverged(candTrack)
                    ) continue;

                    TVector3 candRefPoint = candTrack->getPointAtZ(z);
                    double xCand = candRefPoint.X();
                    double yCand = candRefPoint.Y();

                    if (
                        TMath::Abs(xTag - xCand) > cfg.scifiToDSTrackDistance ||
                        TMath::Abs(yTag - yCand) > cfg.scifiToDSTrackDistance
                    ) continue;

                    getHist2D("x.y", i_candTrackType, true)->Fill(xTag, yTag, weight);
                    getHist1D("x",   i_candTrackType, true)->Fill(xTag, weight);
                    getHist1D("y",   i_candTrackType, true)->Fill(yTag, weight);

                    if (trackIsWithinArea(tagTrack, cfg.zRef.at(i_candTrackType), cfg.xmin, cfg.xmax, cfg.ymin, cfg.ymax)){
                        getHist1D("chi2ndf",   i_candTrackType, true)->Fill(chi2ndf, weight);

                        if (candTrackType==1 or candTrackType==11){
                            getHist1D("NhitsAll",   i_candTrackType, true)->Fill(nSfHitsAll, weight);
                            getHist1D("NhitsPlane",   i_candTrackType, true)->Fill(nSfHitsMax, weight);
                            getHist2D("x.NhitsAll",   i_candTrackType, true)->Fill(xTag, nSfHitsAll, weight);
                            getHist2D("x.NhitsPlane",   i_candTrackType, true)->Fill(xTag, nSfHitsMax, weight);
                            getHist2D("y.NhitsAll",   i_candTrackType, true)->Fill(yTag, nSfHitsAll, weight);
                            getHist2D("y.NhitsPlane",   i_candTrackType, true)->Fill(yTag, nSfHitsMax, weight);
                        } else {
                            getHist1D("NhitsAll",   i_candTrackType, true)->Fill(nDsHitsAll, weight);
                            getHist1D("NhitsPlane",   i_candTrackType, true)->Fill(nDsHitsMax, weight);
                            getHist2D("x.NhitsAll",   i_candTrackType, true)->Fill(xTag, nDsHitsAll, weight);
                            getHist2D("x.NhitsPlane",   i_candTrackType, true)->Fill(xTag, nDsHitsMax, weight);
                            getHist2D("y.NhitsAll",   i_candTrackType, true)->Fill(yTag, nDsHitsAll, weight);
                            getHist2D("y.NhitsPlane",   i_candTrackType, true)->Fill(yTag, nDsHitsMax, weight);
                        }
                    }
                }
            }
        }
    }
    std::cout << "Event loop completed." << std::endl;

    // ----------------------------
    //  Extract TEfficiencies and Save to Output Files
    // ----------------------------
    TFile* fout = nullptr;
    TTree* effTree = nullptr;

    if (!cfg.outFileName.IsNull()) {
        fout = new TFile(cfg.outFileName, "UPDATE");
        effTree = (TTree*)fout->Get("trkeff");  // <-- remove 'TTree*' here!
        if (!effTree) {
            effTree = new TTree("trkeff", "Tracking efficiencies");
        }
    }

    results = createAndSaveTEffs(
        runNum, trackTypes,
        cfg.xmin, cfg.xmax, cfg.ymin, cfg.ymax,
        fout, effTree
    );

    if (fout && effTree) {
        fout->cd();
        effTree->Write("", TObject::kWriteDelete);
        fout->Close();
        std::cout << "Saved all objects to " << fout->GetName() << std::endl;
    }

    // ----------------------------
    //  Print Execution Duration and Return Results
    // ----------------------------
    auto endTime = std::chrono::steady_clock::now();
    double elapsedSec = std::chrono::duration<double>(endTime - startTime).count();
    int h, m, s;
    getSecAsHMS(elapsedSec, h, m, s);
    std::cout << "Elapsed time: "
              << h << "h " << m << "m " << s << "s ("
              << elapsedSec << " seconds)" << std::endl;

    return results;
}



// Nontrivial types (std::map) don't work in python-wrappers
std::vector<double> computeTrackingEfficienciesPy(
    TString inputStr,
    TString geoFile,
    TString outFileName,
    TString histParamsFile,
    double xmin=-48.0, double xmax=-10.0, double ymin=19.0, double ymax=48.0,
    double xzMin = -1e12, double xzMax = 1e12, double yzMin = -1e12, double yzMax = 1e12,
    double zRef1 = 430.0, double zRef11 = 430.0, double zRef3 = 450.0, double zRef13 = 450.0,
    double vetoBarDistance = 3.0,
    double us5BarDistance = 3.0,
    double scifiToDSTrackDistance = 3.0,
    long nBreak = 1e7
) {
    TrkeffConfig cfg;
    cfg.inputStr            = inputStr;
    cfg.geoFile             = geoFile;
    cfg.outFileName         = outFileName;
    cfg.xmin                = xmin;
    cfg.xmax                = xmax;
    cfg.ymin                = ymin;
    cfg.ymax                = ymax;
    cfg.xzMin               = xzMin;
    cfg.xzMax               = xzMax;
    cfg.yzMin               = yzMin;
    cfg.yzMax               = yzMax;
    cfg.zRef[0]             = zRef1;
    cfg.zRef[1]             = zRef11;
    cfg.zRef[2]             = zRef3;
    cfg.zRef[3]             = zRef13;
    cfg.vetoBarDistance     = vetoBarDistance;
    cfg.us5BarDistance      = us5BarDistance;
    cfg.scifiToDSTrackDistance = scifiToDSTrackDistance;
    cfg.nBreak              = nBreak;
    cfg.histParamsFile      = histParamsFile;

    auto effs = computeTrackingEfficiencies(cfg);
    std::vector<double> out;
    for (const auto& [e, ee] : effs) {
        out.push_back(e);
        out.push_back(ee);
    }
    return out;
}




int main(int argc, char** argv)
{
    if (argc < 21) {
        std::cerr << "Too few arguments " << argc << "/20" << std::endl;
        return -1;
    }

    TrkeffConfig cfg;
    cfg.inputStr            = argv[1];
    cfg.geoFile             = argv[2];
    cfg.outFileName         = argv[3];
    cfg.xmin                = atof(argv[4]);
    cfg.xmax                = atof(argv[5]);
    cfg.ymin                = atof(argv[6]);
    cfg.ymax                = atof(argv[7]);
    cfg.xzMin               = atof(argv[8]);
    cfg.xzMax               = atof(argv[9]);
    cfg.yzMin               = atof(argv[10]);
    cfg.yzMax               = atof(argv[11]);
    cfg.zRef[0]             = atof(argv[12]);
    cfg.zRef[1]             = atof(argv[13]);
    cfg.zRef[2]             = atof(argv[14]);
    cfg.zRef[3]             = atof(argv[15]);
    cfg.vetoBarDistance     = atof(argv[16]);
    cfg.us5BarDistance      = atof(argv[17]);
    cfg.scifiToDSTrackDistance = atof(argv[18]);
    cfg.nBreak              = atol(argv[19]);
    cfg.histParamsFile      = argv[20];

    auto effResults = computeTrackingEfficiencies(cfg);

    for (size_t i_tt = 0; i_tt < effResults.size(); ++i_tt) {
        std::cout << "Track type " << i_tt << ":\t"
                  << effResults[i_tt].first << " +/- " << effResults[i_tt].second
                  << std::endl;
    }

    return 0;
}
