#include <array>
#include <map>
#include <string>

#include "TObjArray.h"
#include "TString.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TEnv.h"

using HistPair = std::pair<TObjArray*, TObjArray*>;
using HistCreator = HistPair(*)(int, const std::array<int,4>&, const TEnv& hParams);

std::map<std::string, HistPair> histGroups;

HistPair createXYHists(
    int run, const std::array<int, 4>& trackTypes, const TEnv& hParams
)
{
    auto* histsPassed = new TObjArray(trackTypes.size());
    auto* histsTotal  = new TObjArray(trackTypes.size());

    int nBinsX = hParams.GetValue("xy.xAxis.bins", 40);
    int nBinsY = hParams.GetValue("xy.yAxis.bins", 40);
    double minX = hParams.GetValue("xy.xAxis.min", -70);
    double maxX = hParams.GetValue("xy.xAxis.max", 10);
    double minY = hParams.GetValue("xy.yAxis.min", -10);
    double maxY = hParams.GetValue("xy.yAxis.max", 70);

    for (int trackType : trackTypes) {
        TString namePassed  = TString::Format("h_x.y_%d_%d_passed", trackType, run);
        TString titlePassed = TString::Format("Passed %d %d;x (cm);y (cm)", trackType, run);

        TString nameTotal  = TString::Format("h_x.y_%d_%d_total", trackType, run);
        TString titleTotal = TString::Format("Total %d %d;x (cm);y (cm)", trackType, run);

        histsPassed->Add(
            new TH2D(namePassed, titlePassed, nBinsX, minX, maxX, nBinsY, minY, maxY)
        );
        histsTotal->Add(
            new TH2D(nameTotal, titleTotal, nBinsX, minX, maxX, nBinsY, minY, maxY)
        );
    }

    return {histsPassed, histsTotal};
}



HistPair createXHists(
    int run, const std::array<int, 4>& trackTypes, const TEnv& hParams
)
{
    auto* histsPassed = new TObjArray(trackTypes.size());
    auto* histsTotal  = new TObjArray(trackTypes.size());

    int nBinsX = hParams.GetValue("xy.xAxis.bins", 40);
    double minX = hParams.GetValue("xy.xAxis.min", -70);
    double maxX = hParams.GetValue("xy.xAxis.max", 10);

    for (int trackType : trackTypes) {
        TString namePassed  = TString::Format("h1d_x_%d_%d_passed", trackType, run);
        TString titlePassed = TString::Format("Passed %d %d;x (cm);", trackType, run);

        TString nameTotal  = TString::Format("h1d_x_%d_%d_total", trackType, run);
        TString titleTotal = TString::Format("Total %d %d;x (cm);", trackType, run);

        histsPassed->Add(
            new TH1D(namePassed, titlePassed, nBinsX, minX, maxX)
        );
        histsTotal->Add(
            new TH1D(nameTotal, titleTotal, nBinsX, minX, maxX)
        );
    }

    return {histsPassed, histsTotal};
}




 HistPair createYHists(
    int run, const std::array<int, 4>& trackTypes, const TEnv& hParams
)
{
    auto* histsPassed = new TObjArray(trackTypes.size());
    auto* histsTotal  = new TObjArray(trackTypes.size());

    int nBinsY = hParams.GetValue("xy.yAxis.bins", 40);
    double minY = hParams.GetValue("xy.yAxis.min", -10);
    double maxY = hParams.GetValue("xy.yAxis.max", 70);

    for (int trackType : trackTypes) {
        TString namePassed  = TString::Format("h1d_y_%d_%d_passed", trackType, run);
        TString titlePassed = TString::Format("Passed %d %d;y (cm);", trackType, run);

        TString nameTotal  = TString::Format("h1d_y_%d_%d_total", trackType, run);
        TString titleTotal = TString::Format("Total %d %d;y (cm);", trackType, run);

        histsPassed->Add(
            new TH1D(namePassed, titlePassed, nBinsY, minY, maxY)
        );
        histsTotal->Add(
            new TH1D(nameTotal, titleTotal, nBinsY, minY, maxY)
        );
    }

    return {histsPassed, histsTotal};
}




HistPair createNhitsAllHists(
    int run, const std::array<int, 4>& trackTypes, const TEnv& hParams
)
{
    auto* histsPassed = new TObjArray(trackTypes.size());
    auto* histsTotal  = new TObjArray(trackTypes.size());


    for (int trackType : trackTypes) {
        TString namePassed, titlePassed, nameTotal, titleTotal;
        double minN, maxN;
        int nBinsN;
        if (trackType==1 || trackType==11) {
            nBinsN = hParams.GetValue("N.hits.scifi.all.bins", 70);
            minN = hParams.GetValue("N.hits.scifi.all.min", 0);
            maxN = hParams.GetValue("N.hits.scifi.all.max", 140);

            namePassed  = TString::Format("h1d_NsfHitsAll_%d_%d_passed", trackType, run);
            titlePassed = TString::Format("Passed %d %d;#N_{sfHitsTotal};", trackType, run);

            nameTotal   = TString::Format("h1d_NsfHitsAll_%d_%d_total", trackType, run);
            titleTotal  = TString::Format("Total %d %d;#N_{sfHitsTotal};", trackType, run);
        } else  {
            nBinsN = hParams.GetValue("N.hits.ds.all.bins", 70);
            minN = hParams.GetValue("N.hits.ds.all.min", 0);
            maxN = hParams.GetValue("N.hits.ds.all.max", 140);

            namePassed  = TString::Format("h1d_NdsHitsAll_%d_%d_passed", trackType, run);
            titlePassed = TString::Format("Passed %d %d;#N_{dsHitsAll};", trackType, run);

            nameTotal  = TString::Format("h1d_NdsHitsAll_%d_%d_total", trackType, run);
            titleTotal = TString::Format("Total %d %d;#N_{dsHitsAll};", trackType, run);
        }

        histsPassed->Add(
            new TH1D(namePassed, titlePassed, nBinsN, minN, maxN)
        );
        histsTotal->Add(
            new TH1D(nameTotal, titleTotal, nBinsN, minN, maxN)
        );
    }

    return {histsPassed, histsTotal};
}

HistPair createNhitsPlaneHists(
    int run, const std::array<int, 4>& trackTypes, const TEnv& hParams
)
{
    auto* histsPassed = new TObjArray(trackTypes.size());
    auto* histsTotal  = new TObjArray(trackTypes.size());


    for (int trackType : trackTypes) {
        TString namePassed, titlePassed, nameTotal, titleTotal;
        double minN, maxN;
        int nBinsN;
        if (trackType==1 || trackType==11) {
            nBinsN = hParams.GetValue("N.hits.scifi.plane.bins", 60);
            minN = hParams.GetValue("N.hits.scifi.plane.min", 0);
            maxN = hParams.GetValue("N.hits.scifi.plane.max", 120);

            namePassed  = TString::Format("h1d_NsfHitsPlane_%d_%d_passed", trackType, run);
            titlePassed = TString::Format("Passed %d %d;#N_{sfHitsPlane};", trackType, run);

            nameTotal  = TString::Format("h1d_NsfHitsPlane_%d_%d_total", trackType, run);
            titleTotal = TString::Format("Total %d %d;#N_{sfHitsPlane};", trackType, run);
        } else {
            nBinsN = hParams.GetValue("N.hits.ds.plane.bins", 40);
            minN = hParams.GetValue("N.hits.ds.plane.min", 0);
            maxN = hParams.GetValue("N.hits.ds.plane.max", 80);

            namePassed  = TString::Format("h1d_NdsHitsPlane_%d_%d_passed", trackType, run);
            titlePassed = TString::Format("Passed %d %d;#N_{dsHitsPlane};", trackType, run);

            nameTotal  = TString::Format("h1d_NdsHitsPlane_%d_%d_total", trackType, run);
            titleTotal = TString::Format("Total %d %d;#N_{dsHitsPlane};", trackType, run);
        }

        histsPassed->Add(
            new TH1D(namePassed, titlePassed, nBinsN, minN, maxN)
        );
        histsTotal->Add(
            new TH1D(nameTotal, titleTotal, nBinsN, minN, maxN)
        );
    }

    return {histsPassed, histsTotal};
}



HistPair createXvsNhitsAllHists(
    int run, const std::array<int, 4>& trackTypes, const TEnv& hParams
)
{
    auto* histsPassed = new TObjArray(trackTypes.size());
    auto* histsTotal  = new TObjArray(trackTypes.size());


    for (int trackType : trackTypes) {
        TString namePassed, titlePassed, nameTotal, titleTotal;
        double minN, maxN, minX, maxX;
        int nBinsN, nBinsX;
        nBinsX = hParams.GetValue("xy.xAxis.bins", 60);
        minX = hParams.GetValue("xy.xAxis.min", -70);
        maxX = hParams.GetValue("xy.xAxis.max", -10);
        if (trackType==1 || trackType==11) {
            nBinsN = hParams.GetValue("N.hits.scifi.all.bins", 60);
            minN = hParams.GetValue("N.hits.scifi.all.min", 0);
            maxN = hParams.GetValue("N.hits.scifi.all.max", 120);

            namePassed  = TString::Format("h2d_x.NsfHitsAll_%d_%d_passed", trackType, run);
            titlePassed = TString::Format("Passed %d %d;x (cm);#N_{sfHitsAll};", trackType, run);

            nameTotal  = TString::Format("h2d_x.NsfHitsAll_%d_%d_total", trackType, run);
            titleTotal = TString::Format("Total %d %d;x (cm);#N_{sfHitsAll};", trackType, run);
        } else {
            nBinsN = hParams.GetValue("N.hits.ds.all.bins", 40);
            minN = hParams.GetValue("N.hits.ds.all.min", 0);
            maxN = hParams.GetValue("N.hits.ds.all.max", 80);

            namePassed  = TString::Format("h2d_x.NdsHitsAll%d_%d_passed", trackType, run);
            titlePassed = TString::Format("Passed %d %d;x(cm);#N_{dsHitsAll};", trackType, run);

            nameTotal  = TString::Format("h2d_x.NdsHitsAll_%d_%d_total", trackType, run);
            titleTotal = TString::Format("Total %d %d;x (cm);#N_{dsHitsAll};", trackType, run);
        }

        histsPassed->Add(
            new TH2D(namePassed, titlePassed, nBinsX, minX, maxX, nBinsN, minN, maxN)
        );
        histsTotal->Add(
            new TH2D(nameTotal, titleTotal, nBinsX, minX, maxX, nBinsN, minN, maxN)
        );
    }

    return {histsPassed, histsTotal};
}



HistPair createXvsNhitsPlaneHists(
    int run, const std::array<int, 4>& trackTypes, const TEnv& hParams
)
{
    auto* histsPassed = new TObjArray(trackTypes.size());
    auto* histsTotal  = new TObjArray(trackTypes.size());


    for (int trackType : trackTypes) {
        TString namePassed, titlePassed, nameTotal, titleTotal;
        double minN, maxN, minX, maxX;
        int nBinsN, nBinsX;
        nBinsX = hParams.GetValue("xy.xAxis.bins", 60);
        minX = hParams.GetValue("xy.xAxis.min", -70);
        maxX = hParams.GetValue("xy.xAxis.max", -10);
        if (trackType==1 || trackType==11) {
            nBinsN = hParams.GetValue("N.hits.scifi.plane.bins", 60);
            minN = hParams.GetValue("N.hits.scifi.plane.min", 0);
            maxN = hParams.GetValue("N.hits.scifi.plane.max", 120);

            namePassed  = TString::Format("h2d_x.NsfHitsPlane_%d_%d_passed", trackType, run);
            titlePassed = TString::Format("Passed %d %d;x (cm);#N_{sfHitsPlane};", trackType, run);

            nameTotal  = TString::Format("h2d_x.NsfHitsPlane_%d_%d_total", trackType, run);
            titleTotal = TString::Format("Total %d %d;x (cm);#N_{sfHitsPlane};", trackType, run);
        } else {
            nBinsN = hParams.GetValue("N.hits.ds.plane.bins", 40);
            minN = hParams.GetValue("N.hits.ds.plane.min", 0);
            maxN = hParams.GetValue("N.hits.ds.plane.max", 80);

            namePassed  = TString::Format("h2d_x.NdsHitsPlane_%d_%d_passed", trackType, run);
            titlePassed = TString::Format("Passed %d %d;x(cm);#N_{dsHitsPlane};", trackType, run);

            nameTotal  = TString::Format("h2d_x.NdsHitsPlane_%d_%d_total", trackType, run);
            titleTotal = TString::Format("Total %d %d;x (cm);#N_{dsHitsPlane};", trackType, run);
        }

        histsPassed->Add(
            new TH2D(namePassed, titlePassed, nBinsX, minX, maxX, nBinsN, minN, maxN)
        );
        histsTotal->Add(
            new TH2D(nameTotal, titleTotal, nBinsX, minX, maxX, nBinsN, minN, maxN)
        );
    }

    return {histsPassed, histsTotal};
}


HistPair createYvsNhitsAllHists(
    int run, const std::array<int, 4>& trackTypes, const TEnv& hParams
)
{
    auto* histsPassed = new TObjArray(trackTypes.size());
    auto* histsTotal  = new TObjArray(trackTypes.size());


    for (int trackType : trackTypes) {
        TString namePassed, titlePassed, nameTotal, titleTotal;
        double minN, maxN, minY, maxY;
        int nBinsN, nBinsY;
        nBinsY = hParams.GetValue("xy.yAxis.bins", 60);
        minY = hParams.GetValue("xy.yAxis.min", -70);
        maxY = hParams.GetValue("xy.yAxis.max", -10);
        if (trackType==1 || trackType==11) {
            nBinsN = hParams.GetValue("N.hits.scifi.all.bins", 60);
            minN = hParams.GetValue("N.hits.scifi.all.min", 0);
            maxN = hParams.GetValue("N.hits.scifi.all.max", 120);

            namePassed  = TString::Format("h2d_y.NsfHitsAll_%d_%d_passed", trackType, run);
            titlePassed = TString::Format("Passed %d %d;y (cm);#N_{sfHitsAll};", trackType, run);

            nameTotal  = TString::Format("h2d_y.NsfHitsAll_%d_%d_total", trackType, run);
            titleTotal = TString::Format("Total %d %d;y (cm);#N_{sfHitsAll};", trackType, run);
        } else {
            nBinsN = hParams.GetValue("N.hits.ds.all.bins", 40);
            minN = hParams.GetValue("N.hits.ds.all.min", 0);
            maxN = hParams.GetValue("N.hits.ds.all.max", 80);

            namePassed  = TString::Format("h2d_y.NdsHitsAll_%d_%d_passed", trackType, run);
            titlePassed = TString::Format("Passed %d %d;y(cm);#N_{dsHitsAll};", trackType, run);

            nameTotal  = TString::Format("h2d_y.NdsHitsAll_%d_%d_total", trackType, run);
            titleTotal = TString::Format("Total %d %d;y (cm);#N_{dsHitsAll};", trackType, run);
        }

        histsPassed->Add(
            new TH2D(namePassed, titlePassed, nBinsY, minY, maxY, nBinsN, minN, maxN)
        );
        histsTotal->Add(
            new TH2D(nameTotal, titleTotal, nBinsY, minY, maxY, nBinsN, minN, maxN)
        );
    }

    return {histsPassed, histsTotal};
}



HistPair createYvsNhitsPlaneHists(
    int run, const std::array<int, 4>& trackTypes, const TEnv& hParams
)
{
    auto* histsPassed = new TObjArray(trackTypes.size());
    auto* histsTotal  = new TObjArray(trackTypes.size());


    for (int trackType : trackTypes) {
        TString namePassed, titlePassed, nameTotal, titleTotal;
        double minN, maxN, minY, maxY;
        int nBinsN, nBinsY;
        nBinsY = hParams.GetValue("xy.yAxis.bins", 60);
        minY = hParams.GetValue("xy.yAxis.min", -70);
        maxY = hParams.GetValue("xy.yAxis.max", -10);
        if (trackType==1 || trackType==11) {
            nBinsN = hParams.GetValue("N.hits.scifi.plane.bins", 60);
            minN = hParams.GetValue("N.hits.scifi.plane.min", 0);
            maxN = hParams.GetValue("N.hits.scifi.plane.max", 120);

            namePassed  = TString::Format("h2d_y.NsfHitsPlane_%d_%d_passed", trackType, run);
            titlePassed = TString::Format("Passed %d %d;y (cm);#N_{sfHitsPlane};", trackType, run);

            nameTotal  = TString::Format("h2d_y.NsfHitsPlane_%d_%d_total", trackType, run);
            titleTotal = TString::Format("Total %d %d;y (cm);#N_{sfHitsPlane};", trackType, run);
        } else {
            nBinsN = hParams.GetValue("N.hits.ds.plane.bins", 40);
            minN = hParams.GetValue("N.hits.ds.plane.min", 0);
            maxN = hParams.GetValue("N.hits.ds.plane.max", 80);

            namePassed  = TString::Format("h2d_y.NdsHitsPlane_%d_%d_passed", trackType, run);
            titlePassed = TString::Format("Passed %d %d;y(cm);#N_{dsHitsPlane};", trackType, run);

            nameTotal  = TString::Format("h2d_y.NdsHitsPlane_%d_%d_total", trackType, run);
            titleTotal = TString::Format("Total %d %d;y (cm);#N_{dsHitsPlane};", trackType, run);
        }

        histsPassed->Add(
            new TH2D(namePassed, titlePassed, nBinsY, minY, maxY, nBinsN, minN, maxN)
        );
        histsTotal->Add(
            new TH2D(nameTotal, titleTotal, nBinsY, minY, maxY, nBinsN, minN, maxN)
        );
    }

    return {histsPassed, histsTotal};
}



HistPair createChi2ndfHists(
    int run, const std::array<int, 4>& trackTypes, const TEnv& hParams
)
{
    auto* histsPassed = new TObjArray(trackTypes.size());
    auto* histsTotal  = new TObjArray(trackTypes.size());

    for (int trackType : trackTypes) {
        TString namePassed  = TString::Format("h1d_chi2ndf_%d_%d_passed", trackType, run);
        TString titlePassed = TString::Format("Passed %d %d;Tagging track's #Chi^{2}/ndf;", trackType, run);

        TString nameTotal  = TString::Format("h1d_chi2ndf_%d_%d_total", trackType, run);
        TString titleTotal = TString::Format("Total %d %d;Tagging track's #Chi^{2}/ndf;", trackType, run);

        int nBinsChi2ndf = hParams.GetValue(TString::Format("chi2ndf.%d.bins", trackType), 100);
        double minChi2ndf = hParams.GetValue(TString::Format("chi2ndf.%d.min", trackType), 0);
        double maxChi2ndf = hParams.GetValue(TString::Format("chi2ndf.%d.max", trackType), 100);

        histsPassed->Add(
            new TH1D(namePassed, titlePassed, nBinsChi2ndf, minChi2ndf, maxChi2ndf)
        );
        histsTotal->Add(
            new TH1D(nameTotal, titleTotal, nBinsChi2ndf, minChi2ndf, maxChi2ndf)
        );
    }

    return {histsPassed, histsTotal};
}





std::map<std::string, HistCreator> histFactory = {
    { "x.y",          createXYHists },
    { "x",            createXHists  },
    { "y",            createYHists  },
    { "NhitsAll",     createNhitsAllHists},
    { "NhitsPlane",   createNhitsPlaneHists},
    { "x.NhitsAll",   createXvsNhitsAllHists},
    { "x.NhitsPlane", createXvsNhitsPlaneHists},
    { "y.NhitsAll",   createYvsNhitsAllHists},
    { "y.NhitsPlane", createYvsNhitsPlaneHists},
    { "chi2ndf",      createChi2ndfHists},
};


int getHGroupDim(const std::string& s) {
    return std::count(s.begin(), s.end(), '.') + 1;
}

TH2D* getHist2D(const std::string& key, int idx, bool passed) {
    auto& pair = histGroups.at(key);
    if (passed) {
        return dynamic_cast<TH2D*>(pair.first->At(idx));
    } else {
        return dynamic_cast<TH2D*>(pair.second->At(idx));
    }
}

TH1D* getHist1D(const std::string& key, int idx, bool passed) {
    auto& pair = histGroups.at(key);
    if (passed) {
        return dynamic_cast<TH1D*>(pair.first->At(idx));
    } else {
        return dynamic_cast<TH1D*>(pair.second->At(idx));
    }
}
