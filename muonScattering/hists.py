import ROOT
from sndUtils import sys, alg, system, algorithm, att

h_ThetaRms = ROOT.TH1F("h_ThetaRMS_sim", ";Entries;Outgoing angle root mean square $theta_{RMS} [mrad]", 30, 0, 10)

pr_xz_thetaRms = ROOT.TProfile("pr_ThetaXZ.ThetaRMS-XZ_sim", ";#theta_{{XZ}} [mrad];#theta_{{XZ}} _{{RMS}} [mrad]", 20, -20, 20)
pr_yz_thetaRms = ROOT.TProfile("pr_ThetaYZ.ThetaRMS-YZ_sim", ";#theta_{{YZ}} [mrad];#theta_{{YZ}} _{{RMS}} [mrad]", 20, -20, 20)
h_xz_thetaRms  = ROOT.TH2F("h_ThetaXZ.ThetaRMS-XZ_sim", ";#theta_{{XZ}} [mrad];#theta_{{XZ}} _{{RMS}} [mrad]", 20, -20, 20, 20, 1e-6, 1e-5)
h_yz_thetaRms  = ROOT.TH2F("h_ThetaYZ.ThetaRMS-YZ_sim", ";#theta_{{YZ}} [mrad];#theta_{{YZ}} _{{RMS}} [mrad]", 20, -20, 20, 20, 1e-6, 1e-5)


pr_xz_scat = {
    tt: ROOT.TProfile(
        f"pr_xz.scat_{tt}_sim",
        f"{system(tt)} {algorithm(tt)};#theta_{{XZ}} [mrad];Probability of misdirection",
        20, -20, 20
    ) for tt in (3, 13)
}
pr_yz_scat = {
    tt: ROOT.TProfile(
        f"pr_yz.scat_{tt}_sim",
        f"{system(tt)} {algorithm(tt)};#theta_{{YZ}} [mrad];Probability of misdirection",
        20, -20, 20
    ) for tt in (3, 13)
}
pr_xz_yz_scat = {
    tt: ROOT.TProfile2D(
        f"pr_xz.yz.scat_{tt}_sim",
        f"{system(tt)} {algorithm(tt)};#theta_{{XZ}} [mrad];#theta_{{YZ}} [mrad];Probability of misdirection",
        20, -20, 20, 20, -20, 20
    ) for tt in (3, 13)
}
h_xz_scat  = {
    tt: ROOT.TH2F(
        f"h_xz.scat_{tt}_sim",
        f"{system(tt)} {algorithm(tt)};#theta_{{XZ}} [mrad];Probability of misdirection",
        20, -20, 20, 20, 0, 1
    ) for tt in (3, 13)
}
h_yz_scat = {
    tt: ROOT.TH2F(
        f"h_yz.scat_{tt}_sim",
        f"{system(tt)} {algorithm(tt)};#theta_{{YZ}} [mrad];Probability of misdirection",
        20, -20, 20, 20, 0, 1
    ) for tt in (3, 13)
}

pr_x_scat = {
    tt: ROOT.TProfile(
        f"pr_x.scat_{tt}_sim",
        f"{system(tt)} {algorithm(tt)};X_{{ref}} [cm];Probability of misdirection",
        20, -60, 0
    ) for tt in (3, 13)
}
pr_y_scat = {
    tt: ROOT.TProfile(
        f"pr_y.scat_{tt}_sim",
        f"{system(tt)} {algorithm(tt)};Y_{{ref}} [cm];Probability of misdirection",
        20, 0, 60
    ) for tt in (3, 13)
}
pr_x_y_scat = {
    tt: ROOT.TProfile2D(
        f"pr_x.y.scat_{tt}_sim",
        f"{system(tt)} {algorithm(tt)};X_{{ref}} [cm];Y_{{ref}} [cm];Probability of misdirection",
        20, -60, 0, 20, 0, 60
    ) for tt in (3, 13)
}
h_x_scat  = {
    tt: ROOT.TH2F(
        f"h_x.scat_{tt}_sim",
        f"{system(tt)} {algorithm(tt)};X_{{ref}} [cm];Probability of misdirection",
        20, -60, 0, 20, 0, 1
    ) for tt in (3, 13)
}
h_y_scat = {
    tt: ROOT.TH2F(
        f"h_y.scat_{tt}_sim",
        f"{system(tt)} {algorithm(tt)};Y_{{ref}} [cm];Probability of misdirection",
        20, 0, 60, 20, 0, 1
    ) for tt in (3, 13)
}
