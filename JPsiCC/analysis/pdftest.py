# First run . env.sh. It is in the RareHiggsRun3 directory.

from kytools import rootpdf
import ROOT

mcbkg0_rdf = ROOT.RDataFrame('Events', 'snapshot_MC_BKG_2018_GF_v202406_20240719.root')
mcsig_rdf = ROOT.RDataFrame('Events', 'snapshot_MC_SIG_2018_GF_v202406_20240719.root')

mcbkg0_rdf = mcbkg0_rdf.Define('mJpsi', 'Jpsi_kin_mass[0]')
mcsig_rdf = mcsig_rdf.Define('mJpsi', 'Jpsi_kin_mass[0]')
mcbkg1_rdf = mcbkg0_rdf.Filter('sample=="BToJpsi_JPsiToMuMu_BMuonFilter_HardQCD_TuneCP5_13TeV-pythia8-evtgen+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1+MINIAODSIM"')
mcbkg2_rdf = mcbkg0_rdf.Filter('sample=="JpsiToMuMu_JpsiPt8_TuneCP5_13TeV-pythia8+RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2+MINIAODSIM"')

fit = rootpdf.FittingTool(2018, 'GF', '202406', 'mJpsi', 'm_{#mu^{+}#mu^{-}}', 2., 4.)
fit.makeSignalPDF('gaussian')
fit.makeBackgroundPDF('gaussian')
fit.loadRooDataSet(mcsig_rdf, 'MC_SIG')
fit.loadRooDataHist(mcbkg0_rdf, 'MC_BKG0')
fit.loadRooDataHist(mcbkg1_rdf, 'MC_BKG1')
fit.loadRooDataHist(mcbkg2_rdf, 'MC_BKG2')

fit.fit('MC_SIG', 'gaussian', binned=False)
fit.fit('MC_BKG0', 'gaussian', binned=True)

fit.plot(mask=['MC_BKG0'])
fit.saveWorkspace()