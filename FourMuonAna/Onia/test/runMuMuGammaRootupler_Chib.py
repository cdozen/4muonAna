import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')			#for 2016 
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v6', '')			#for 2017 Rereco
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Prompt_v4', '')			#for 2018 PromptReco
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Sep2018Rereco_v1', '')			#for 2018 Rereco

process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'file:/eos/uscms/store/user/lpcmuon/fourmuonMC/H0ToUps1SMuMu_m18p5_TuneCUEP8M1_13TeV-pythia8/BPHSkim-v10/180308_232212/0000/BPHSkim_1.root',
#'file:/eos/uscms/store/user/l1upgrades/Run2017/fourmuon/MuOnia/BPHSkim-v4-Run2017B-17Nov2017-v1/180314_061950/0000/BPHSkim_12.root'
#'file:/uscms_data/d3/huiwang/CMSSW_10_2_1/src/HeavyFlavorAnalysis/Onia2MuMu/test/BPHSkim_2018A_Rereco.root'		#A sample for 2018 Rereco

#'file:/eos/uscms/store/user/lpcbphy/hui/MuOnia/BPHSkim-2018A_Rereco-Run2018A-17Sep2018-v1/190118_180123/0000/BPHSkim_390.root')		#A sample for 2018 Rereco
#'file:/eos/uscms/store/user/lpcbphy/hui/MuOnia/BPHSkim-2018B_Rereco-Run2018B-17Sep2018-v1/190104_025327/0000/BPHSkim_426.root'),

#'file:/eos/uscms/store/user/lpcbphy/hui/MuOnia/BPHSkim-2018A_Rereco-Run2018A-17Sep2018-v1/190118_180123/0000/BPHSkim_390.root'		#A sample for 2018 Rereco
#'file:/eos/uscms/store/user/lpcbphy/hui/MuOnia/BPHSkim-2018B_Rereco-Run2018B-17Sep2018-v1/190104_025327/0000/BPHSkim_426.root'),
#'file:/eos/uscms/store/user/cdozen/FourMuon_Analysis/MuOnia/2016_v2/MuOnia/BPHSkim--Run2016D-07Aug17-v1/190620_113122/0000/BPHSkim_2016_95.root'  #new skim for 2016
'file:/eos/uscms/store/user/muahmad/FourMuon_Analysis/MuOnia/2017_v2/MuOnia/BPHSkim--Run2017B-17Nov2017-v1/190609_091355/0000/BPHSkim_2017_10.root') #new skim for 2017
#eventsToProcess = cms.untracked.VEventRange('317641:1352331850-317641:1352331850')
)


process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Rootuple.root'),
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load('FourMuonAna.Onia.MuMuGammaRootupler_cfi')
process.p = cms.Path(process.rootuple)

process.rootuple.upsilon_mass = cms.double(9.4603)
process.rootuple.triggerCuts = cms.uint32(36)
process.rootuple.isMC = cms.bool(False)                 # is mc?
process.rootuple.onia_mass_cuts = cms.vdouble(0.1,500)    # you may need to adjust this
#process.rootuple.SecondSource.fileNames = cms.untracked.vstring(
#    'file:/eos/cms/store/user/jblee/SpyFEDemulated234824.root'.
#     'file:/eos/uscms/store/user/l1upgrades/Run2017/fourmuon/ZeroBias/BPHSkim-v4-Run2017F-17Nov2017-v1/180314_055941/0000/BPHSkim_425.root',
#	  'file:/eos/uscms/store/user/l1upgrades/Run2017/fourmuon/ZeroBias/BPHSkim-v4-Run2017F-17Nov2017-v1/180314_055941/0000/BPHSkim_426.root'
#	  )



