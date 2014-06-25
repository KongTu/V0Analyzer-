import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.ana = cms.EDAnalyzer('V0Analyzer',
                          trackPtMin = cms.untracked.double(0.4),
                          simTrackPtMin = cms.untracked.double(0.4),
                          vertexSrc = cms.string('offlinePrimaryVertices'),
                          trackSrc = cms.InputTag('generalTracks'),
                          pfCandSrc = cms.InputTag('particleFlow'),
			  beamSpotSrc = cms.untracked.InputTag('offlineBeamSpot'),
                          doPFMatching = cms.untracked.bool(False),
                          doSimTrack = cms.untracked.bool(False),
                          doSimVertex = cms.untracked.bool(False),                          
                          useQuality = cms.untracked.bool(False),
                          qualityString = cms.untracked.string('highPurity'),
                          tpFakeSrc = cms.untracked.InputTag('mergedtruth','MergedTrackTruth'),
                          tpEffSrc = cms.untracked.InputTag('mergedtruth','MergedTrackTruth'),
                           
                          jetSrc = cms.string('ak5PFJets'),             
                          generalV0_ks = cms.InputTag('generalV0Candidates:Kshort'),
                          generalV0_la = cms.InputTag('generalV0Candidates:Lambda'),

                          doOfflinePrimaryVertex = cms.untracked.bool(True), 
                          doVertex = cms.untracked.bool(False),
                          doV0_Kshort = cms.untracked.bool(True),
                          doV0_Lambda = cms.untracked.bool(True),
                          doPFJet = cms.untracked.bool(True),
                          isMatchByHitsOrChi2 = cms.untracked.bool(True)
)

### standard includes
process.load("Configuration.StandardSequences.Digi_cff")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")

### conditions
#Global Tag:give access to different data info, needs to be changed for different data.
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#for MC data:
#process.GlobalTag.globaltag = 'STARTHI53_V17::All'

#globaltag for pPb data
process.GlobalTag.globaltag = 'GR_P_V43F::All'

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('Demo')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options   = cms.untracked.PSet( wantSummary =
cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1001_1_NXZ.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1022_1_URF.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1004_1_Fop.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1005_1_GTp.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1006_1_oRo.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1007_1_lcp.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1008_1_RDf.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1009_1_3ZV.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1010_1_CIp.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1011_1_RJF.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1012_1_fQ7.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1013_1_8cm.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1014_1_8kV.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1015_1_qAF.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1016_1_BAp.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1017_1_hsm.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1018_1_7A4.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1019_1_nuh.root',

#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_101_1_w46.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1020_1_5Mv.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1021_1_rwD.root',
#'root://xrootd1.cmsaf.mit.edu//store/user/vzhukova/HIJING_MB/HIJING_MB_RECO_v3/13a591fee6315e7fb1e99e6ba8e52eaa/reco_hijing_1022_1_URF.root'

'root://xrootd3.cmsaf.mit.edu//store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Reverse_v12/9d480a1aca4172fd4f159e4efe1592ae/pPb_HM_100_1_quK.root'
     

   )
)

process.TFileService = cms.Service("TFileService",fileName = cms.string("V0reco_pPb_jets.root"))


process.p = cms.Path(process.ana)