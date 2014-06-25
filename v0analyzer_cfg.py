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
                           
                          jetSrc = cms.InputTag('ak5PFJets'),             
                          generalV0_ks = cms.InputTag('generalV0Candidates:Kshort'),
                          generalV0_la = cms.InputTag('generalV0Candidates:Lambda'),
                          genParticleSrc = cms.InputTag('genParticles'),
                          genJetSrc = cms.InputTag('iterativeCone5GenJets'),
                          caloJetSrc = cms.InputTag('iterativeCone5CaloJets'),
                          doGenJet = cms.untracked.bool(True), 
               		        doCaloJet = cms.untracked.bool(True),
                          doGenParticle = cms.untracked.bool(True),
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

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options   = cms.untracked.PSet( wantSummary =
cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

#'root://xrootd3.cmsaf.mit.edu//store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_TrkHM_Gplus_Reverse_v12/9d480a1aca4172fd4f159e4efe1592ae/pPb_HM_100_1_quK.root'
 'root://xrootd.unl.edu//store/user/davidlw/PAHighPt/PA2013_FlowCorr_PromptReco_Jet_Gplus_v12/d2a71c3e10bdac33a8236fd290776517/pPb_Jet_100_1_LKn.root'    

   )
)

process.TFileService = cms.Service("TFileService",fileName = cms.string("V0reco_pPb_jetTrig.root"))


process.p = cms.Path(process.ana)
