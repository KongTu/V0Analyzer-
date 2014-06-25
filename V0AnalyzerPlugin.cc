// -*- C++ -*-
//
// Package:    V0Analyzer
// Class:      V0Analyzer
// 
/**\class V0Analyzer V0Analyzer.cc V0Analyzer/V0Analyzer/src/V0AnalyzerPlugin.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhoudunming Tu,,,
//         Created:  Mon Feb  3 01:26:00 CET 2014
// $Id$
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <map>


#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TNtuple.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

//////////////////////////////////////////////
// CMSSW user include files
#include "DataFormats/Common/interface/DetSetAlgorithm.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
//#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerLayerIdAccessor.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"

#include "L1Trigger/GlobalTrigger/interface/L1GlobalTrigger.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// Heavyion
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"

// Particle Flow
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"

// Vertex significance
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

//Jet:
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

// Root include files
#include "TTree.h"
//
// Track Matching and fake rate calculations     
//#include "RiceHIG/V0Analysis/interface/V0Validator.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
//
// class declaration
//
using namespace std;
using namespace edm;
using namespace reco;
//

#define PI 3.14159265358979

#define MAXCOUNTS 50000
#define MAXVTX 100
#define MAXQUAL 5
#define MAXNUMBER 1000
//
// declare struct //
//

struct V0Candidate{

  // event information
  int nRun;
  int nEv;
  int nLumi;
  int nBX;
  
  int N; // multiplicity variable 
  int N1;
  int nJets;
  int nTracks;
  int EtaPtCutnTracks;
  int nGen;
  int nGenJet;
  int nCaloJet;

  // Vertex information
  int nv;
  float vx[MAXVTX];
  float vy[MAXVTX];
  float vz[MAXVTX];
  float vxError[MAXVTX];
  float vyError[MAXVTX];
  float vzError[MAXVTX];
  int nDaughter[MAXVTX];

  float bestvz;
  float bestvx;
  float bestvy;
  float bestvxError;
  float bestvyError;
  float bestvzError;

  //// Multiple vtx information
  int nVtx;

  int nTrkVtx[MAXVTX];
  float normChi2Vtx[MAXVTX];
  float sumPtVtx[MAXVTX];
  //int nTrkVtxHard[MAXVTX];
  int maxVtx;
  //int maxVtxHard;

  float xVtx[MAXVTX];
  float yVtx[MAXVTX];
  float zVtx[MAXVTX];
  float xVtxErr[MAXVTX];
  float yVtxErr[MAXVTX];
  float zVtxErr[MAXVTX];

  float vtxDist2D[MAXVTX];
  float vtxDist2DErr[MAXVTX];
  float vtxDist2DSig[MAXVTX];
  float vtxDist3D[MAXVTX];
  float vtxDist3DErr[MAXVTX];
  float vtxDist3DSig[MAXVTX];


  //track information
  int nTrk;

  unsigned int trkVtxIndex[MAXNUMBER];

//--------------------------------------
//Kshort candidate
//--------------------------------------

  float ks_eta[MAXNUMBER];
  float ks_phi[MAXNUMBER];
  float ks_pt[MAXNUMBER];
  float ks_mass[MAXNUMBER];
  float ks_p[MAXNUMBER];
  float ks_px[MAXNUMBER];
  float ks_py[MAXNUMBER];
  float ks_pz[MAXNUMBER];
  float ks_charge[MAXNUMBER];
  float ks_vtxChi2[MAXNUMBER];
  float ks_agl[MAXNUMBER];
  
  float ks_qT[MAXNUMBER];
  float ks_alpha[MAXNUMBER];

  float ks_dl[MAXNUMBER];
  float ks_dlerror[MAXNUMBER];
  float ks_dlos[MAXNUMBER];

  float ks_vx[MAXNUMBER];
  float ks_vy[MAXNUMBER];
  float ks_vz[MAXNUMBER];
  float ks_vxError[MAXNUMBER];
  float ks_vyError[MAXNUMBER];
  float ks_vzError[MAXNUMBER];

  float ks_dau1_eta[MAXNUMBER];
  float ks_dau1_phi[MAXNUMBER];
  float ks_dau1_pt[MAXNUMBER];
  float ks_dau1_px[MAXNUMBER];
  float ks_dau1_py[MAXNUMBER];
  float ks_dau1_pz[MAXNUMBER];
  float ks_dau1_charge[MAXNUMBER];
  float ks_dau1_p[MAXNUMBER];
  
  float ks_dau1_nhit[MAXNUMBER];
  float ks_dau1_dzbest[MAXNUMBER];
  float ks_dau1_dxybest[MAXNUMBER];
  float ks_dau1_dzerror[MAXNUMBER];
  float ks_dau1_dxyerror[MAXNUMBER];
  float ks_dau1_dzos[MAXNUMBER];
  float ks_dau1_dxyos[MAXNUMBER];
  float ks_dau1_dedx[MAXNUMBER];
  
  float ks_dau2_eta[MAXNUMBER];
  float ks_dau2_phi[MAXNUMBER];
  float ks_dau2_pt[MAXNUMBER];
  float ks_dau2_px[MAXNUMBER];
  float ks_dau2_py[MAXNUMBER];
  float ks_dau2_pz[MAXNUMBER];
  float ks_dau2_charge[MAXNUMBER];
  float ks_dau2_p[MAXNUMBER];
 
  float ks_dau2_nhit[MAXNUMBER];
  float ks_dau2_dzbest[MAXNUMBER];
  float ks_dau2_dxybest[MAXNUMBER];
  float ks_dau2_dzerror[MAXNUMBER];
  float ks_dau2_dxyerror[MAXNUMBER];
  float ks_dau2_dzos[MAXNUMBER];
  float ks_dau2_dxyos[MAXNUMBER];
  float ks_dau2_dedx[MAXNUMBER];

//---------------------------------------
//Lambda candidate;
//---------------------------------------

  float la_eta[MAXNUMBER];
  float la_phi[MAXNUMBER];
  float la_pt[MAXNUMBER];
  float la_mass[MAXNUMBER];
  float la_p[MAXNUMBER];
  float la_px[MAXNUMBER];
  float la_py[MAXNUMBER];
  float la_pz[MAXNUMBER];
  float la_charge[MAXNUMBER];
  float la_vtxChi2[MAXNUMBER];
  float la_agl[MAXNUMBER];
  
  float la_qT[MAXNUMBER];
  float la_alpha[MAXNUMBER];

  float la_dl[MAXNUMBER];
  float la_dlerror[MAXNUMBER];
  float la_dlos[MAXNUMBER];

  float la_vx[MAXNUMBER];
  float la_vy[MAXNUMBER];
  float la_vz[MAXNUMBER];

  float la_vxError[MAXNUMBER];
  float la_vyError[MAXNUMBER];
  float la_vzError[MAXNUMBER];

  float la_dau1_eta[MAXNUMBER];
  float la_dau1_phi[MAXNUMBER];
  float la_dau1_pt[MAXNUMBER];
  float la_dau1_px[MAXNUMBER];
  float la_dau1_py[MAXNUMBER];
  float la_dau1_pz[MAXNUMBER];
  float la_dau1_charge[MAXNUMBER];
  float la_dau1_p[MAXNUMBER];
  
  float la_dau1_nhit[MAXNUMBER];
  float la_dau1_dzbest[MAXNUMBER];
  float la_dau1_dxybest[MAXNUMBER];
  float la_dau1_dzerror[MAXNUMBER];
  float la_dau1_dxyerror[MAXNUMBER];
  float la_dau1_dzos[MAXNUMBER];
  float la_dau1_dxyos[MAXNUMBER];
  float la_dau1_dedx[MAXNUMBER];
  
  float la_dau2_eta[MAXNUMBER];
  float la_dau2_phi[MAXNUMBER];
  float la_dau2_pt[MAXNUMBER];
  float la_dau2_px[MAXNUMBER];
  float la_dau2_py[MAXNUMBER];
  float la_dau2_pz[MAXNUMBER];
  float la_dau2_charge[MAXNUMBER];
  float la_dau2_p[MAXNUMBER];
 
  float la_dau2_nhit[MAXNUMBER];
  float la_dau2_dzbest[MAXNUMBER];
  float la_dau2_dxybest[MAXNUMBER];
  float la_dau2_dzerror[MAXNUMBER];
  float la_dau2_dxyerror[MAXNUMBER];
  float la_dau2_dzos[MAXNUMBER];
  float la_dau2_dxyos[MAXNUMBER];
  float la_dau2_dedx[MAXNUMBER];

  //--------------------
  //PFJets candidates;
  //--------------------
   
  int jet_numberOfDaughters[MAXNUMBER]; 
  
  float jet_vx[MAXNUMBER];
  float jet_vy[MAXNUMBER];
  float jet_vz[MAXNUMBER];
  float jet_pt[MAXNUMBER];
  float jet_px[MAXNUMBER];
  float jet_py[MAXNUMBER];
  float jet_pz[MAXNUMBER];
  float jet_eta[MAXNUMBER];
  float jet_phi[MAXNUMBER];
  float jet_p[MAXNUMBER];
  
  float jet_dau_pt[MAXNUMBER][50];
  float jet_dau_px[MAXNUMBER][50];
  float jet_dau_py[MAXNUMBER][50];
  float jet_dau_pz[MAXNUMBER][50];
  float jet_dau_p[MAXNUMBER][50];
  float jet_dau_eta[MAXNUMBER][50];
  float jet_dau_phi[MAXNUMBER][50];

  float jet_chargedEmEnergy[MAXNUMBER];
  float jet_maxDistance[MAXNUMBER];
  float jet_jetArea[MAXNUMBER];

  //---------------------
  //CaloJet:
  //---------------------

  float calojet_pt[MAXNUMBER];
  float calojet_eta[MAXNUMBER];
  float calojet_phi[MAXNUMBER];
  float calojet_px[MAXNUMBER];
  float calojet_py[MAXNUMBER];
  float calojet_pz[MAXNUMBER];

  //---------------------
  //General Tracks
  //---------------------
  float track_pt[MAXNUMBER];
  float track_eta[MAXNUMBER];
  float track_phi[MAXNUMBER];
  int track_charge[MAXNUMBER];
  float track_p[MAXNUMBER];

  //---------------------
  //GenParticle
  //---------------------
  
  int pdg[MAXNUMBER];
  float genP_pt[MAXNUMBER];
  float genP_eta[MAXNUMBER];
  float genP_phi[MAXNUMBER];
  int genP_mid[MAXNUMBER];


  //----------------------
  //GenJet
  //----------------------

  float genpt[MAXNUMBER];
  float geneta[MAXNUMBER];
  float genphi[MAXNUMBER];
  float geny[MAXNUMBER];


};



class V0AnalyzerPlugin : public edm::EDAnalyzer {
   public:
      explicit V0AnalyzerPlugin(const edm::ParameterSet&);
      ~V0AnalyzerPlugin();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      void V0_Kshort(const edm::Event& iEvent);

      void V0_Lambda(const edm::Event& iEvent);

      void fillVertices(const edm::Event& iEvent);

      void OfflinePrimaryVertex(const edm::Event& iEvent);

      void PF_Jet(const edm::Event& iEvent, const edm::EventSetup& iSetup);

      void Calo_Jet(const edm::Event& iEvent, const edm::EventSetup& iSetup);

      void GeneralTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup);

      void GenParticle(const edm::Event& iEvent, const edm::EventSetup& iSetup);

      void GenJet(const edm::Event& iEvent, const edm::EventSetup& iSetup);


      // ----------member data ---------------------------
      edm::InputTag trackSrc_;
      edm::InputTag simVertexSrc_;
      edm::InputTag generalV0_ks_;
      edm::InputTag generalV0_la_;
      edm::InputTag genParticleSrc_;
      edm::InputTag caloJetSrc_;
      edm::InputTag genJetSrc_;


      std::string vertexSrc_;
      std::string jetSrc_;

      bool doSimVertex_;
      bool doVertex_;
      bool doOfflinePrimaryVertex_;
      bool doV0_Kshort_;
      bool doV0_Lambda_;
      bool doPFJet_;
      bool doGeneralTracks_;
      bool doGenParticle_;
      bool doCaloJet_;
      bool doGenJet_;



      // Root Object

      TTree* v0Tree_;
      TTree* v0Tree_VTX_;
      TTree* v0_Kshort_;
      TTree* v0_Lambda_;
      TTree* PF_Jet_;
      TTree* GeneralTracks_;
      TTree* GenParticle_;
      TTree* Calo_Jet_;
      TTree* GenJet_;


      V0Candidate event_;


};



V0AnalyzerPlugin::V0AnalyzerPlugin(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  trackSrc_ = iConfig.getParameter<edm::InputTag>("trackSrc");
  vertexSrc_ = iConfig.getParameter<std::string>("vertexSrc");
  simVertexSrc_ =  iConfig.getUntrackedParameter<edm::InputTag>("tpVtxSrc",edm::InputTag("mergedtruth","MergedTrackTruth"));
  generalV0_ks_ = iConfig.getParameter<edm::InputTag>("generalV0_ks");
  generalV0_la_ = iConfig.getParameter<edm::InputTag>("generalV0_la");
  jetSrc_ = iConfig.getParameter<std::string>("jetSrc");
  genParticleSrc_ = iConfig.getParameter<edm::InputTag>("genParticleSrc");
  caloJetSrc_ = iConfig.getParameter<edm::InputTag>( "caloJetSrc" );
  genJetSrc_ = iConfig.getParameter<edm::InputTag>("genJetSrc");


  doPFJet_ = iConfig.getUntrackedParameter<bool> ("doPFJet", true);
  doSimVertex_ = iConfig.getUntrackedParameter<bool>  ("doSimVertex",false);
  doVertex_ = iConfig.getUntrackedParameter<bool> ("doVertex",false);
  doOfflinePrimaryVertex_ = iConfig.getUntrackedParameter<bool>("doOfflinePrimaryVertex",true);
  doV0_Kshort_ = iConfig.getUntrackedParameter<bool>("doV0_Kshort",true);
  doV0_Lambda_ = iConfig.getUntrackedParameter<bool>("doV0_Lambda",true);
  doGeneralTracks_ = iConfig.getUntrackedParameter<bool>("doGeneralTracks",true);
  doGenParticle_ = iConfig.getUntrackedParameter<bool>("doGenParticle",false);
  doCaloJet_ = iConfig.getUntrackedParameter<bool>("doCaloJet",false);
  doGenJet_ = iConfig.getUntrackedParameter<bool>("doGenJet",false);

}


V0AnalyzerPlugin::~V0AnalyzerPlugin()
{
 
   // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

double Angle(double px, double py, double pz, double px_d, double py_d, double pz_d)
{
       double pd = px*px_d+py*py_d+pz*pz_d;
       double p = sqrt(px*px+py*py+pz*pz);
       double d = sqrt(px_d*px_d+py_d*py_d+pz_d*pz_d);
       double pd_module = p*d;
       double temp = (pd)/(pd_module);
       double angle = acos(temp);
       return angle;

}

//
// member functions
//

// ------------ method called for each event  ------------
void
V0AnalyzerPlugin::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   event_.nEv = (int)iEvent.id().event();
   event_.nRun = (int)iEvent.id().run();
   event_.nLumi = (int)iEvent.luminosityBlock();
   event_.nBX = (int)iEvent.bunchCrossing();

   event_.N = 0;
   event_.N1 = 0;
   event_.nv = 0;
   event_.nTrk = 0;
   event_.nJets = 0;
   event_.nTracks = 0;
   event_.EtaPtCutnTracks = 0;
   event_.nGen = 0;
   event_.nCaloJet = 0;
   event_.nGenJet = 0;
   

   if(doV0_Kshort_) V0_Kshort(iEvent);

   if(doV0_Lambda_) V0_Lambda(iEvent);
   //Fill greatest vertices(primary vertices);

   if(doVertex_) fillVertices(iEvent);

   //Fill offlinePrimaryVertex

   if(doOfflinePrimaryVertex_) OfflinePrimaryVertex(iEvent);

   if (doPFJet_) PF_Jet( iEvent, iSetup );

   if (doGeneralTracks_) GeneralTracks( iEvent, iSetup );

   if (doGenParticle_) GenParticle( iEvent, iSetup );

   if (doCaloJet_) Calo_Jet(iEvent, iSetup);

   if (doGenJet_) GenJet( iEvent, iSetup );
  
//Fill in the trees:

   v0_Kshort_->Fill();
   v0_Lambda_->Fill();
   v0Tree_->Fill();
   v0Tree_VTX_->Fill();
   PF_Jet_->Fill();
   GeneralTracks_->Fill();
   GenParticle_->Fill();
   Calo_Jet_->Fill();
   GenJet_->Fill();

}


void 
V0AnalyzerPlugin::GenJet(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace edm;
  using namespace std;
  using namespace reco;
  using std::vector;

  edm::Handle<vector<reco::GenJet> > genjets;
  iEvent.getByLabel(genJetSrc_, genjets);

   for(unsigned int igen = 0; igen < genjets->size(); ++igen){
     
     const reco::GenJet & genjet = (*genjets)[igen];

     float genjet_pt = genjet.pt();

     // threshold to reduce size of output in minbias PbPb
     if( genjet_pt > 20 ){

       event_.genpt[event_.nGenJet] = genjet_pt;
       event_.geneta[event_.nGenJet] = genjet.eta();
       event_.genphi[event_.nGenJet] = genjet.phi();
       event_.geny[event_.nGenJet] = genjet.eta();

      }

      event_.nGenJet++;
    }
}

void 
V0AnalyzerPlugin::Calo_Jet(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace edm;
  using namespace std;
  using namespace reco;
  using std::vector;
 
  edm::Handle<vector<reco::CaloJet> > Calojets;
  iEvent.getByLabel( caloJetSrc_ , Calojets);

  for (unsigned it = 0; it < Calojets->size(); it++){

    const reco::CaloJet & calojets = (*Calojets)[it];

    event_.calojet_pt[event_.nCaloJet] = calojets.pt();
    event_.calojet_eta[event_.nCaloJet] = calojets.eta();
    event_.calojet_phi[event_.nCaloJet] = calojets.phi();
    event_.calojet_px[event_.nCaloJet] = calojets.px();
    event_.calojet_py[event_.nCaloJet] = calojets.py();
    event_.calojet_pz[event_.nCaloJet] = calojets.pz();

    event_.nCaloJet++;

  }

}

void
V0AnalyzerPlugin::GenParticle(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace edm;
  using namespace std;
  using namespace reco;
  using std::vector;

  edm::Handle<reco::GenParticleCollection> genParticleCollection;
  iEvent.getByLabel(genParticleSrc_, genParticleCollection);

  for(unsigned it=0; it<genParticleCollection->size(); ++it) {

    const reco::GenParticle & genCand = (*genParticleCollection)[it];
    int id = genCand.pdgId();
    int status = genCand.status();

  if ( status == 1 ){

        event_.pdg[event_.nGen] = genCand.pdgId();
        event_.genP_pt[event_.nGen] = genCand.pt();
        event_.genP_eta[event_.nGen] = genCand.eta();
        event_.genP_phi[event_.nGen] = genCand.phi();

//Finding mother:

    int mid = 0;
    if( TMath::Abs(id) == 3122 ){

      if(genCand.numberOfMothers()==1){
        const reco::Candidate * mom = genCand.mother();
        mid = mom->pdgId();
        if(mom->numberOfMothers()==1){
          const reco::Candidate * mom1 = mom->mother();
          mid = mom1->pdgId();
        }
      }

      event_.genP_mid[event_.nGen] = mid;

    }
 
      event_.nGen++;  
   }
  }

}

void 
V0AnalyzerPlugin::GeneralTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace edm;
  using namespace std;
  using namespace reco;
  using std::vector;

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vertexSrc_,vertices);
  double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
  double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
  const reco::Vertex & vtx = (*vertices)[0];
  bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
  bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(trackSrc_, tracks);

  for(unsigned it = 0; it < tracks->size(); it++){

     const reco::Track & trk = (*tracks)[it];

     math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
        
        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt()>0.10) continue;
        if(fabs(dzvtx/dzerror) > 3) continue;
        if(fabs(dxyvtx/dxyerror) > 3) continue;
        
        event_.track_eta[event_.nTracks] = trk.eta();
        event_.track_phi[event_.nTracks] = trk.phi();
        event_.track_pt[event_.nTracks] = trk.pt();
        event_.track_charge[event_.nTracks] = trk.charge();
        event_.track_p[event_.nTracks] = trk.p();

  event_.nTracks++;

  //ntrack selection:

        if ( fabs(trk.eta()) > 2.4 || trk.pt() < 0.4  ) continue;

        event_.EtaPtCutnTracks++;
        
  }
}

void 
V0AnalyzerPlugin::PF_Jet(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  using namespace edm;
  using namespace std;
  using namespace reco;
  using std::vector;

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByLabel( jetSrc_ , jets);
   // jetSrc = "ak5PFJets" inputTag (string);

  event_.nJets = 0;

  for(unsigned it = 0; it < jets->size(); ++it){

    const reco::PFJet & pfj = (*jets)[it];
  
    vector<const reco::Candidate*> myDaughters;

//can't include any daughter information for jet in this jet trigger data.
    
/*
        for ( unsigned int j = 0; j < pfj.numberOfDaughters(); j++){

            //const reco::Candidate * d1 = pfj.daughter(i);
            myDaughters.push_back( pfj.daughter( j ) );
            
            event_.jet_dau_pt[event_.nJets][j] = myDaughters[j]->pt();
            event_.jet_dau_pt[event_.nJets][j] = myDaughters[j]->pt();
            event_.jet_dau_px[event_.nJets][j] = myDaughters[j]->px();
            event_.jet_dau_py[event_.nJets][j] = myDaughters[j]->py();
            event_.jet_dau_pz[event_.nJets][j] = myDaughters[j]->pz();
            event_.jet_dau_p[event_.nJets][j] = myDaughters[j]->p();
            event_.jet_dau_eta[event_.nJets][j] = myDaughters[j]->eta();
            event_.jet_dau_phi[event_.nJets][j] = myDaughters[j]->phi();
          }
    
*/
    


    //event_.jet_numberOfDaughters[event_.nJets] = pfj.numberOfDaughters();
    event_.jet_chargedEmEnergy[event_.nJets] = pfj.chargedEmEnergy();
    //event_.jet_maxDistance[event_.nJets] = pfj.maxDistance();

    //vertex
    event_.jet_vx[event_.nJets] = pfj.vx();
    event_.jet_vy[event_.nJets] = pfj.vy();
    event_.jet_vz[event_.nJets] = pfj.vz();

    //momentum
    event_.jet_pt[event_.nJets] = pfj.pt();
    event_.jet_px[event_.nJets] = pfj.px();
    event_.jet_py[event_.nJets] = pfj.py();
    event_.jet_pz[event_.nJets] = pfj.pz();
   
    event_.jet_p[event_.nJets] = pfj.p();
    event_.jet_eta[event_.nJets] = pfj.eta();
    event_.jet_phi[event_.nJets] = pfj.phi();
   
    event_.nJets++;  


  }

}


void 
V0AnalyzerPlugin::OfflinePrimaryVertex(const edm::Event& iEvent){

  using namespace edm;
  using std::vector;

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vertexSrc_,vertices);

  const reco::Vertex & vtx = (*vertices)[0];
  
  event_.bestvx = vtx.x();
  event_.bestvy = vtx.y();
  event_.bestvz = vtx.z();
  event_.bestvxError = vtx.xError();
  event_.bestvyError = vtx.yError();
  event_.bestvzError = vtx.zError();

}

void
V0AnalyzerPlugin::V0_Kshort(const edm::Event& iEvent){

    using namespace edm;
    using std::vector;

    //Primary Vertex//
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByLabel(vertexSrc_,vertices);

    const reco::Vertex & vtx = (*vertices)[0];
      
      double bestvx = vtx.x();
      double bestvy = vtx.y();
      double bestvz = vtx.z();
      double bestvxError = vtx.xError();
      double bestvyError = vtx.yError();
      double bestvzError = vtx.zError();

    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_ks;
    iEvent.getByLabel( generalV0_ks_ ,v0candidates_ks);
    if(!v0candidates_ks.isValid()) return;
/*
    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle;
    iEvent.getByLabel("dedxTruncated40", dEdxHandle);
*/
     event_.N = 0;

        for(unsigned it = 0; it < v0candidates_ks->size(); ++it){

          const reco::VertexCompositeCandidate & trk = (*v0candidates_ks)[it];

          const reco::Candidate * d1 = trk.daughter(0);
          const reco::Candidate * d2 = trk.daughter(1);

          auto dau1 = d1->get<reco::TrackRef>();
          auto dau2 = d2->get<reco::TrackRef>();

          //defining variables for Kshort and daughters

          event_.ks_eta[event_.N] = trk.eta();
          event_.ks_phi[event_.N] = trk.phi();
          event_.ks_pt[event_.N] = trk.pt();
          event_.ks_mass[event_.N] = trk.mass();
          event_.ks_p[event_.N] = trk.p();
          event_.ks_px[event_.N] = trk.px();
          event_.ks_py[event_.N] = trk.py();
          event_.ks_pz[event_.N] = trk.pz();
          event_.ks_charge[event_.N] = trk.charge();
          event_.ks_vtxChi2[event_.N] = trk.vertexChi2();
          event_.ks_vx[event_.N] = trk.vx();
          event_.ks_vy[event_.N] = trk.vy();
          event_.ks_vz[event_.N] = trk.vz();

          event_.ks_dau1_eta[event_.N] = dau1->eta();
          event_.ks_dau1_phi[event_.N] = dau1->phi();
          event_.ks_dau1_pt[event_.N] = dau1->pt();
          event_.ks_dau1_px[event_.N] = dau1->px();
          event_.ks_dau1_py[event_.N] = dau1->py();
          event_.ks_dau1_pz[event_.N] = dau1->pz();
          event_.ks_dau1_charge[event_.N] = dau1->charge();
          event_.ks_dau1_nhit[event_.N] = dau1->numberOfValidHits();
          event_.ks_dau1_p[event_.N] = dau1->p();

          event_.ks_dau2_eta[event_.N] = dau2->eta();
          event_.ks_dau2_phi[event_.N] = dau2->phi();
          event_.ks_dau2_pt[event_.N] = dau2->pt();
          event_.ks_dau2_px[event_.N] = dau2->px();
          event_.ks_dau2_py[event_.N] = dau2->py();
          event_.ks_dau2_pz[event_.N] = dau2->pz();
          event_.ks_dau2_charge[event_.N] = dau2->charge();
          event_.ks_dau2_nhit[event_.N] = dau2->numberOfValidHits();
          event_.ks_dau2_p[event_.N] = dau2->p();

          //PAngle 

          TVector3 ptosvec(event_.ks_vx[event_.N]-bestvx,event_.ks_vy[event_.N]-bestvy,event_.ks_vz[event_.N]-bestvz);
          TVector3 secvec(event_.ks_px[event_.N],event_.ks_py[event_.N],event_.ks_pz[event_.N]);
                
          double agl = cos(secvec.Angle(ptosvec));

          event_.ks_agl[event_.N] = agl;

          //Decay length

          typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
          typedef ROOT::Math::SVector<double, 3> SVector3;
          
          SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
          SVector3 distanceVector(event_.ks_vx[event_.N]-bestvx,event_.ks_vy[event_.N]-bestvy,event_.ks_vz[event_.N]-bestvz);
          
          double dl = ROOT::Math::Mag(distanceVector);
          double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;

          double dlos = dl/dlerror;

          event_.ks_dl[event_.N] = dl;
          event_.ks_dlerror[event_.N] = dlerror;
          event_.ks_dlos[event_.N] = dlos;

          //DCA
          math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
          
          double dzbest1 = dau1->dz(bestvtx);
          double dxybest1 = dau1->dxy(bestvtx);
          double dzerror1 = sqrt(dau1->dzError()*dau1->dzError()+bestvzError*bestvzError);
          double dxyerror1 = sqrt(dau1->d0Error()*dau1->d0Error()+bestvxError*bestvyError);

          double dzos1 = dzbest1/dzerror1;
          double dxyos1 = dxybest1/dxyerror1;
          
          double dzbest2 = dau2->dz(bestvtx);
          double dxybest2 = dau2->dxy(bestvtx);
          double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
          double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
          
          double dzos2 = dzbest2/dzerror2;
          double dxyos2 = dxybest2/dxyerror2;

          event_.ks_dau1_dzbest[event_.N] = dzbest1;
          event_.ks_dau1_dxybest[event_.N] = dxybest1;
          event_.ks_dau1_dzerror[event_.N] = dzerror1;
          event_.ks_dau1_dxyerror[event_.N] = dxyerror1;
          event_.ks_dau1_dzos[event_.N] = dzos1;
          event_.ks_dau1_dxyos[event_.N] = dxyos1;

          event_.ks_dau2_dzbest[event_.N] = dzbest2;
          event_.ks_dau2_dxybest[event_.N] = dxybest2;
          event_.ks_dau2_dzerror[event_.N] = dzerror2;
          event_.ks_dau2_dxyerror[event_.N] = dxyerror2;
          event_.ks_dau2_dzos[event_.N] = dzos2;
          event_.ks_dau2_dxyos[event_.N] = dxyos2;

          //qT AND alpha;

           double A2 = Angle(event_.ks_px[event_.N],event_.ks_py[event_.N],event_.ks_pz[event_.N],event_.ks_dau1_px[event_.N],event_.ks_dau1_py[event_.N],event_.ks_dau1_pz[event_.N]);
           //stop here!!02/21/2014:03:56pm
           double P2 = event_.ks_dau1_p[event_.N];
           double qT2 = P2*sin(A2);
           double qL1 = event_.ks_dau1_p[event_.N]*cos(A2);
           double qL2 = event_.ks_dau2_p[event_.N]*cos(A2);

           double alpha = (qL1-qL2)/(qL1+qL2);

           event_.ks_qT[event_.N] = qT2;
           event_.ks_alpha[event_.N] = alpha;

          //dEdx
/*
           if(dEdxHandle->size()){
                const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle.product();
                event_.ks_dau1_dedx[event_.N] = dEdxTrack[dau1].dEdx();
                event_.ks_dau2_dedx[event_.N] = dEdxTrack[dau2].dEdx();
            }
*/
          event_.N++;//inside of Kshort loop
        }

}

void 
V0AnalyzerPlugin::V0_Lambda(const edm::Event& iEvent){

    using namespace edm;
    using std::vector;

    //Primary Vertex//
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByLabel(vertexSrc_,vertices);

    const reco::Vertex & vtx = (*vertices)[0];
      
      double bestvx = vtx.x();
      double bestvy = vtx.y();
      double bestvz = vtx.z();
      double bestvxError = vtx.xError();
      double bestvyError = vtx.yError();
      double bestvzError = vtx.zError();

    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_la;
    iEvent.getByLabel( generalV0_la_ ,v0candidates_la);
    if(!v0candidates_la.isValid()) return;
/*
    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle;
    iEvent.getByLabel("dedxTruncated40", dEdxHandle);
*/

     //number of V0(kshort)
     event_.N1 = 0;

     for(unsigned it = 0; it < v0candidates_la->size(); ++it){

          const reco::VertexCompositeCandidate & trk = (*v0candidates_la)[it];

          const reco::Candidate * d1 = trk.daughter(0);
          const reco::Candidate * d2 = trk.daughter(1);

          auto dau1 = d1->get<reco::TrackRef>();
          auto dau2 = d2->get<reco::TrackRef>();

          //defining variables for Lambda and daughters

          event_.la_eta[event_.N1] = trk.eta();
          event_.la_phi[event_.N1] = trk.phi();
          event_.la_pt[event_.N1] = trk.pt();
          event_.la_mass[event_.N1] = trk.mass();
          event_.la_p[event_.N] = trk.p();
          event_.la_px[event_.N1] = trk.px();
          event_.la_py[event_.N1] = trk.py();
          event_.la_pz[event_.N1] = trk.pz();
          event_.la_charge[event_.N1] = trk.charge();
          event_.la_vtxChi2[event_.N1] = trk.vertexChi2();
          event_.la_vx[event_.N1] = trk.vx();
          event_.la_vy[event_.N1] = trk.vy();
          event_.la_vz[event_.N1] = trk.vz();

          event_.la_dau1_eta[event_.N1] = dau1->eta();
          event_.la_dau1_phi[event_.N1] = dau1->phi();
          event_.la_dau1_pt[event_.N1] = dau1->pt();
          event_.la_dau1_px[event_.N1] = dau1->px();
          event_.la_dau1_py[event_.N1] = dau1->py();
          event_.la_dau1_pz[event_.N1] = dau1->pz();
          event_.la_dau1_charge[event_.N1] = dau1->charge();
          event_.la_dau1_nhit[event_.N1] = dau1->numberOfValidHits();
          event_.la_dau1_p[event_.N1] = dau1->p();

          event_.la_dau2_eta[event_.N1] = dau2->eta();
          event_.la_dau2_phi[event_.N1] = dau2->phi();
          event_.la_dau2_pt[event_.N1] = dau2->pt();
          event_.la_dau2_px[event_.N1] = dau2->px();
          event_.la_dau2_py[event_.N1] = dau2->py();
          event_.la_dau2_pz[event_.N1] = dau2->pz();
          event_.la_dau2_charge[event_.N1] = dau2->charge();
          event_.la_dau2_nhit[event_.N1] = dau2->numberOfValidHits();
          event_.la_dau2_p[event_.N1] = dau2->p();

          //PAngle 

          TVector3 ptosvec(event_.la_vx[event_.N1]-bestvx,event_.la_vy[event_.N1]-bestvy,event_.la_vz[event_.N1]-bestvz);
          TVector3 secvec(event_.la_px[event_.N1],event_.la_py[event_.N1],event_.la_pz[event_.N1]);
                
          double agl = cos(secvec.Angle(ptosvec));

          event_.la_agl[event_.N1] = agl;

          //Decay length

          typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
          typedef ROOT::Math::SVector<double, 3> SVector3;
          
          SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
          SVector3 distanceVector(event_.la_vx[event_.N1]-bestvx,event_.la_vy[event_.N1]-bestvy,event_.la_vz[event_.N1]-bestvz);
          
          double dl = ROOT::Math::Mag(distanceVector);
          double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;

          double dlos = dl/dlerror;

          event_.la_dl[event_.N1] = dl;
          event_.la_dlerror[event_.N1] = dlerror;
          event_.la_dlos[event_.N1] = dlos;

          //DCA
          math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
          
          double dzbest1 = dau1->dz(bestvtx);
          double dxybest1 = dau1->dxy(bestvtx);
          double dzerror1 = sqrt(dau1->dzError()*dau1->dzError()+bestvzError*bestvzError);
          double dxyerror1 = sqrt(dau1->d0Error()*dau1->d0Error()+bestvxError*bestvyError);

          double dzos1 = dzbest1/dzerror1;
          double dxyos1 = dxybest1/dxyerror1;
          
          double dzbest2 = dau2->dz(bestvtx);
          double dxybest2 = dau2->dxy(bestvtx);
          double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
          double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
          
          double dzos2 = dzbest2/dzerror2;
          double dxyos2 = dxybest2/dxyerror2;

          event_.la_dau1_dzbest[event_.N1] = dzbest1;
          event_.la_dau1_dxybest[event_.N1] = dxybest1;
          event_.la_dau1_dzerror[event_.N1] = dzerror1;
          event_.la_dau1_dxyerror[event_.N1] = dxyerror1;
          event_.la_dau1_dzos[event_.N1] = dzos1;
          event_.la_dau1_dxyos[event_.N1] = dxyos1;

          event_.la_dau2_dzbest[event_.N1] = dzbest2;
          event_.la_dau2_dxybest[event_.N1] = dxybest2;
          event_.la_dau2_dzerror[event_.N1] = dzerror2;
          event_.la_dau2_dxyerror[event_.N1] = dxyerror2;
          event_.la_dau2_dzos[event_.N1] = dzos2;
          event_.la_dau2_dxyos[event_.N1] = dxyos2;

          //qT AND alpha;

           double A2 = Angle(event_.la_px[event_.N1],event_.la_py[event_.N1],event_.la_pz[event_.N1],event_.la_dau1_px[event_.N1],event_.la_dau1_py[event_.N1],event_.la_dau1_pz[event_.N1]);
           //stop here!!02/21/2014:03:56pm
           double P2 = event_.la_dau1_p[event_.N1];
           double qT2 = P2*sin(A2);
           double qL1 = event_.la_dau1_p[event_.N1]*cos(A2);
           double qL2 = event_.la_dau2_p[event_.N1]*cos(A2);

           double alpha = (qL1-qL2)/(qL1+qL2);

           event_.la_qT[event_.N1] = qT2;
           event_.la_alpha[event_.N1] = alpha;

          //dEdx
/*
           if(dEdxHandle->size()){
                const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle.product();
                event_.la_dau1_dedx[event_.N1] = dEdxTrack[dau1].dEdx();
                event_.la_dau2_dedx[event_.N1] = dEdxTrack[dau2].dEdx();
            }
*/
          event_.N1++;

        }

}
void
V0AnalyzerPlugin::fillVertices(const edm::Event& iEvent){

  // Vertex 0 : event_vz[0] MC information from TrackingVertexCollection
  // Vertex 1 - n : Reconstructed Vertex from various of algorithms
  event_.vx[0]=0;
  event_.vy[0]=0;
  event_.vz[0]=0;
  event_.vxError[0]=0;
  event_.vyError[0]=0;
  event_.vzError[0]=0;
  event_.nDaughter[0]=0;

  if (doSimVertex_){
    Handle<TrackingVertexCollection> vertices;
    iEvent.getByLabel(simVertexSrc_, vertices);
    int greatestvtx = -1;
    for (unsigned int i = 0 ; i< vertices->size(); ++i){
      unsigned int daughter = (*vertices)[i].nDaughterTracks();
      if (greatestvtx==-1) {
         greatestvtx = i;
      } else {
         if( daughter > (*vertices)[greatestvtx].nDaughterTracks() && fabs((*vertices)[i].position().z()) < 30000) 
          greatestvtx = i;
      }
    }

    if(greatestvtx != -1){
      event_.vz[event_.nv] = (*vertices)[greatestvtx].position().z();
      }
    else{
      event_.vz[event_.nv] =  -99; 
    } 
    

  event_.nv++;

  }

  // Fill reconstructed vertices.   
  for(unsigned int iv = 0; iv < vertexSrc_.size(); ++iv){
    const reco::VertexCollection * recoVertices;
    edm::Handle<reco::VertexCollection> vertexCollection;
    //cout <<vertexSrc_[iv]<<endl;
    iEvent.getByLabel(vertexSrc_,vertexCollection);
    recoVertices = vertexCollection.product();
    unsigned int daughter = 0;
    int nVertex = 0;
    unsigned int greatestvtx = 0;

    nVertex = recoVertices->size();
    event_.nVtx = nVertex;
    for (unsigned int i = 0 ; i< recoVertices->size(); ++i){
      daughter = (*recoVertices)[i].tracksSize();
      if( daughter > (*recoVertices)[greatestvtx].tracksSize()) greatestvtx = i;

      event_.xVtx[i] = (*recoVertices)[i].position().x();
      event_.yVtx[i] = (*recoVertices)[i].position().y();
      event_.zVtx[i] = (*recoVertices)[i].position().z();
      event_.xVtxErr[i] = (*recoVertices)[i].xError();
      event_.yVtxErr[i] = (*recoVertices)[i].yError();
      event_.zVtxErr[i] = (*recoVertices)[i].zError();
      event_.nTrkVtx[i] = (*recoVertices)[i].tracksSize();
      event_.normChi2Vtx[i] = (*recoVertices)[i].normalizedChi2();

    
      float vtxSumPt=0.;
      for (reco::Vertex::trackRef_iterator it = (*recoVertices)[i].tracks_begin(); it != (*recoVertices)[i].tracks_end(); it++) {
  vtxSumPt += (**it).pt();

  Handle<vector<Track> > etracks;
  iEvent.getByLabel(trackSrc_, etracks);

  for(unsigned itrack=0; itrack<etracks->size(); ++itrack){
    reco::TrackRef trackRef=reco::TrackRef(etracks,itrack);
    //cout<<" trackRef.key() "<<trackRef.key()<< " it->key() "<<it->key()<<endl;
    if(trackRef.key()==it->key()){
      event_.trkVtxIndex[itrack] = i+1;  // note that index starts from 1 
      //cout<< " matching track "<<itrack<<endl;
    }
  }
      }

      event_.sumPtVtx[i] = vtxSumPt;
    
    }

    event_.maxVtx = greatestvtx;

    //loop over vertices again to get the significance wrt the leading vertex -Matt
    for (unsigned int i = 0 ; i< recoVertices->size(); ++i){
      if(i==greatestvtx) continue;
      GlobalVector direction = GlobalVector(event_.xVtx[i]-event_.xVtx[greatestvtx],event_.xVtx[i]-event_.xVtx[greatestvtx],event_.xVtx[i]-event_.xVtx[greatestvtx]);
      Measurement1D vtxDist2D = reco::SecondaryVertex::computeDist2d((*recoVertices)[greatestvtx], (*recoVertices)[i], direction, true);
      Measurement1D vtxDist3D = reco::SecondaryVertex::computeDist3d((*recoVertices)[greatestvtx], (*recoVertices)[i], direction, true);
      event_.vtxDist2D[i]=vtxDist2D.value();
      event_.vtxDist2DErr[i]=vtxDist2D.error();
      event_.vtxDist2DSig[i]=vtxDist2D.significance();
      event_.vtxDist3D[i]=vtxDist3D.value();
      event_.vtxDist3DErr[i]=vtxDist3D.error();
      event_.vtxDist3DSig[i]=vtxDist3D.significance();
    }

    if(recoVertices->size()>0){
      event_.vx[event_.nv] = (*recoVertices)[greatestvtx].position().x();
      event_.vy[event_.nv] = (*recoVertices)[greatestvtx].position().y();
      event_.vz[event_.nv] = (*recoVertices)[greatestvtx].position().z();
      event_.vxError[event_.nv] = (*recoVertices)[greatestvtx].xError();
      event_.vyError[event_.nv] = (*recoVertices)[greatestvtx].yError();
      event_.vzError[event_.nv] = (*recoVertices)[greatestvtx].zError();
      event_.nDaughter[event_.nv] = (*recoVertices)[greatestvtx].tracksSize();
    }else{
      event_.vx[event_.nv] =  -99;
      event_.vy[event_.nv] =  -99;
      event_.vz[event_.nv] =  -99;
      event_.vxError[event_.nv] =  -99;
      event_.vyError[event_.nv] =  -99;
      event_.vzError[event_.nv] =  -99;
      event_.nDaughter[event_.nv] = -99;
      
    }
    event_.nv++;
    //cout <<event_.nv<<endl;
  }


}

// ------------ method called once each job just before starting event loop  ------------
void 
V0AnalyzerPlugin::beginJob()
{
  edm::Service<TFileService> fs;
    
  TH1D::SetDefaultSumw2();

  v0Tree_ = fs->make<TTree>("v0Tree","v1");
  v0Tree_VTX_ = fs->make<TTree>("v0Tree_VTX","v4");
  v0_Kshort_ = fs->make<TTree>("v0_Kshort","v2");
  v0_Lambda_ = fs->make<TTree>("v0_Lambda","v3");
  PF_Jet_ = fs->make<TTree>("PFJet","v4");
  GeneralTracks_ = fs->make<TTree>("GeneralTracks","v5");
  GenParticle_ = fs->make<TTree>("GenParticle","v6");
  GenJet_ = fs->make<TTree>("GenJet","v7");
  Calo_Jet_ = fs->make<TTree>("Calo_Jet","v8");

   //v0Tree//----------------------------------------------
  v0Tree_->Branch("nEv",&event_.nEv,"nEv/I");
  v0Tree_->Branch("nLumi",&event_.nLumi,"nLumi/I");
  v0Tree_->Branch("nBX",&event_.nBX,"nBX/I");
  v0Tree_->Branch("nRun",&event_.nRun,"nRun/I");
  v0Tree_->Branch("N",&event_.N,"N/I");
  v0Tree_->Branch("N1",&event_.N1,"N1/I");
  v0Tree_->Branch("nJets", &event_.nJets, "nJets/I");
  v0Tree_->Branch("nTracks",&event_.nTracks,"nTracks/I");
  v0Tree_->Branch("EtaPtCutnTracks",&event_.EtaPtCutnTracks,"EtaPtCutnTracks/I");
  v0Tree_->Branch("nGen",&event_.nGen,"nGen/I");
  v0Tree_->Branch("nCaloJet",&event_.nCaloJet,"nCaloJet/I");
  v0Tree_->Branch("nGenJet", &event_.nGenJet,"nGenJet/I");


  //OfflinePrimaryVertex
  
  //v0Tree_->Branch("bestvx",&event_.bestvx,"bestvx/F");
  //v0Tree_->Branch("bestvy",&event_.bestvy,"bestvy/F");
  //v0Tree_->Branch("bestvz",&event_.bestvz,"bestvz/F");

  //GenJet:

  if(doGenJet_){

    GenJet_->Branch("nGenJet", &event_.nGenJet,"nGenJet/I");
    GenJet_->Branch("genpt", &event_.genpt,"genpt[nGenJet]/F");
    GenJet_->Branch("geneta", &event_.geneta,"geneta[nGenJet]/F");
    GenJet_->Branch("genphi", &event_.genphi,"genphi[nGenJet]/F");
    GenJet_->Branch("geny", &event_.geny,"geny[nGenJet]/F");
  
  }

  //CaloJet:

  if(doCaloJet_){

    Calo_Jet_->Branch("nCaloJet",&event_.nCaloJet,"nCaloJet/I");
    Calo_Jet_->Branch("calojet_pt",&event_.calojet_pt,"calojet_pt[nCaloJet]/F");
    Calo_Jet_->Branch("calojet_eta",&event_.calojet_eta,"calojet_eta[nCaloJet]/F");
    Calo_Jet_->Branch("calojet_phi",&event_.calojet_phi,"calojet_phi[nCaloJet]/F");
    Calo_Jet_->Branch("calojet_px",&event_.calojet_px,"calojet_px[nCaloJet]/F");
    Calo_Jet_->Branch("calojet_py",&event_.calojet_py,"calojet_py[nCaloJet]/F");
    Calo_Jet_->Branch("calojet_pz",&event_.calojet_pz,"calojet_pz[nCaloJet]/F");

  }

  //GenParticle:

  if(doGenParticle_){

    GenParticle_->Branch("nGen",&event_.nGen,"nGen/I");
    GenParticle_->Branch("genP_pt",&event_.genP_pt,"genP_pt[nGen]/F");
    GenParticle_->Branch("genP_eta",&event_.genP_eta,"genP_eta[nGen]/F");
    GenParticle_->Branch("genP_phi",&event_.genP_phi,"genP_phi[nGen]/F");
    GenParticle_->Branch("genP_mid",&event_.genP_mid,"genP_mid[nGen]/I");

    GenParticle_->Branch("pdg",&event_.pdg,"pdg[nGen]/I");

  }

  //GeneralTracks:

  if(doGeneralTracks_){
    
    GeneralTracks_->Branch("nTracks",&event_.nTracks,"nTracks/I");
    GeneralTracks_->Branch("track_eta",&event_.track_eta, "track_eta[nTracks]/F");
    GeneralTracks_->Branch("track_phi",&event_.track_phi, "track_phi[nTracks]/F");
    GeneralTracks_->Branch("track_pt",&event_.track_pt, "track_pt[nTracks]/F");
    GeneralTracks_->Branch("track_p",&event_.track_p, "track_p[nTracks]/F");
    GeneralTracks_->Branch("track_charge",&event_.track_charge, "track_charge/I");

  }

  //PF_Jet

  if(doPFJet_){
    
    PF_Jet_->Branch("nJets", &event_.nJets, "nJets/I");
  //  PF_Jet_->Branch("jet_chargedEmEnergy", &event_.jet_chargedEmEnergy, "jet_chargedEmEnergy[nJets]/F");
  // PF_Jet_->Branch("jet_maxDistance", &event_.jet_maxDistance, "jet_maxDistance[nJets]/F");
  //  PF_Jet_->Branch("jet_numberOfDaughters", &event_.jet_numberOfDaughters, "jet_numberOfDaughters/I");
    PF_Jet_->Branch("jet_vx", &event_.jet_vx, "jet_vx[nJets]/F");
    PF_Jet_->Branch("jet_vy", &event_.jet_vy, "jet_vy[nJets]/F");
    PF_Jet_->Branch("jet_vz", &event_.jet_vz, "jet_vz[nJets]/F");
    PF_Jet_->Branch("jet_pt", &event_.jet_pt, "jet_pt[nJets]/F");
    PF_Jet_->Branch("jet_px", &event_.jet_px, "jet_px[nJets]/F");
    PF_Jet_->Branch("jet_py", &event_.jet_py, "jet_py[nJets]/F");
    PF_Jet_->Branch("jet_pz", &event_.jet_pz, "jet_pz[nJets]/F");

    PF_Jet_->Branch("jet_p", &event_.jet_p, "jet_p[nJets]/F");
    PF_Jet_->Branch("jet_eta", &event_.jet_eta, "jet_eta[nJets]/F");
    PF_Jet_->Branch("jet_phi", &event_.jet_phi, "jet_phi[nJets]/F");


    }
  //v0_Kshort//----------------------------------------------
    //V0 candidate vertex;

  if(doV0_Kshort_){

    v0_Kshort_->Branch("N",&event_.N,"N/I");
  
    v0_Kshort_->Branch("ks_eta",&event_.ks_eta,"ks_eta[N]/F");
    v0_Kshort_->Branch("ks_phi",&event_.ks_phi,"ks_phi[N]/F");
    v0_Kshort_->Branch("ks_pt",&event_.ks_pt,"ks_pt[N]/F");
    v0_Kshort_->Branch("ks_mass",&event_.ks_mass,"ks_mass[N]/F");
    v0_Kshort_->Branch("ks_px",&event_.ks_px,"ks_px[N]/F");
    v0_Kshort_->Branch("ks_py",&event_.ks_py,"ks_py[N]/F");
    v0_Kshort_->Branch("ks_pz",&event_.ks_pz,"ks_pz[N]/F");
    v0_Kshort_->Branch("ks_charge",&event_.ks_charge,"ks_charge[N]/F");
    v0_Kshort_->Branch("ks_vtxChi2",&event_.ks_vtxChi2,"ks_vtxChi2[N]/F");
    v0_Kshort_->Branch("ks_agl",&event_.ks_agl,"ks_agl[N]/F");
    v0_Kshort_->Branch("ks_dl",&event_.ks_dl,"ks_dl[N]/F");
    v0_Kshort_->Branch("ks_dlerror",&event_.ks_dlerror,"ks_dlerror[N]/F");
    v0_Kshort_->Branch("ks_dlos",&event_.ks_dlos,"ks_dlos[N]/F");
 //   v0_Kshort_->Branch("ks_qT",&event_.ks_qT,"ks_qT[N]/F");
//    v0_Kshort_->Branch("ks_alpha",&event_.ks_alpha,"ks_alpha[N]/F");

    v0_Kshort_->Branch("ks_vx",&event_.ks_vx,"ks_vx[N]/F");
    v0_Kshort_->Branch("ks_vy",&event_.ks_vy,"ks_vy[N]/F");
    v0_Kshort_->Branch("ks_vz",&event_.ks_vz,"ks_vz[N]/F");
  //  v0_Kshort_->Branch("ks_vxError",&event_.ks_vxError,"ks_vxError[N]/F");
  //  v0_Kshort_->Branch("ks_vyError",&event_.ks_vyError,"ks_vyError[N]/F");
  //  v0_Kshort_->Branch("ks_vzError",&event_.ks_vzError,"ks_vzError[N]/F");

    //V0 candidate daughter

    v0_Kshort_->Branch("ks_dau1_eta",&event_.ks_dau1_eta,"ks_dau1_eta[N]/F");
    v0_Kshort_->Branch("ks_dau1_phi",&event_.ks_dau1_phi,"ks_dau1_phi[N]/F");
    v0_Kshort_->Branch("ks_dau1_pt",&event_.ks_dau1_pt,"ks_dau1_pt[N]/F");
    v0_Kshort_->Branch("ks_dau1_px",&event_.ks_dau1_px,"ks_dau1_px[N]/F");
    v0_Kshort_->Branch("ks_dau1_py",&event_.ks_dau1_py,"ks_dau1_py[N]/F");
    v0_Kshort_->Branch("ks_dau1_pz",&event_.ks_dau1_pz,"ks_dau1_pz[N]/F");
    v0_Kshort_->Branch("ks_dau1_charge",&event_.ks_dau1_charge,"ks_dau1_charge[N]/F");
    v0_Kshort_->Branch("ks_dau1_nhit",&event_.ks_dau1_nhit,"ks_dau1_nhit[N]/F");

    v0_Kshort_->Branch("ks_dau1_dzbest",&event_.ks_dau1_dzbest,"ks_dau1_dzbest[N]/F");
    v0_Kshort_->Branch("ks_dau1_dxybest",&event_.ks_dau1_dxybest,"ks_dau1_dxybest[N]/F");
    v0_Kshort_->Branch("ks_dau1_dzerror",&event_.ks_dau1_dzerror,"ks_dau1_dzerror[N]/F");
    v0_Kshort_->Branch("ks_dau1_dxyerror",&event_.ks_dau1_dxyerror,"ks_dau1_dxyerror[N]/F");
    v0_Kshort_->Branch("ks_dau1_dzos",&event_.ks_dau1_dzos,"ks_dau1_dzos[N]/F");
    v0_Kshort_->Branch("ks_dau1_dxyos",&event_.ks_dau1_dxyos,"ks_dau1_dxyos[N]/F");
    v0_Kshort_->Branch("ks_dau1_dedx",&event_.ks_dau1_dedx,"ks_dau1_dedx[N]/F");


    v0_Kshort_->Branch("ks_dau2_eta",&event_.ks_dau2_eta,"ks_dau2_eta[N]/F");
    v0_Kshort_->Branch("ks_dau2_phi",&event_.ks_dau2_phi,"ks_dau2_phi[N]/F");
    v0_Kshort_->Branch("ks_dau2_pt",&event_.ks_dau2_pt,"ks_dau2_pt[N]/F");
    v0_Kshort_->Branch("ks_dau2_px",&event_.ks_dau2_px,"ks_dau2_px[N]/F");
    v0_Kshort_->Branch("ks_dau2_py",&event_.ks_dau2_py,"ks_dau2_py[N]/F");
    v0_Kshort_->Branch("ks_dau2_pz",&event_.ks_dau2_pz,"ks_dau2_pz[N]/F");
    v0_Kshort_->Branch("ks_dau2_charge",&event_.ks_dau2_charge,"ks_dau2_charge[N]/F");
    v0_Kshort_->Branch("ks_dau2_nhit",&event_.ks_dau2_nhit,"ks_dau2_nhit[N]/F");

    v0_Kshort_->Branch("ks_dau2_dzbest",&event_.ks_dau2_dzbest,"ks_dau2_dzbest[N]/F");
    v0_Kshort_->Branch("ks_dau2_dxybest",&event_.ks_dau2_dxybest,"ks_dau2_dxybest[N]/F");
    v0_Kshort_->Branch("ks_dau2_dzerror",&event_.ks_dau2_dzerror,"ks_dau2_dzerror[N]/F");
    v0_Kshort_->Branch("ks_dau2_dxyerror",&event_.ks_dau2_dxyerror,"ks_dau2_dxyerror[N]/F");
    v0_Kshort_->Branch("ks_dau2_dzos",&event_.ks_dau2_dzos,"ks_dau2_dzos[N]/F");
    v0_Kshort_->Branch("ks_dau2_dxyos",&event_.ks_dau2_dxyos,"ks_dau2_dxyos[N]/F");
    v0_Kshort_->Branch("ks_dau2_dedx",&event_.ks_dau2_dedx,"ks_dau2_dedx[N]/F");

  }

  //v0_Lambda//------------------------------------------------------------

  if(doV0_Lambda_){

    v0_Lambda_->Branch("N1",&event_.N1,"N1/I");

    v0_Lambda_->Branch("la_eta",&event_.la_eta,"la_eta[N1]/F");
    v0_Lambda_->Branch("la_phi",&event_.la_phi,"la_phi[N1]/F");
    v0_Lambda_->Branch("la_pt",&event_.la_pt,"la_pt[N1]/F");
    v0_Lambda_->Branch("la_mass",&event_.la_mass,"la_mass[N1]/F");
    v0_Lambda_->Branch("la_px",&event_.la_px,"la_px[N1]/F");
    v0_Lambda_->Branch("la_py",&event_.la_py,"la_py[N1]/F");
    v0_Lambda_->Branch("la_pz",&event_.la_pz,"la_pz[N1]/F");
    v0_Lambda_->Branch("la_charge",&event_.la_charge,"la_charge[N1]/F");
    v0_Lambda_->Branch("la_vtxChi2",&event_.la_vtxChi2,"la_vtxChi2[N1]/F");
    v0_Lambda_->Branch("la_agl",&event_.la_agl,"la_agl[N1]/F");
    v0_Lambda_->Branch("la_dl",&event_.la_dl,"la_dl[N1]/F");
    v0_Lambda_->Branch("la_dlerror",&event_.la_dlerror,"la_dlerror[N1]/F");
    v0_Lambda_->Branch("la_dlos",&event_.la_dlos,"la_dlos[N1]/F");
  //  v0_Lambda_->Branch("la_qT",&event_.la_qT,"la_qT[N1]/F");
  //  v0_Lambda_->Branch("la_alpha",&event_.la_alpha,"la_alpha[N1]/F");

    v0_Lambda_->Branch("la_vx",&event_.la_vx,"la_vx[N1]/F");
    v0_Lambda_->Branch("la_vy",&event_.la_vy,"la_vy[N1]/F");
    v0_Lambda_->Branch("la_vz",&event_.la_vz,"la_vz[N1]/F");
  //  v0_Lambda_->Branch("la_vxError",&event_.la_vxError,"la_vxError[N1]/F");
  //  v0_Lambda_->Branch("la_vyError",&event_.la_vyError,"la_vyError[N1]/F");
   // v0_Lambda_->Branch("la_vzError",&event_.la_vzError,"la_vzError[N1]/F");

    //V0 candidate daughter

    v0_Lambda_->Branch("la_dau1_eta",&event_.la_dau1_eta,"la_dau1_eta[N1]/F");
    v0_Lambda_->Branch("la_dau1_phi",&event_.la_dau1_phi,"la_dau1_phi[N1]/F");
    v0_Lambda_->Branch("la_dau1_pt",&event_.la_dau1_pt,"la_dau1_pt[N1]/F");
    v0_Lambda_->Branch("la_dau1_px",&event_.la_dau1_px,"la_dau1_px[N1]/F");
    v0_Lambda_->Branch("la_dau1_py",&event_.la_dau1_py,"la_dau1_py[N1]/F");
    v0_Lambda_->Branch("la_dau1_pz",&event_.la_dau1_pz,"la_dau1_pz[N1]/F");
    v0_Lambda_->Branch("la_dau1_charge",&event_.la_dau1_charge,"la_dau1_charge[N1]/F");
    v0_Lambda_->Branch("la_dau1_nhit",&event_.la_dau1_nhit,"la_dau1_nhit[N1]/F");

    v0_Lambda_->Branch("la_dau1_dzbest",&event_.la_dau1_dzbest,"la_dau1_dzbest[N1]/F");
    v0_Lambda_->Branch("la_dau1_dxybest",&event_.la_dau1_dxybest,"la_dau1_dxybest[N1]/F");
    v0_Lambda_->Branch("la_dau1_dzerror",&event_.la_dau1_dzerror,"la_dau1_dzerror[N1]/F");
    v0_Lambda_->Branch("la_dau1_dxyerror",&event_.la_dau1_dxyerror,"la_dau1_dxyerror[N1]/F");
    v0_Lambda_->Branch("la_dau1_dzos",&event_.la_dau1_dzos,"la_dau1_dzos[N1]/F");
    v0_Lambda_->Branch("la_dau1_dxyos",&event_.la_dau1_dxyos,"la_dau1_dxyos[N1]/F");
    v0_Lambda_->Branch("la_dau1_dedx",&event_.la_dau1_dedx,"la_dau1_dedx[N1]/F");


    v0_Lambda_->Branch("la_dau2_eta",&event_.la_dau2_eta,"la_dau2_eta[N1]/F");
    v0_Lambda_->Branch("la_dau2_phi",&event_.la_dau2_phi,"la_dau2_phi[N1]/F");
    v0_Lambda_->Branch("la_dau2_pt",&event_.la_dau2_pt,"la_dau2_pt[N1]/F");
    v0_Lambda_->Branch("la_dau2_px",&event_.la_dau2_px,"la_dau2_px[N1]/F");
    v0_Lambda_->Branch("la_dau2_py",&event_.la_dau2_py,"la_dau2_py[N1]/F");
    v0_Lambda_->Branch("la_dau2_pz",&event_.la_dau2_pz,"la_dau2_pz[N1]/F");
    v0_Lambda_->Branch("la_dau2_charge",&event_.la_dau2_charge,"la_dau2_charge[N1]/F");
    v0_Lambda_->Branch("la_dau2_nhit",&event_.la_dau2_nhit,"la_dau2_nhit[N1]/F");

    v0_Lambda_->Branch("la_dau2_dzbest",&event_.la_dau2_dzbest,"la_dau2_dzbest[N1]/F");
    v0_Lambda_->Branch("la_dau2_dxybest",&event_.la_dau2_dxybest,"la_dau2_dxybest[N1]/F");
    v0_Lambda_->Branch("la_dau2_dzerror",&event_.la_dau2_dzerror,"la_dau2_dzerror[N1]/F");
    v0_Lambda_->Branch("la_dau2_dxyerror",&event_.la_dau2_dxyerror,"la_dau2_dxyerror[N1]/F");
    v0_Lambda_->Branch("la_dau2_dzos",&event_.la_dau2_dzos,"la_dau2_dzos[N1]/F");
    v0_Lambda_->Branch("la_dau2_dxyos",&event_.la_dau2_dxyos,"la_dau2_dxyos[N1]/F");
    v0_Lambda_->Branch("la_dau2_dedx",&event_.la_dau2_dedx,"la_dau2_dedx[N1]/F");

  }
    
    //vertex----------------------------------------------------------------------------
    if(doVertex_){

    v0Tree_VTX_->Branch("nv",&event_.nv,"nv/I");
    v0Tree_VTX_->Branch("vx",event_.vx,"vx[nv]/F");
    v0Tree_VTX_->Branch("vy",event_.vy,"vy[nv]/F");
    v0Tree_VTX_->Branch("vz",event_.vz,"vz[nv]/F");
    v0Tree_VTX_->Branch("vxErr",event_.vxError,"vxErr[nv]/F"); 
    v0Tree_VTX_->Branch("vyErr",event_.vyError,"vyErr[nv]/F"); 
    v0Tree_VTX_->Branch("vzErr",event_.vzError,"vzErr[nv]/F");
    v0Tree_VTX_->Branch("nDaughter",event_.nDaughter,"nDaughter[nv]/I");

    v0Tree_VTX_->Branch("nVtx",&event_.nVtx,"nVtx/I");
    v0Tree_VTX_->Branch("maxVtx",&event_.maxVtx,"maxVtx/I");
    //  trackTree_->Branch("maxVtxHard",&pev_.maxVtxHard,"maxVtxHard/I");

    v0Tree_VTX_->Branch("nTrkVtx",event_.nTrkVtx,"nTrkVtx[nVtx]/I");
    v0Tree_VTX_->Branch("normChi2Vtx",event_.normChi2Vtx,"normChi2Vtx[nVtx]/F");
    v0Tree_VTX_->Branch("sumPtVtx",event_.sumPtVtx,"sumPtVtx[nVtx]/F");
    //  trackTree_->Branch("nTrkVtxHard",pev_.nTrkVtxHard,"nTrkVtxHard[nVtx]/I");

    v0Tree_VTX_->Branch("xVtx",event_.xVtx,"xVtx[nVtx]/F");
    v0Tree_VTX_->Branch("yVtx",event_.yVtx,"yVtx[nVtx]/F");
    v0Tree_VTX_->Branch("zVtx",event_.zVtx,"zVtx[nVtx]/F");
    v0Tree_VTX_->Branch("xVtxErr",event_.xVtxErr,"xVtxErr[nVtx]/F");
    v0Tree_VTX_->Branch("yVtxErr",event_.yVtxErr,"yVtxErr[nVtx]/F");
    v0Tree_VTX_->Branch("zVtxErr",event_.zVtxErr,"zVtxErr[nVtx]/F");

    v0Tree_VTX_->Branch("vtxDist2D",event_.vtxDist2D,"vtxDist2D[nVtx]/F");
    v0Tree_VTX_->Branch("vtxDist2DErr",event_.vtxDist2DErr,"vtxDist2DErr[nVtx]/F");
    v0Tree_VTX_->Branch("vtxDist2DSig",event_.vtxDist2DSig,"vtxDist2DSig[nVtx]/F");
    v0Tree_VTX_->Branch("vtxDist2D",event_.vtxDist3D,"vtxDist3D[nVtx]/F");
    v0Tree_VTX_->Branch("vtxDist3DErr",event_.vtxDist3DErr,"vtxDist3DErr[nVtx]/F");
    v0Tree_VTX_->Branch("vtxDist3DSig",event_.vtxDist3DSig,"vtxDist3DSig[nVtx]/F");
   
}
 
 //-------------------------------------------------------------------------------------
  


}

// ------------ method called once each job just after ending the event loop  ------------
void 
V0AnalyzerPlugin::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
V0AnalyzerPlugin::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
V0AnalyzerPlugin::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
V0AnalyzerPlugin::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
V0AnalyzerPlugin::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
V0AnalyzerPlugin::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(V0AnalyzerPlugin);
