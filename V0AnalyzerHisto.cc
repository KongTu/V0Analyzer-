//cription: [one line class summary]

// Implementation:
  //   [Notes on implementation]

//
// Original Author:  Zhoudunming Tu,
//         Created:  Mon Jun 13 20:56:30 CEST 2011
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
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
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

// Root include files
#include "TTree.h"
//
// Track Matching and fake rate calculations     
//#include "RiceHIG/V0Analysis/interface/V0Validator.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//
// class decleration
//

class V0AnalyzerHisto : public edm::EDAnalyzer {
public:
  explicit V0AnalyzerHisto(const edm::ParameterSet&);
  ~V0AnalyzerHisto();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;



  // ----------member data ---------------------------

  edm::InputTag trackSrc_;
  edm::InputTag simVertexSrc_;
  edm::InputTag generalV0_ks_;
  edm::InputTag generalV0_la_;
  std::string vertexSrc_;
  std::string jetSrc_;
  edm::InputTag genParticleSrc_;

  TH1D* InvMass_la[10];
  TH1D* InvMass_ks_underlying[10];
  TH1D* InvMass_la_underlying[10];
  TH1D* InvMass_ks[10];

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

V0AnalyzerHisto::V0AnalyzerHisto(const edm::ParameterSet& iConfig)
 
{
  trackSrc_ = iConfig.getParameter<edm::InputTag>("trackSrc");
  vertexSrc_ = iConfig.getParameter<std::string>("vertexSrc");
  simVertexSrc_ =  iConfig.getUntrackedParameter<edm::InputTag>("tpVtxSrc",edm::InputTag("mergedtruth","MergedTrackTruth"));
  generalV0_ks_ = iConfig.getParameter<edm::InputTag>("generalV0_ks");
  generalV0_la_ = iConfig.getParameter<edm::InputTag>("generalV0_la");
  jetSrc_ = iConfig.getParameter<std::string>("jetSrc");
  genParticleSrc_ = iConfig.getParameter<edm::InputTag>("genParticleSrc");
  
}


V0AnalyzerHisto::~V0AnalyzerHisto()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//==================
// member functions
//==================

double Mass_ks(double px_1,double py_1,double pz_1,double px_2,double py_2,double pz_2)
{
       
  double temp = 0.0;
        double E1 = sqrt((px_1*px_1+py_1*py_1+pz_1*pz_1)+(0.13957*0.13957));
        double E2 = sqrt((px_2*px_2+py_2*py_2+pz_2*pz_2)+(0.93827*0.93827));
        double E_tot = E1+E2;
  temp = (E_tot*E_tot) - ((px_1+px_2)*(px_1+px_2)+(py_1+py_2)*(py_1+py_2)+(pz_1+pz_2)*(pz_1+pz_2));
  return sqrt(temp);
}

double Mass_la(double px_1,double py_1,double pz_1,double px_2,double py_2,double pz_2)
{
       
  double temp = 0.0;
        double E1 = sqrt((px_1*px_1+py_1*py_1+pz_1*pz_1)+(0.13957*0.13957));
        double E2 = sqrt((px_2*px_2+py_2*py_2+pz_2*pz_2)+(0.13957*0.13957));
        double E_tot = E1+E2;
  temp = (E_tot*E_tot) - ((px_1+px_2)*(px_1+px_2)+(py_1+py_2)*(py_1+py_2)+(pz_1+pz_2)*(pz_1+pz_2));
  return sqrt(temp);
}

double Mass_e(double px_1,double py_1,double pz_1,double px_2,double py_2,double pz_2)
{
        double temp = 0.0;
        double E1 = sqrt((px_1*px_1+py_1*py_1+pz_1*pz_1)+(0.000511*0.000511));
        double E2 = sqrt((px_2*px_2+py_2*py_2+pz_2*pz_2)+(0.000511*0.000511));
        double E_tot = E1+E2;
        temp = (E_tot*E_tot) - ((px_1+px_2)*(px_1+px_2)+(py_1+py_2)*(py_1+py_2)+(pz_1+pz_2)*(pz_1+pz_2));
        return sqrt(temp);
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

// ------------ method called to for each event  ------------
void
V0AnalyzerHisto::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vertexSrc_,vertices);
  double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
  double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
  const reco::Vertex & vtx = (*vertices)[0];
  bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
  bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();

  Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(trackSrc_, tracks);

  int nTracks = 0;

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
      
  //ntrack selection:

        if ( fabs(trk.eta()) > 2.4 || trk.pt() < 0.4  ) continue;

        nTracks++;
        
  } 

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByLabel( jetSrc_ , jets);
   // jetSrc = "ak5PFJets" inputTag (string);

  edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_ks;
    iEvent.getByLabel(generalV0_ks_,v0candidates_ks);
    if(!v0candidates_ks.isValid()) return;
    
  edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_la;
    iEvent.getByLabel(generalV0_la_,v0candidates_la);
    if(!v0candidates_la.isValid()) return;
 
//first selection; vertices
    if(bestvz < -15.0 || bestvz>15.0) return;

//multiplicity bins:
    //if (nTracks < 220 || nTracks > 260) return;


  for(unsigned it=0; it<v0candidates_ks->size(); ++it){     
  
          const reco::VertexCompositeCandidate & trk = (*v0candidates_ks)[it];
          const reco:: Candidate * d1 = trk.daughter(0);
          const reco:: Candidate * d2 = trk.daughter(1); 

          auto dau1 = d1->get<reco::TrackRef>();
          auto dau2 = d2->get<reco::TrackRef>();
   
          double eta_dau1 = dau1->eta();
          double phi_dau1 = dau1->phi();
          double pt_dau1 = dau1->pt();
          double px_dau1 = dau1->px();
          double py_dau1 = dau1->py();
          double pz_dau1 = dau1->pz();
          double charge_dau1 = dau1->charge();  
          double p_dau1 = dau1->p();        

          double eta_dau2 = dau2->eta();
          double phi_dau2 = dau2->phi();
          double pt_dau2 = dau2->pt();
          double px_dau2 = dau2->px();
          double py_dau2 = dau2->py();
          double pz_dau2 = dau2->pz();    
          double charge_dau2 = dau2->charge();
          double p_dau2 = dau2->p();

          double ks_mass = trk.mass();
          double ks_pt = trk.pt();
          double ks_px = trk.px();
          double ks_py = trk.py();
          double ks_pz = trk.pz();
          double ks_eta = trk.eta();
          double ks_phi = trk.phi();

          //PAngle
          double secvz = trk.vz();
          double secvx = trk.vx();
          double secvy = trk.vy();
          TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
          TVector3 secvec(ks_px,ks_py,ks_pz);

          double agl = cos(secvec.Angle(ptosvec));

         //Decay length
          typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
          typedef ROOT::Math::SVector<double, 3> SVector3;

          SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
          SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);

          double dl = ROOT::Math::Mag(distanceVector);
          double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
          double dlos = dl/dlerror;
          
          //NumberofValidHits for two daughters"
          double dau1_Nhits = dau1->numberOfValidHits();
          double dau2_Nhits = dau2->numberOfValidHits();

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

//underlying K0s:

          if ( dau1_Nhits > 3 && dau2_Nhits > 3 && dlos > 5 && agl > 0.999 && TMath::Abs(dzos1) > 1 && TMath::Abs(dzos2) > 1 && TMath::Abs(dxyos1) > 1 && TMath::Abs(dxyos2) > 1)
          {

              if (ks_pt > 0.7 && ks_pt < 1.0)
                { 
                    double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 1.13568 && temp > 1.09568 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_ks_underlying[0]->Fill(ks_mass);
                }

              if (ks_pt > 1.0 && ks_pt < 1.4)
                {
                    double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 1.13568 && temp > 1.09568 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_ks_underlying[1]->Fill(ks_mass);
                }

              if (ks_pt > 1.4 && ks_pt < 1.8)
                {
                    double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 1.13568 && temp > 1.09568 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_ks_underlying[2]->Fill(ks_mass);
                }

              if (ks_pt > 1.8 && ks_pt < 2.2)
                {
                    double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                     if ( temp < 1.13568 && temp > 1.09568 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_ks_underlying[3]->Fill(ks_mass);
                }

              if (ks_pt > 2.2 && ks_pt < 2.8)
                {
                    double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 1.13568 && temp > 1.09568 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_ks_underlying[4]->Fill(ks_mass);
                }


              if (ks_pt > 2.8 && ks_pt < 3.6)
                {
                    double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 1.13568 && temp > 1.09568 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_ks_underlying[5]->Fill(ks_mass);
                }

              if (ks_pt > 3.6 && ks_pt < 4.6)
                {
                    double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 1.13568 && temp > 1.09568 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_ks_underlying[6]->Fill(ks_mass);
                }

              if (ks_pt > 4.6 && ks_pt < 6.0)
                {
                    double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 1.13568 && temp > 1.09568 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_ks_underlying[7]->Fill(ks_mass);
                }

              if (ks_pt > 6.0 && ks_pt < 9.0)
                {
                    double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 1.13568 && temp > 1.09568 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_ks_underlying[8]->Fill(ks_mass);
                }

              if (ks_pt > 9.0 && ks_pt < 12.0)
                {
                    double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 1.13568 && temp > 1.09568 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_ks_underlying[9]->Fill(ks_mass);
                }

          }


//in jet K0s:

          for( unsigned it = 0; it < jets->size(); ++it ){

              const reco::PFJet & pfj = (*jets)[it];

              double jet_pt = pfj.pt();
              double jet_eta = pfj.eta();
              double jet_phi = pfj.phi();
              double jet_p = pfj.p();

              if ( jet_pt < 20 || TMath::Abs(jet_eta) > 2.0 ) continue;

              double delta_ks_eta = (jet_eta) - (ks_eta);
              double delta_ks_phi = (jet_phi) - (ks_phi);
              double KSconeSize = 0.0;

              if ( delta_ks_phi > 3.14 ){

                 KSconeSize = sqrt((6.28 - delta_ks_phi)*(6.28 - delta_ks_phi)+(delta_ks_eta)*(delta_ks_eta));
              }
              else if ( delta_ks_phi < -3.14 ){

                 KSconeSize = sqrt((6.28 + delta_ks_phi)*(6.28 + delta_ks_phi)+(delta_ks_eta)*(delta_ks_eta));

              }
              else{

                 KSconeSize = sqrt((delta_ks_phi)*(delta_ks_phi)+(delta_ks_eta)*(delta_ks_eta));

              }

              if ( KSconeSize < 0.3 ){

          
                if ( dau1_Nhits > 3 && dau2_Nhits > 3 && dlos > 5 && agl > 0.999 && TMath::Abs(dzos1) > 1 && TMath::Abs(dzos2) > 1 && TMath::Abs(dxyos1) > 1 && TMath::Abs(dxyos2) > 1)
                {
      //ss
                      if (ks_pt > 0.7 && ks_pt < 1.0)
                        {
                            double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                            double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              if ( temp < 1.13568 && temp > 1.09568 ) continue;
                                if ( temp_e < 0.015511) continue;

                                  InvMass_ks[0]->Fill(ks_mass);
                        }

                      if (ks_pt > 1.0 && ks_pt < 1.4)
                        {
                            double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                            double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              if ( temp < 1.13568 && temp > 1.09568 ) continue;
                                if ( temp_e < 0.015511) continue;

                                  InvMass_ks[1]->Fill(ks_mass);
                        }

                      if (ks_pt > 1.4 && ks_pt < 1.8)
                        {
                            double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                            double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              if ( temp < 1.13568 && temp > 1.09568 ) continue;
                                if ( temp_e < 0.015511) continue;

                                  InvMass_ks[2]->Fill(ks_mass);
                        }

                      if (ks_pt > 1.8 && ks_pt < 2.2)
                        {
                          double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                          double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                           if ( temp < 1.13568 && temp > 1.09568 ) continue;
                              if ( temp_e < 0.015511) continue;

                                InvMass_ks[3]->Fill(ks_mass);
                        }

                      if (ks_pt > 2.2 && ks_pt < 2.8)
                        {
                          double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                          double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                            if ( temp < 1.13568 && temp > 1.09568 ) continue;
                              if ( temp_e < 0.015511) continue;

                                InvMass_ks[4]->Fill(ks_mass);
                        }

                      if (ks_pt > 2.8 && ks_pt < 3.6)
                        {
                            double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                            double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              if ( temp < 1.13568 && temp > 1.09568 ) continue;
                                if ( temp_e < 0.015511) continue;

                                  InvMass_ks[5]->Fill(ks_mass);
                        }

                      if (ks_pt > 3.6 && ks_pt < 4.6)
                        {
                            double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                            double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              if ( temp < 1.13568 && temp > 1.09568 ) continue;
                                if ( temp_e < 0.015511) continue;

                                  InvMass_ks[6]->Fill(ks_mass);
                        }

                      if (ks_pt > 4.6 && ks_pt < 6.0)
                        {
                            double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                            double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              if ( temp < 1.13568 && temp > 1.09568 ) continue;
                                if ( temp_e < 0.015511) continue;

                                  InvMass_ks[7]->Fill(ks_mass);
                        }

                      if (ks_pt > 6.0 && ks_pt < 9.0)
                        {
                            double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                            double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              if ( temp < 1.13568 && temp > 1.09568 ) continue;
                                if ( temp_e < 0.015511) continue;

                                  InvMass_ks[8]->Fill(ks_mass);
                        }

                      if (ks_pt > 9.0 && ks_pt < 12.0)
                        {
                            double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                            double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              if ( temp < 1.13568 && temp > 1.09568 ) continue;
                                if ( temp_e < 0.015511) continue;

                                  InvMass_ks[9]->Fill(ks_mass);
                        }

                  
                }
              }  
          }

    }

//la V0 candidate:

    for(unsigned it=0; it<v0candidates_la->size(); ++it){

          const reco::VertexCompositeCandidate & trk = (*v0candidates_la)[it];
          const reco:: Candidate * d1 = trk.daughter(0);
          const reco:: Candidate * d2 = trk.daughter(1);

          auto dau1 = d1->get<reco::TrackRef>();
          auto dau2 = d2->get<reco::TrackRef>();

          double eta_dau1 = dau1->eta();
          double phi_dau1 = dau1->phi();
          double pt_dau1 = dau1->pt();
          double px_dau1 = dau1->px();
          double py_dau1 = dau1->py();
          double pz_dau1 = dau1->pz();
          double charge_dau1 = dau1->charge();
          double p_dau1 = dau1->p();

          double eta_dau2 = dau2->eta();
          double phi_dau2 = dau2->phi();
          double pt_dau2 = dau2->pt();
          double px_dau2 = dau2->px();
          double py_dau2 = dau2->py();
          double pz_dau2 = dau2->pz();     
          double charge_dau2 = dau2->charge();
          double p_dau2 = dau2->p();  
  
          double la_mass = trk.mass();
          double la_pt = trk.pt();
          double la_px = trk.px();
          double la_py = trk.py();
          double la_pz = trk.pz();
          double la_eta = trk.eta();
          double la_phi = trk.phi();

          //PAngle
          double secvz = trk.vz();
          double secvx = trk.vx();
          double secvy = trk.vy();
          TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
          TVector3 secvec(la_px,la_py,la_pz);

          double agl = cos(secvec.Angle(ptosvec));

         //Decay length
          typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
          typedef ROOT::Math::SVector<double, 3> SVector3;

          SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
          SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);

          double dl = ROOT::Math::Mag(distanceVector);
          double dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
          double dlos = dl/dlerror;
          
          //NumberofValidHits for two daughters"
          double dau1_Nhits = dau1->numberOfValidHits();
          double dau2_Nhits = dau2->numberOfValidHits();

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

//underlying Lambdas:

          if ( dau1_Nhits > 3 && dau2_Nhits > 3 && dlos > 5 && agl > 0.999 && TMath::Abs(dzos1) > 1 && TMath::Abs(dzos2) > 1 && TMath::Abs(dxyos1) > 1 && TMath::Abs(dxyos2) > 1)
          {

              if (la_pt > 0.7 && la_pt < 1.0)
                { 
                    double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 0.5076 && temp > 0.4876 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_la_underlying[0]->Fill(la_mass);
                }

              if (la_pt > 1.0 && la_pt < 1.4)
                {
                    double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 0.5076 && temp > 0.4876 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_la_underlying[1]->Fill(la_mass);
                }

              if (la_pt > 1.4 && la_pt < 1.8)
                {
                    double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 0.5076 && temp > 0.4876 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_la_underlying[2]->Fill(la_mass);
                }

              if (la_pt > 1.8 && la_pt < 2.2)
                {
                    double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                     if ( temp < 0.5076 && temp > 0.4876 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_la_underlying[3]->Fill(la_mass);
                }

              if (la_pt > 2.2 && la_pt < 2.8)
                {
                    double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 0.5076 && temp > 0.4876 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_la_underlying[4]->Fill(la_mass);
                }


              if (la_pt > 2.8 && la_pt < 3.6)
                {
                    double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 0.5076 && temp > 0.4876 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_la_underlying[5]->Fill(la_mass);
                }

              if (la_pt > 3.6 && la_pt < 4.6)
                {
                    double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 0.5076 && temp > 0.4876 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_la_underlying[6]->Fill(la_mass);
                }

              if (la_pt > 4.6 && la_pt < 6.0)
                {
                    double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 0.5076 && temp > 0.4876 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_la_underlying[7]->Fill(la_mass);
                }

              if (la_pt > 6.0 && la_pt < 9.0)
                {
                    double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 0.5076 && temp > 0.4876 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_la_underlying[8]->Fill(la_mass);
                }

              if (la_pt > 9.0 && la_pt < 12.0)
                {
                    double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                    double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                      if ( temp < 0.5076 && temp > 0.4876 ) continue;
                        if ( temp_e < 0.015511) continue;

                          InvMass_la_underlying[9]->Fill(la_mass);
                }

          }

//in jet K0s:

          for( unsigned it = 0; it < jets->size(); ++it ){

              const reco::PFJet & pfj = (*jets)[it];

              double jet_pt = pfj.pt();
              double jet_eta = pfj.eta();
              double jet_phi = pfj.phi();
              double jet_p = pfj.p();

              if ( jet_pt < 20 || TMath::Abs(jet_eta) > 2.0 ) continue;

              double delta_la_eta = (jet_eta) - (la_eta);
              double delta_la_phi = (jet_phi) - (la_phi);
              double LAconeSize = 0.0;

              if ( delta_la_phi > 3.14 ){

                 LAconeSize = sqrt((6.28 - delta_la_phi)*(6.28 - delta_la_phi)+(delta_la_eta)*(delta_la_eta));
              }
              else if ( delta_la_phi < -3.14 ){

                 LAconeSize = sqrt((6.28 + delta_la_phi)*(6.28 + delta_la_phi)+(delta_la_eta)*(delta_la_eta));

              }
              else{

                 LAconeSize = sqrt((delta_la_phi)*(delta_la_phi)+(delta_la_eta)*(delta_la_eta));

              }

              if( LAconeSize < 0.3 ){

                 if (dau1_Nhits > 3 && dau2_Nhits > 3 && dlos > 5 && agl > 0.999 && TMath::Abs(dzos1) > 1 && TMath::Abs(dzos2) > 1 && TMath::Abs(dxyos1) > 1 && TMath::Abs(dxyos2) > 1 ){

                      if (la_pt > 0.7 && la_pt < 1.0)
                        {
                          
                            double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                            double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              if ( temp < 0.5076 && temp > 0.4876 ) continue;
                                if ( temp_e < 0.015511) continue;

                                  InvMass_la[0]->Fill(la_mass);
                        }

                      if (la_pt > 1.0 && la_pt < 1.4)
                        {
                          
                            double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                            double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              if ( temp < 0.5076 && temp > 0.4876 ) continue;
                                if ( temp_e < 0.015511) continue;

                                  InvMass_la[1]->Fill(la_mass);
                        }

                      if (la_pt > 1.4 && la_pt < 1.8)
                        {
                          
                            double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                            double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              if ( temp < 0.5076 && temp > 0.4876 ) continue;
                                if ( temp_e < 0.015511) continue;

                                  InvMass_la[2]->Fill(la_mass);
                        }

                        if (la_pt > 1.8 && la_pt < 2.2)
                        {
                          
                            double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                            double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                             if ( temp < 0.5076 && temp > 0.4876 ) continue;
                                if ( temp_e < 0.015511) continue;

                                  InvMass_la[3]->Fill(la_mass);
                        }

                        if (la_pt > 2.2 && la_pt < 2.8)
                          {
                           
                            double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                            double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              if ( temp < 0.5076 && temp > 0.4876 ) continue;
                                if ( temp_e < 0.015511) continue;

                                  InvMass_la[4]->Fill(la_mass);
                          }


                        if (la_pt > 2.8 && la_pt < 3.6)
                          {
                            
                              double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                                if ( temp < 0.5076 && temp > 0.4876 ) continue;
                                  if ( temp_e < 0.015511) continue;

                                    InvMass_la[5]->Fill(la_mass);
                          }

                        if (la_pt > 3.6 && la_pt < 4.6)
                          {
                            
                              double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                                if ( temp < 0.5076 && temp > 0.4876 ) continue;
                                  if ( temp_e < 0.015511) continue;

                                    InvMass_la[6]->Fill(la_mass);
                          }

                        if (la_pt > 4.6 && la_pt < 6.0)
                          {
                            
                              double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                                if ( temp < 0.5076 && temp > 0.4876 ) continue;
                                  if ( temp_e < 0.015511) continue;

                                    InvMass_la[7]->Fill(la_mass);
                          }

                        if (la_pt > 6.0 && la_pt < 9.0)
                          {
                            
                              double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                                if ( temp < 0.5076 && temp > 0.4876 ) continue;
                                  if ( temp_e < 0.015511) continue;

                                    InvMass_la[8]->Fill(la_mass);
                          }

                        if (la_pt > 9.0 && la_pt < 12.0)
                          {
                            
                              double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                              double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
                                if ( temp < 0.5076 && temp > 0.4876 ) continue;
                                  if ( temp_e < 0.015511) continue;

                                    InvMass_la[9]->Fill(la_mass);
                          }
                    }
              }
          }
    }
}
// ------------ method called once each job just before starting event loop  ------------
void 
V0AnalyzerHisto::beginJob()
{
  edm::Service<TFileService> fs;
    
  TH1D::SetDefaultSumw2();

    for(int i = 0; i<10; i++)
     {
       InvMass_ks_underlying[i] = fs->make<TH1D>(Form("InvMass_ks_underlying%d",i),";ks_mass(GeV/c^{2});#Events",360,0.44,0.56);
       InvMass_ks[i] = fs->make<TH1D>(Form("InvMass_ks%d",i),";ks_mass(GeV);#Events",360,0.44,0.56);
       InvMass_la_underlying[i] = fs->make<TH1D>(Form("InvMass_la_underlying%d",i),";la_mass(GeV/c^{2});#Events",360,1.08,1.16);
       InvMass_la[i] = fs->make<TH1D>(Form("InvMass_la%d",i),";la_mass(GeV);#Events",360,1.08,1.16);
     }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
V0AnalyzerHisto::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(V0AnalyzerHisto);