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
#include <sstream>


#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
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
  edm::InputTag generalV0_xi_;
  edm::InputTag genParticleSrc_;
  std::string vertexSrc_;
  std::string jetSrc_;
  int multmin_;
  int multmax_;

  bool doGenParticle_;
  
  TH3D* InvMass_ks_underlying;
  TH3D* InvMass_la_underlying;

  TH3D* XiDaughter;

  TH3D* genKS_underlying;
  TH3D* genLA_underlying;

  TH1D* vertexDistZ;
  TH1D* xiMass;

  TH2D* ks_ptRapidity;
  TH2D* la_ptRapidity;

  TH1D* multiDist;
  TH1D* etaDist;


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
  generalV0_xi_ = iConfig.getParameter<edm::InputTag>("generalV0_xi");

  jetSrc_ = iConfig.getParameter<std::string>("jetSrc");
  genParticleSrc_ = iConfig.getParameter<edm::InputTag>("genParticleSrc");
  multmin_ = iConfig.getUntrackedParameter<int>("multmin", 120);
  multmax_ = iConfig.getUntrackedParameter<int>("multmax", 150); 

  doGenParticle_ = iConfig.getUntrackedParameter<bool>("doGenParticle",false);
  
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
        double E1 = sqrt((px_1*px_1+py_1*py_1+pz_1*pz_1)+(0.93827203*0.93827203));
        double E2 = sqrt((px_2*px_2+py_2*py_2+pz_2*pz_2)+(0.13957018*0.13957018));
        double E_tot = E1+E2;
  temp = (E_tot*E_tot) - ((px_1+px_2)*(px_1+px_2)+(py_1+py_2)*(py_1+py_2)+(pz_1+pz_2)*(pz_1+pz_2));
  return sqrt(temp);
}

double Mass_la(double px_1,double py_1,double pz_1,double px_2,double py_2,double pz_2)
{
       
  double temp = 0.0;
        double E1 = sqrt((px_1*px_1+py_1*py_1+pz_1*pz_1)+(0.13957018*0.13957018));
        double E2 = sqrt((px_2*px_2+py_2*py_2+pz_2*pz_2)+(0.13957018*0.13957018));
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
  using namespace std;

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByLabel(vertexSrc_,vertices);
  double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
  double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
  const reco::Vertex & vtx = (*vertices)[0];
  bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
  bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();


  vertexDistZ->Fill( vtx.z() );
  

  //first selection; vertices
    if(bestvz < -15.0 || bestvz>15.0) return;

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

        etaDist->Fill( trk.eta() );

        nTracks++;
        
  } 

  if( doGenParticle_ ){

      edm::Handle<reco::GenParticleCollection> genParticleCollection;
      iEvent.getByLabel(genParticleSrc_, genParticleCollection);

      for(unsigned it=0; it<genParticleCollection->size(); ++it) {

        const reco::GenParticle & genCand = (*genParticleCollection)[it];
        int id = genCand.pdgId();
        int status = genCand.status();
        double rpy_cm = genCand.rapidity()-0.47;

      if ( rpy_cm < -2.87 || rpy_cm > 1.93 ) continue;

      if ( status == 1 ){

        if( id == 310 ){

            genKS_underlying->Fill(rpy_cm, genCand.pt(), genCand.mass());
        }

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

          if (TMath::Abs(mid) != 3322 && TMath::Abs(mid) != 3312 && TMath::Abs(mid) != 3324 && TMath::Abs(mid) != 3314 && TMath::Abs(mid) != 3334){

            genLA_underlying->Fill( rpy_cm, genCand.pt(), genCand.mass() );
          }
        }
     
       }
      }

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

  edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates_xi;
    iEvent.getByLabel(generalV0_xi_,v0candidates_xi);
    if(!v0candidates_xi.isValid()) return;
 
/*
Filling multiplicity into 1D histogram:
 */

multiDist->Fill( nTracks );

//multiplicity bins:
//

if ( nTracks > multmin_ && nTracks < multmax_ ){

    for(unsigned it=0; it<v0candidates_ks->size(); ++it){     
    
            const reco::VertexCompositeCandidate & trk = (*v0candidates_ks)[it];
            const reco:: Candidate * d1 = trk.daughter(0);
            const reco:: Candidate * d2 = trk.daughter(1); 

            auto dau1 = d1->get<reco::TrackRef>();
            auto dau2 = d2->get<reco::TrackRef>();
     
            double eta_dau1 = d1->eta();
            double phi_dau1 = d1->phi();
            double pt_dau1 = d1->pt();
            double px_dau1 = d1->px();
            double py_dau1 = d1->py();
            double pz_dau1 = d1->pz();
            double charge_dau1 = d1->charge();  
            double p_dau1 = d1->p();        

            double eta_dau2 = d2->eta();
            double phi_dau2 = d2->phi();
            double pt_dau2 = d2->pt();
            double px_dau2 = d2->px();
            double py_dau2 = d2->py();
            double pz_dau2 = d2->pz();    
            double charge_dau2 = d2->charge();
            double p_dau2 = d2->p();

            double ks_mass = trk.mass();
            double ks_pt = trk.pt();
            double ks_px = trk.px();
            double ks_py = trk.py();
            double ks_pz = trk.pz();
            double ks_eta = trk.eta();
            double ks_phi = trk.phi();
            double ks_y = trk.rapidity();

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

            ks_y = -ks_y - 0.47;

            ks_ptRapidity->Fill(ks_y,ks_pt);

            //InvMass_ks_underlying->Fill(ks_eta,ks_pt,ks_mass);
            
            if (dau1_Nhits > 3 && dau2_Nhits > 3 && ks_y > -2.87 && ks_y < 1.93 && dlos > 5 && agl > 0.999 && TMath::Abs(dzos1) > 1 && 
              TMath::Abs(dzos2) > 1 && TMath::Abs(dxyos1) > 1 && TMath::Abs(dxyos2) > 1)
            {

              double temp = Mass_ks(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
              double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
              double temp_reverse = Mass_ks(px_dau2,py_dau2,pz_dau2,px_dau1,py_dau1,pz_dau1);
                  if ( (temp < 1.125683 && temp > 1.105683) )continue;
                  if ((temp_reverse < 1.125683 && temp_reverse > 1.105683)) continue;
                  if ( temp_e < 0.015) continue;

                  InvMass_ks_underlying->Fill(ks_y,ks_pt,ks_mass);

                
            }

        }


    for(unsigned it=0; it<v0candidates_la->size(); ++it){     
    
            const reco::VertexCompositeCandidate & trk = (*v0candidates_la)[it];
            const reco:: Candidate * d1 = trk.daughter(0);
            const reco:: Candidate * d2 = trk.daughter(1); 

            auto dau1 = d1->get<reco::TrackRef>();
            auto dau2 = d2->get<reco::TrackRef>();
     
            double eta_dau1 = d1->eta();
            double phi_dau1 = d1->phi();
            double pt_dau1 = d1->pt();
            double px_dau1 = d1->px();
            double py_dau1 = d1->py();
            double pz_dau1 = d1->pz();
            double charge_dau1 = d1->charge();  
            double p_dau1 = d1->p();        

            double eta_dau2 = d2->eta();
            double phi_dau2 = d2->phi();
            double pt_dau2 = d2->pt();
            double px_dau2 = d2->px();
            double py_dau2 = d2->py();
            double pz_dau2 = d2->pz();    
            double charge_dau2 = d2->charge();
            double p_dau2 = d2->p();

            double la_mass = trk.mass();
            double la_pt = trk.pt();
            double la_px = trk.px();
            double la_py = trk.py();
            double la_pz = trk.pz();
            double la_eta = trk.eta();
            double la_phi = trk.phi();
            double la_p = trk.p();
            double la_y = trk.rapidity();

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

            la_y = -la_y - 0.47;

            la_ptRapidity->Fill(la_y,la_pt);

            //InvMass_la_underlying->Fill(la_eta,la_pt,la_mass);

            
            if (dau1_Nhits > 3 && dau2_Nhits > 3 && la_y > -2.87 && la_y < 1.93 && dlos > 5 && agl > 0.999 && TMath::Abs(dzos1) > 1 && 
              TMath::Abs(dzos2) > 1 && TMath::Abs(dxyos1) > 1 && TMath::Abs(dxyos2) > 1)
            {

              double temp = Mass_la(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
              double temp_e = Mass_e(px_dau1,py_dau1,pz_dau1,px_dau2,py_dau2,pz_dau2);
              if ( (temp < 0.517614 && temp > 0.477614) ) continue;
                  if ( temp_e < 0.015) continue;

                  InvMass_la_underlying->Fill(la_y,la_pt,la_mass);

                  for(unsigned it=0; it<v0candidates_xi->size(); ++it){

                    const reco::VertexCompositeCandidate & trk = (*v0candidates_xi)[it];
                    const reco:: Candidate * d1 = trk.daughter(0);
                    const reco:: Candidate * d2 = trk.daughter(1);

                    //PAngle
                    double secvz = trk.vz();
                    double secvx = trk.vx();
                    double secvy = trk.vy();
                    TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
                    TVector3 secvec(trk.px(),trk.py(),trk.pz());

                    double agl = cos(secvec.Angle(ptosvec));

                    double mass = d1->mass();
                    double pt1 = d1->pt();

                    if ( mass == la_mass && pt1 == la_pt ){

                      if ( agl > 0.999 ){

                        xiMass->Fill( trk.mass() );

                          if ( trk.mass() > 1.31486 && trk.mass() < 1.33486 ){

                            XiDaughter->Fill(la_y,la_pt,la_mass);

                    }
                      }
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
    
  TH3D::SetDefaultSumw2();

  InvMass_ks_underlying = fs->make<TH3D>("InvMass_ks_underlying",";eta;pT(GeV/c);mass(GeV/c^{2})",700,-3.5,3.5,120,0,12,360,0.44,0.56);
  InvMass_la_underlying = fs->make<TH3D>("InvMass_la_underlying",";eta;pT(GeV/c);mass(GeV/c^{2})",700,-3.5,3.5,120,0,12,360,1.08,1.16);
  
  if(doGenParticle_){

    genKS_underlying = fs->make<TH3D>("genKS_underlying",";eta;pT(GeV/c);mass(GeV/c^{2})",700,-3.5,3.5,120,0,12,360,0.44,0.56);
    genLA_underlying = fs->make<TH3D>("genLA_underlying",";eta;pT(GeV/c);mass(GeV/c^{2})",700,-3.5,3.5,120,0,12,360,1.08,1.16);
  
  }
  
  XiDaughter = fs->make<TH3D>("XiDaughter",";eta;pT(GeV/c);mass(GeV/c^{2})",700,-3.5,3.5,120,0,12,360,1.08,1.16);

  vertexDistZ = fs->make<TH1D>("vertexDistZ",";Vz;#Events",100,-15,15);
  multiDist = fs->make<TH1D>("multiDist",";mult;#Events",300,0,300);
  xiMass = fs->make<TH1D>("xiMass",";mass",360,1.28,1.36);
  ks_ptRapidity = fs->make<TH2D>("ks_ptRapidity",";rapidity;pT(GeV/c)",700,-3.5,3.5,120,0,12);
  la_ptRapidity = fs->make<TH2D>("la_ptRapidity",";rapidity;pT(GeV/c)",700,-3.5,3.5,120,0,12);
  etaDist = fs->make<TH1D>("etaDist",";eta",600,-3,3);


}

// ------------ method called once each job just after ending the event loop  ------------
void 
V0AnalyzerHisto::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(V0AnalyzerHisto);