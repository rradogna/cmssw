// -*- C++ -*-
//
// Package:    TestGEMCSCSegmentAnalyzer
// Class:      TestGEMCSCSegmentAnalyzer
// 
/**\class TestGEMCSCSegmentAnalyzer TestGEMCSCSegmentAnalyzer.cc MyAnalyzers/TestGEMCSCSegmentAnalyzer/src/TestGEMCSCSegmentAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Raffaella Radogna

// system include files
#include <memory>
#include <fstream>
#include <sys/time.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>


// root include files
#include "TFile.h"
#include "TH1F.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include <DataFormats/GEMRecHit/interface/GEMCSCSegmentCollection.h>

#include <Geometry/CSCGeometry/interface/CSCChamberSpecs.h>//?
#include <Geometry/CSCGeometry/interface/CSCChamber.h>//?
#include <Geometry/CSCGeometry/interface/CSCLayer.h>//?
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include <Geometry/GEMGeometry/interface/GEMGeometry.h>
#include <Geometry/GEMGeometry/interface/GEMEtaPartition.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <DataFormats/MuonDetId/interface/GEMDetId.h>

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
//
// class declaration
//

class TestGEMCSCSegmentAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TestGEMCSCSegmentAnalyzer(const edm::ParameterSet&);
      ~TestGEMCSCSegmentAnalyzer();



   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      //virtual void endRun(edm::Run const&, edm::EventSetup const&);
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
  //edm::ESHandle<GEMCSCGeometry> gemcscGeom;
    edm::ESHandle<GEMGeometry> gemGeom;
    edm::ESHandle<CSCGeometry> cscGeom;
    

  std::string rootFileName;
  TFile * outputfile;
  TH1F * CSC_fitchi2;
  TH1F * GEMCSC_fitchi2;
  TH1F * GEMCSC_fitchi2_odd;
  TH1F * GEMCSC_fitchi2_even;
  TH1F * GEMCSC_NumGEMRH;
  TH1F * GEMCSC_NumCSCRH;
  TH1F * GEMCSC_NumGEMCSCRH;
  TH1F * GEMCSC_NumGEMCSCSeg;
  TH1F * GEMCSC_Residuals_x;
  TH1F * GEMCSC_Residuals_gem_x;
  TH1F * GEMCSC_Residuals_csc_x;
  TH1F * GEMCSC_Residuals_gem_odd_x;
  TH1F * GEMCSC_Residuals_csc_odd_x;
  TH1F * GEMCSC_Residuals_gem_even_x;
  TH1F * GEMCSC_Residuals_csc_even_x;
  TH1F * GEMCSC_Residuals_cscl1_x;
  TH1F * GEMCSC_Residuals_cscl2_x;
  TH1F * GEMCSC_Residuals_cscl3_x;
  TH1F * GEMCSC_Residuals_cscl4_x;
  TH1F * GEMCSC_Residuals_cscl5_x;
  TH1F * GEMCSC_Residuals_cscl6_x;
  TH1F * GEMCSC_Residuals_geml1_x;
  TH1F * GEMCSC_Residuals_geml2_x;
  TH1F * GEMCSC_Pool_x;
  TH1F * GEMCSC_Pool_gem_x;
  TH1F * GEMCSC_Pool_csc_x;
  TH1F * GEMCSC_Pool_gem_odd_x;
  TH1F * GEMCSC_Pool_csc_odd_x;
  TH1F * GEMCSC_Pool_gem_even_x;
  TH1F * GEMCSC_Pool_csc_even_x;
  TH1F * GEMCSC_Pool_cscl1_x;
  TH1F * GEMCSC_Pool_cscl2_x;
  TH1F * GEMCSC_Pool_cscl3_x;
  TH1F * GEMCSC_Pool_cscl4_x;
  TH1F * GEMCSC_Pool_cscl5_x;
  TH1F * GEMCSC_Pool_cscl6_x;
  TH1F * GEMCSC_Pool_geml1_x;
  TH1F * GEMCSC_Pool_geml2_x;
  TH1F * GEMCSC_Residuals_y;
  TH1F * GEMCSC_Residuals_gem_y;
  TH1F * GEMCSC_Residuals_csc_y;
  TH1F * GEMCSC_Residuals_gem_odd_y;
  TH1F * GEMCSC_Residuals_csc_odd_y;
  TH1F * GEMCSC_Residuals_gem_even_y;
  TH1F * GEMCSC_Residuals_csc_even_y;
  TH1F * GEMCSC_Residuals_cscl1_y;
  TH1F * GEMCSC_Residuals_cscl2_y;
  TH1F * GEMCSC_Residuals_cscl3_y;
  TH1F * GEMCSC_Residuals_cscl4_y;
  TH1F * GEMCSC_Residuals_cscl5_y;
  TH1F * GEMCSC_Residuals_cscl6_y;
  TH1F * GEMCSC_Residuals_geml1_y;
  TH1F * GEMCSC_Residuals_geml2_y;
  TH1F * GEMCSC_Pool_y;
  TH1F * GEMCSC_Pool_gem_y;
  TH1F * GEMCSC_Pool_csc_y;
  TH1F * GEMCSC_Pool_gem_odd_y;
  TH1F * GEMCSC_Pool_csc_odd_y;
  TH1F * GEMCSC_Pool_gem_even_y;
  TH1F * GEMCSC_Pool_csc_even_y;
  TH1F * GEMCSC_Pool_cscl1_y;
  TH1F * GEMCSC_Pool_cscl2_y;
  TH1F * GEMCSC_Pool_cscl3_y;
  TH1F * GEMCSC_Pool_cscl4_y;
  TH1F * GEMCSC_Pool_cscl5_y;
  TH1F * GEMCSC_Pool_cscl6_y;
  TH1F * GEMCSC_Pool_geml1_y;
  TH1F * GEMCSC_Pool_geml2_y;
};

//
// constants, enums and typedefs
//
// constructors and destructor
//
TestGEMCSCSegmentAnalyzer::TestGEMCSCSegmentAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
  rootFileName  = iConfig.getUntrackedParameter<std::string>("RootFileName");

  outputfile = new TFile(rootFileName.c_str(), "RECREATE" );
  CSC_fitchi2 = new TH1F("ReducedChi2_csc","ReducedChi2_csc",160,0.,4.);
    
  GEMCSC_fitchi2 = new TH1F("ReducedChi2_gemcsc","ReducedChi2_gemcsc",160,0.,4.);
  GEMCSC_fitchi2_odd = new TH1F("ReducedChi2_odd_gemcsc","ReducedChi2_odd_gemcsc",160,0.,4.);
  GEMCSC_fitchi2_even = new TH1F("ReducedChi2_even_gemcsc","ReducedChi2_even_gemcsc",160,0.,4.);
  GEMCSC_NumGEMRH = new TH1F("NumGEMRH","NumGEMRH",20,0.,20);
  GEMCSC_NumCSCRH = new TH1F("NumCSCRH","NumCSCRH",20,0.,20);
  GEMCSC_NumGEMCSCRH = new TH1F("NumGEMCSCRH","NumGEMCSCRH",20,0.,20);
  GEMCSC_NumGEMCSCSeg = new TH1F("NumGMCSCSeg","NumGEMCSCSeg",20,0.,20);

  GEMCSC_Residuals_x    = new TH1F("xGEMCSCRes","xGEMCSCRes",100,-0.5,0.5);
  GEMCSC_Residuals_gem_x    = new TH1F("xGEMRes","xGEMRes",100,-0.5,0.5);
  GEMCSC_Residuals_csc_x    = new TH1F("xCSCRes","xCSCRes",100,-0.5,0.5);
  GEMCSC_Residuals_gem_even_x    = new TH1F("xGEMRes_even","xGEMRes even",100,-0.5,0.5);
  GEMCSC_Residuals_csc_even_x    = new TH1F("xCSCRes_even","xCSCRes even",100,-0.5,0.5);
  GEMCSC_Residuals_gem_odd_x    = new TH1F("xGEMRes_odd","xGEMRes odd",100,-0.5,0.5);
  GEMCSC_Residuals_csc_odd_x    = new TH1F("xCSCRes_odd","xCSCRes odd",100,-0.5,0.5);
  GEMCSC_Residuals_cscl1_x = new TH1F("xGEMCSCRes_cscl1","xGEMCSCRes_cscl1",100,-0.5,0.5);
  GEMCSC_Residuals_cscl2_x = new TH1F("xGEMCSCRes_cscl2","xGEMCSCRes_cscl2",100,-0.5,0.5);
  GEMCSC_Residuals_cscl3_x = new TH1F("xGEMCSCRes_cscl3","xGEMCSCRes_cscl3",100,-0.5,0.5);
  GEMCSC_Residuals_cscl4_x = new TH1F("xGEMCSCRes_cscl4","xGEMCSCRes_cscl4",100,-0.5,0.5);
  GEMCSC_Residuals_cscl5_x = new TH1F("xGEMCSCRes_cscl5","xGEMCSCRes_cscl5",100,-0.5,0.5);
  GEMCSC_Residuals_cscl6_x = new TH1F("xGEMCSCRes_cscl6","xGEMCSCRes_cscl6",100,-0.5,0.5);
  GEMCSC_Residuals_geml1_x = new TH1F("xGEMCSCRes_geml1","xGEMCSCRes_geml1",100,-0.5,0.5);
  GEMCSC_Residuals_geml2_x = new TH1F("xGEMCSCRes_geml2","xGEMCSCRes_geml2",100,-0.5,0.5);
  GEMCSC_Pool_x    = new TH1F("xGEMCSCPool","xGEMCSCPool",100,-5.,5.);
  GEMCSC_Pool_gem_x    = new TH1F("xGEMPool","xGEMPool",100,-5.,5.);
  GEMCSC_Pool_csc_x    = new TH1F("xCSCPool","xCSCPool",100,-5.,5.);
  GEMCSC_Pool_gem_even_x    = new TH1F("xGEMPool_even","xGEMPool even",100,-5.,5.);
  GEMCSC_Pool_csc_even_x    = new TH1F("xCSCPool_even","xCSCPool even",100,-5.,5.);
  GEMCSC_Pool_gem_odd_x    = new TH1F("xGEMPool_odd","xGEMPool odd",100,-5.,5.);
  GEMCSC_Pool_csc_odd_x    = new TH1F("xCSCPool_odd","xCSCPool odd",100,-5.,5.);
  GEMCSC_Pool_cscl1_x = new TH1F("xGEMCSCPool_cscl1","xGEMCSCPool_cscl1",100,-5.,5.);
  GEMCSC_Pool_cscl2_x = new TH1F("xGEMCSCPool_cscl2","xGEMCSCPool_cscl2",100,-5.,5.);
  GEMCSC_Pool_cscl3_x = new TH1F("xGEMCSCPool_cscl3","xGEMCSCPool_cscl3",100,-5.,5.);
  GEMCSC_Pool_cscl4_x = new TH1F("xGEMCSCPool_cscl4","xGEMCSCPool_cscl4",100,-5.,5.);
  GEMCSC_Pool_cscl5_x = new TH1F("xGEMCSCPool_cscl5","xGEMCSCPool_cscl5",100,-5.,5.);
  GEMCSC_Pool_cscl6_x = new TH1F("xGEMCSCPool_cscl6","xGEMCSCPool_cscl6",100,-5.,5.);
  GEMCSC_Pool_geml1_x = new TH1F("xGEMCSCPool_geml1","xGEMCSCPool_geml1",100,-5.,5.);
  GEMCSC_Pool_geml2_x = new TH1F("xGEMCSCPool_geml2","xGEMCSCPool_geml2",100,-5.,5.);
  GEMCSC_Residuals_y    = new TH1F("yGEMCSCRes","yGEMCSCRes",100,-5.,5.);
  GEMCSC_Residuals_gem_y    = new TH1F("yGEMRes","yGEMRes",100,-5.,5.);
  GEMCSC_Residuals_csc_y    = new TH1F("yCSCRes","yCSCRes",100,-5.,5.);
  GEMCSC_Residuals_gem_even_y    = new TH1F("yGEMRes_even","yGEMRes even",100,-5.,5.);
  GEMCSC_Residuals_csc_even_y    = new TH1F("yCSCRes_even","yCSCRes even",100,-5.,5.);
  GEMCSC_Residuals_gem_odd_y    = new TH1F("yGEMRes_odd","yGEMRes odd",100,-5.,5.);
  GEMCSC_Residuals_csc_odd_y    = new TH1F("yCSCRes_odd","yCSCRes odd",100,-5.,5.);
  GEMCSC_Residuals_cscl1_y = new TH1F("yGEMCSCRes_cscl1","yGEMCSCRes_cscl1",100,-5.,5.);
  GEMCSC_Residuals_cscl2_y = new TH1F("yGEMCSCRes_cscl2","yGEMCSCRes_cscl2",100,-5.,5.);
  GEMCSC_Residuals_cscl3_y = new TH1F("yGEMCSCRes_cscl3","yGEMCSCRes_cscl3",100,-5.,5.);
  GEMCSC_Residuals_cscl4_y = new TH1F("yGEMCSCRes_cscl4","yGEMCSCRes_cscl4",100,-5.,5.);
  GEMCSC_Residuals_cscl5_y = new TH1F("yGEMCSCRes_cscl5","yGEMCSCRes_cscl5",100,-5.,5.);
  GEMCSC_Residuals_cscl6_y = new TH1F("yGEMCSCRes_cscl6","yGEMCSCRes_cscl6",100,-5.,5.);
  GEMCSC_Residuals_geml1_y = new TH1F("yGEMCSCRes_geml1","yGEMCSCRes_geml1",100,-5.,5.);
  GEMCSC_Residuals_geml2_y = new TH1F("yGEMCSCRes_geml2","yGEMCSCRes_geml2",100,-5.,5.);
  GEMCSC_Pool_y    = new TH1F("yGEMCSCPool","yGEMCSCPool",100,-5.,5.);
  GEMCSC_Pool_gem_y    = new TH1F("yGEMPool","yGEMPool",100,-5.,5.);
  GEMCSC_Pool_csc_y    = new TH1F("yCSCPool","yCSCPool",100,-5.,5.);
  GEMCSC_Pool_gem_even_y    = new TH1F("yGEMPool_even","yGEMPool even",100,-5.,5.);
  GEMCSC_Pool_csc_even_y    = new TH1F("yCSCPool_even","yCSCPool even",100,-5.,5.);
  GEMCSC_Pool_gem_odd_y    = new TH1F("yGEMPool_odd","yGEMPool odd",100,-5.,5.);
  GEMCSC_Pool_csc_odd_y    = new TH1F("yCSCPool_odd","yCSCPool odd",100,-5.,5.);
  GEMCSC_Pool_cscl1_y = new TH1F("yGEMCSCPool_cscl1","yGEMCSCPool_cscl1",100,-5.,5.);
  GEMCSC_Pool_cscl2_y = new TH1F("yGEMCSCPool_cscl2","yGEMCSCPool_cscl2",100,-5.,5.);
  GEMCSC_Pool_cscl3_y = new TH1F("yGEMCSCPool_cscl3","yGEMCSCPool_cscl3",100,-5.,5.);
  GEMCSC_Pool_cscl4_y = new TH1F("yGEMCSCPool_cscl4","yGEMCSCPool_cscl4",100,-5.,5.);
  GEMCSC_Pool_cscl5_y = new TH1F("yGEMCSCPool_cscl5","yGEMCSCPool_cscl5",100,-5.,5.);
  GEMCSC_Pool_cscl6_y = new TH1F("yGEMCSCPool_cscl6","yGEMCSCPool_cscl6",100,-5.,5.);
  GEMCSC_Pool_geml1_y = new TH1F("yGEMCSCPool_geml1","yGEMCSCPool_geml1",100,-5.,5.);
  GEMCSC_Pool_geml2_y = new TH1F("yGEMCSCPool_geml2","yGEMCSCPool_geml2",100,-5.,5.);
}


TestGEMCSCSegmentAnalyzer::~TestGEMCSCSegmentAnalyzer()
{
  CSC_fitchi2->Write();
  GEMCSC_fitchi2->Write();
  GEMCSC_fitchi2_odd->Write();
  GEMCSC_fitchi2_even->Write();
  GEMCSC_NumGEMRH->Write();
  GEMCSC_NumCSCRH->Write();
  GEMCSC_NumGEMCSCRH->Write();
  GEMCSC_NumGEMCSCSeg->Write();
  GEMCSC_Residuals_x->Write();
  GEMCSC_Residuals_gem_x->Write();
  GEMCSC_Residuals_csc_x->Write();
  GEMCSC_Residuals_gem_even_x->Write();
  GEMCSC_Residuals_csc_even_x->Write();
  GEMCSC_Residuals_gem_odd_x->Write();
  GEMCSC_Residuals_csc_odd_x->Write();
  GEMCSC_Residuals_cscl1_x->Write();
  GEMCSC_Residuals_cscl2_x->Write();
  GEMCSC_Residuals_cscl3_x->Write();
  GEMCSC_Residuals_cscl4_x->Write();
  GEMCSC_Residuals_cscl5_x->Write();
  GEMCSC_Residuals_cscl6_x->Write();
  GEMCSC_Residuals_geml1_x->Write();
  GEMCSC_Residuals_geml2_x->Write();
  GEMCSC_Pool_x->Write();
  GEMCSC_Pool_gem_x->Write();
  GEMCSC_Pool_csc_x->Write();
  GEMCSC_Pool_gem_even_x->Write();
  GEMCSC_Pool_csc_even_x->Write();
  GEMCSC_Pool_gem_odd_x->Write();
  GEMCSC_Pool_csc_odd_x->Write();
  GEMCSC_Pool_cscl1_x->Write();
  GEMCSC_Pool_cscl2_x->Write();
  GEMCSC_Pool_cscl3_x->Write();
  GEMCSC_Pool_cscl4_x->Write();
  GEMCSC_Pool_cscl5_x->Write();
  GEMCSC_Pool_cscl6_x->Write();
  GEMCSC_Pool_geml1_x->Write();
  GEMCSC_Pool_geml2_x->Write();
  GEMCSC_Residuals_y->Write();
  GEMCSC_Residuals_gem_y->Write();
  GEMCSC_Residuals_csc_y->Write();
  GEMCSC_Residuals_gem_even_y->Write();
  GEMCSC_Residuals_csc_even_y->Write();
  GEMCSC_Residuals_gem_odd_y->Write();
  GEMCSC_Residuals_csc_odd_y->Write();
  GEMCSC_Residuals_cscl1_y->Write();
  GEMCSC_Residuals_cscl2_y->Write();
  GEMCSC_Residuals_cscl3_y->Write();
  GEMCSC_Residuals_cscl4_y->Write();
  GEMCSC_Residuals_cscl5_y->Write();
  GEMCSC_Residuals_cscl6_y->Write();
  GEMCSC_Residuals_geml1_y->Write();
  GEMCSC_Residuals_geml2_y->Write();
  GEMCSC_Pool_y->Write();
  GEMCSC_Pool_gem_y->Write();
  GEMCSC_Pool_csc_y->Write();
  GEMCSC_Pool_gem_even_y->Write();
  GEMCSC_Pool_csc_even_y->Write();
  GEMCSC_Pool_gem_odd_y->Write();
  GEMCSC_Pool_csc_odd_y->Write();
  GEMCSC_Pool_cscl1_y->Write();
  GEMCSC_Pool_cscl2_y->Write();
  GEMCSC_Pool_cscl3_y->Write();
  GEMCSC_Pool_cscl4_y->Write();
  GEMCSC_Pool_cscl5_y->Write();
  GEMCSC_Pool_cscl6_y->Write();
  GEMCSC_Pool_geml1_y->Write();
  GEMCSC_Pool_geml2_y->Write();
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TestGEMCSCSegmentAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    iSetup.get<MuonGeometryRecord>().get(gemGeom);
    iSetup.get<MuonGeometryRecord>().get(cscGeom);
    const CSCGeometry* cscGeom_ = &*cscGeom;
 
 // ================
 // Sim Tracks
 // ================
    edm::Handle<edm::SimTrackContainer> simTracks;
    iEvent.getByLabel("g4SimHits",simTracks);
    edm::SimTrackContainer::const_iterator simTrack;
    

    
  // ================
  // CSC Segments
  // ================
  edm::Handle<CSCSegmentCollection> cscSegment;
  iEvent.getByLabel("cscSegments","",cscSegment);

  std::cout <<"CSC: Number of Segments "<<cscSegment->size()<<std::endl;
  for (auto cscs = cscSegment->begin(); cscs != cscSegment->end(); cscs++) {

    CSCDetId CSCId = cscs->cscDetId();
    if(!(CSCId.station() == 1 && CSCId.ring() == 1)) continue;
    std::cout <<"CSC: Original CSCDetID "<<CSCId<<std::endl;
    //const CSCChamber* cscChamber = cscGeom_->chamber(CSCId);
    //std::cout<<cscChamber<<std::endl;
    auto cscrhs = cscs->specificRecHits();

    std::cout <<"CSC: Number of RecHits "<<cscrhs.size()<<std::endl;

    std::cout <<"CSC: Chi2: "<<cscs->chi2()<<" ndof "<<cscs->degreesOfFreedom()<<std::endl;

  }

  std::cout <<"------------------------------------------------------------------------------"<<std::endl;

  // ================
  // GEMCSC Segments
  // ================
  edm::Handle<GEMCSCSegmentCollection> gemcscSegment;
  iEvent.getByLabel("gemcscSegments","",gemcscSegment);

  if(gemcscSegment->size()!=0)GEMCSC_NumGEMCSCSeg->Fill(gemcscSegment->size());
  std::cout <<"GEMCSC: Number of Segments "<<gemcscSegment->size()<<std::endl;


    for (auto gemcscs = gemcscSegment->begin(); gemcscs != gemcscSegment->end(); gemcscs++) {
      std::cout<<"++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
      std::cout<<"GEMCSC segment"<<std::endl;
      
      auto gemrhs_if = gemcscs->gemRecHits();
      if(gemrhs_if.size()==0)continue;
        
        
        ///////// GEMCSC seg //////////////////////////////////////
        CSCDetId gemcscId = gemcscs->cscDetId();
        std::cout <<" Original GEMCSCDetID "<<gemcscId<<std::endl;
        const CSCChamber* cscChamber = cscGeom_->chamber(gemcscId);
        std::cout<<cscChamber<<std::endl;
        //const CSCChamber* chamber = cscGeom->chamber(gemcscId); // common reference //non funziona, da modificare in data formats?
        //const CSCChamber* chamber = cscGeom->chamber(id.rawId()); // common reference

        auto gemcscsegLP = gemcscs->localPosition();
        auto gemcscsegLD = gemcscs->localDirection();
        std::cout <<" Local Position x = "<<gemcscsegLP.x()<<" y= "<<gemcscsegLP.y()<<" z= "<<gemcscsegLP.z()<<std::endl;
        std::cout <<" Local Direction theta = "<<gemcscsegLD.theta()<<" phi="<<gemcscsegLD.phi()<<std::endl;
        std::cout <<"GEMCSC: Chi2: "<<gemcscs->chi2()<<" ndof "<<gemcscs->degreesOfFreedom()<<" chi2/ndof "<< gemcscs->chi2()/gemcscs->degreesOfFreedom()<<std::endl;
        GEMCSC_fitchi2->Fill(gemcscs->chi2()/gemcscs->degreesOfFreedom());
        
        CSCDetId id((*gemcscs).geographicalId());
        int chamber = id.chamber();
        if(chamber%2!=0){ //std::cout<<"camera dispari"<<chamber<<std::endl;
          GEMCSC_fitchi2_odd->Fill(gemcscs->chi2()/gemcscs->degreesOfFreedom());}
        else{GEMCSC_fitchi2_even->Fill(gemcscs->chi2()/gemcscs->degreesOfFreedom());}
      
        GEMCSC_NumGEMCSCRH->Fill(gemcscs->cscSegment().specificRecHits().size()+gemcscs->gemRecHits().size());

      
        /*auto roll = cscGeom->etaPartition(id);
         std::cout <<"Global Segment Position "<<  roll->toGlobal(gemcscs->localPosition())<<std::endl;
         auto segLP = gemcscs->localPosition();
         auto segLD = gemcscs->localDirection();
         std::cout <<" Global Direction theta = "<<segLD.theta()<<" phi="<<segLD.phi()<<std::endl;*/
      
        ///////// CSC seg /////////////////////////////////////
        CSCSegment cscSeg = gemcscs->cscSegment();
        std::cout<<"-----------------CSC segment"<<std::endl;
        CSCDetId CSCId_new = cscSeg.cscDetId();
        //if(!(CSCId_new.station() == 1 && CSCId_new.ring() == 1)) continue;
        std::cout <<"GEMCSC: CSC segment CSCDetID "<<CSCId_new<<std::endl;
        
        auto cscrhs = cscSeg.specificRecHits();
        GEMCSC_NumCSCRH->Fill(cscrhs.size());
        std::cout <<"GEMCSC: CSC segment Number of RecHits "<<cscrhs.size()<<std::endl;
        
        auto cscsegLP = cscSeg.localPosition();
        auto cscsegLD = cscSeg.localDirection();
        std::cout <<" Local Position x = "<<cscsegLP.x()<<" y= "<<cscsegLP.y()<<" z= "<<cscsegLP.z()<<std::endl;
        std::cout <<" Local Direction theta = "<<cscsegLD.theta()<<" phi="<<cscsegLD.phi()<<std::endl;
        std::cout <<"GEMCSC: CSC segment Chi2 "<<cscSeg.chi2()<<" ndof "<<cscSeg.degreesOfFreedom()<<" chi2/ndof "<<cscSeg.chi2()/cscSeg.degreesOfFreedom()<<std::endl;
        CSC_fitchi2->Fill(cscSeg.chi2()/cscSeg.degreesOfFreedom());
      
        ////////////////////////
        for (auto rh = cscrhs.begin(); rh!= cscrhs.end(); rh++){
            std::cout<<"-----------------CSC RecHits"<<std::endl;
            //CSCDetId cscrhId = rh.cscDetId();
            CSCDetId cscrhId = (CSCDetId)(*rh).cscDetId();
            const CSCLayer* cscrhRef = cscGeom_->layer( cscrhId );
          
            auto cscrhLP = rh->localPosition();
            auto cscrhLEP = rh->localPositionError();
            std::cout <<" Local Position x = "<<cscrhLP.x()<<" y= "<<cscrhLP.y()<<" z= "<<cscrhLP.z()<<std::endl;
            auto cscrhGP = cscrhRef->toGlobal(cscrhLP);
            auto cscrhLP_inSegmRef = cscChamber->toLocal(cscrhGP);
            float xe  = gemcscsegLP.x()+gemcscsegLD.x()*cscrhLP_inSegmRef.z()/gemcscsegLD.z();
            float ye  = gemcscsegLP.y()+gemcscsegLD.y()*cscrhLP_inSegmRef.z()/gemcscsegLD.z();
            float ze = cscrhLP_inSegmRef.z();
            LocalPoint extrPoint(xe,ye,ze); // in segment rest frame
          
            auto extSegm = cscrhRef->toLocal(cscChamber->toGlobal(extrPoint)); // in layer restframe
          
            std::cout <<" CSC Layer Id "<<cscrhId<<"  error on the local point "<<  cscrhLEP
            <<"\n-> Ensamble Rest Frame  RH local  position "<<cscrhLP_inSegmRef<<"  Segment extrapolation "<<extrPoint
            <<"\n-> Layer Rest Frame  RH local  position "<<cscrhLP<<"  Segment extrapolation "<<extSegm
            <<std::endl;
            GEMCSC_Residuals_x->Fill(cscrhLP.x()-extSegm.x());
            GEMCSC_Residuals_y->Fill(cscrhLP.y()-extSegm.y());
            GEMCSC_Pool_x->Fill((cscrhLP.x()-extSegm.x())/sqrt(cscrhLEP.xx()));
            GEMCSC_Pool_y->Fill((cscrhLP.y()-extSegm.y())/sqrt(cscrhLEP.yy()));
            GEMCSC_Residuals_csc_x->Fill(cscrhLP.x()-extSegm.x());
            GEMCSC_Residuals_csc_y->Fill(cscrhLP.y()-extSegm.y());
            GEMCSC_Pool_csc_x->Fill((cscrhLP.x()-extSegm.x())/sqrt(cscrhLEP.xx()));
            GEMCSC_Pool_csc_y->Fill((cscrhLP.y()-extSegm.y())/sqrt(cscrhLEP.yy()));

            CSCDetId id((*rh).geographicalId());
            int chamber = id.chamber();
            if(chamber%2!=0){
              //std::cout<<"camera dispari"<<chamber<<std::endl;
              GEMCSC_Residuals_csc_odd_x->Fill(cscrhLP.x()-extSegm.x());
              GEMCSC_Residuals_csc_odd_y->Fill(cscrhLP.y()-extSegm.y());
              GEMCSC_Pool_csc_odd_x->Fill((cscrhLP.x()-extSegm.x())/sqrt(cscrhLEP.xx()));
              GEMCSC_Pool_csc_odd_y->Fill((cscrhLP.y()-extSegm.y())/sqrt(cscrhLEP.yy()));
            }
            else {
              GEMCSC_Residuals_csc_even_x->Fill(cscrhLP.x()-extSegm.x());
              GEMCSC_Residuals_csc_even_y->Fill(cscrhLP.y()-extSegm.y());
              GEMCSC_Pool_csc_even_x->Fill((cscrhLP.x()-extSegm.x())/sqrt(cscrhLEP.xx()));
              GEMCSC_Pool_csc_even_y->Fill((cscrhLP.y()-extSegm.y())/sqrt(cscrhLEP.yy()));

            }
          
            switch (cscrhId.layer()){
              case 1:
                  std::cout<<"l1"<<std::endl;
                  GEMCSC_Residuals_cscl1_x->Fill(cscrhLP.x()-extSegm.x());
                  GEMCSC_Residuals_cscl1_y->Fill(cscrhLP.y()-extSegm.y());
                  GEMCSC_Pool_cscl1_x->Fill((cscrhLP.x()-extSegm.x())/sqrt(cscrhLEP.xx()));
                  GEMCSC_Pool_cscl1_y->Fill((cscrhLP.x()-extSegm.y())/sqrt(cscrhLEP.yy()));
                  break;
              case 2:
                  std::cout<<"l2"<<std::endl;
                  GEMCSC_Residuals_cscl2_x->Fill(cscrhLP.x()-extSegm.x());
                  GEMCSC_Residuals_cscl2_y->Fill(cscrhLP.y()-extSegm.y());
                  GEMCSC_Pool_cscl2_x->Fill((cscrhLP.x()-extSegm.x())/sqrt(cscrhLEP.xx()));
                  GEMCSC_Pool_cscl2_y->Fill((cscrhLP.x()-extSegm.y())/sqrt(cscrhLEP.yy()));
                  break;
                  
              case 3:
                  std::cout<<"l3"<<std::endl;
                  GEMCSC_Residuals_cscl3_x->Fill(cscrhLP.x()-extSegm.x());
                  GEMCSC_Residuals_cscl3_y->Fill(cscrhLP.y()-extSegm.y());
                  GEMCSC_Pool_cscl3_x->Fill((cscrhLP.x()-extSegm.x())/sqrt(cscrhLEP.xx()));
                  GEMCSC_Pool_cscl3_y->Fill((cscrhLP.x()-extSegm.y())/sqrt(cscrhLEP.yy()));
                  break;
              case 4:
                  std::cout<<"l4"<<std::endl;
                  GEMCSC_Residuals_cscl4_x->Fill(cscrhLP.x()-extSegm.x());
                  GEMCSC_Residuals_cscl4_y->Fill(cscrhLP.y()-extSegm.y());
                  GEMCSC_Pool_cscl4_x->Fill((cscrhLP.x()-extSegm.x())/sqrt(cscrhLEP.xx()));
                  GEMCSC_Pool_cscl4_y->Fill((cscrhLP.x()-extSegm.y())/sqrt(cscrhLEP.yy()));
                  break;
              case 5:
                  std::cout<<"l5"<<std::endl;
                  GEMCSC_Residuals_cscl5_x->Fill(cscrhLP.x()-extSegm.x());
                  GEMCSC_Residuals_cscl5_y->Fill(cscrhLP.y()-extSegm.y());
                  GEMCSC_Pool_cscl5_x->Fill((cscrhLP.x()-extSegm.x())/sqrt(cscrhLEP.xx()));
                  GEMCSC_Pool_cscl5_y->Fill((cscrhLP.x()-extSegm.y())/sqrt(cscrhLEP.yy()));
                  break;
              case 6:
                  std::cout<<"l6"<<std::endl;
                  GEMCSC_Residuals_cscl6_x->Fill(cscrhLP.x()-extSegm.x());
                  GEMCSC_Residuals_cscl6_y->Fill(cscrhLP.y()-extSegm.y());
                  GEMCSC_Pool_cscl6_x->Fill((cscrhLP.x()-extSegm.x())/sqrt(cscrhLEP.xx()));
                  GEMCSC_Pool_cscl6_y->Fill((cscrhLP.x()-extSegm.y())/sqrt(cscrhLEP.yy()));
                  break;
              default:
                  std::cout <<" Unphysical GEMCSC layer "<<cscrhId<<std::endl;
          }
      }
        //////////////////
      
        //////// GEM recHits ////////////////////
        auto gemrhs = gemcscs->gemRecHits();
        std::cout<<"-----------------GEM RecHits"<<std::endl;
        GEMCSC_NumGEMRH->Fill(gemrhs.size());
        //std::cout <<"GEMCSC Ensamble Det Id "<<id<<"  Numbero of RecHits "<<gemcscrhs.size()<<std::endl;
        std::cout <<"GEMCSC: Number of GEM RecHits "<<gemrhs.size()<<std::endl;
        //loop on rechits.... take layer local position -> global -> ensamble local position same frame as segment
      
        for (auto rh = gemrhs.begin(); rh!= gemrhs.end(); rh++){
            //GEMDetId gemid(rh->geographicalId());
            auto gemrhId = rh->gemId();
            const GEMEtaPartition* gemrhRef  = gemGeom->etaPartition(gemrhId);
            auto gemrhLP = rh->localPosition();
            auto gemrhLEP = rh->localPositionError();
            std::cout <<" Local Position x = "<<gemrhLP.x()<<" y= "<<gemrhLP.y()<<" z= "<<gemrhLP.z()<<std::endl;
            auto gemrhGP = gemrhRef->toGlobal(gemrhLP);
            auto gemrhLP_inSegmRef = cscChamber->toLocal(gemrhGP);
            float xe  = gemcscsegLP.x()+gemcscsegLD.x()*gemrhLP_inSegmRef.z()/gemcscsegLD.z();
            float ye  = gemcscsegLP.y()+gemcscsegLD.y()*gemrhLP_inSegmRef.z()/gemcscsegLD.z();
            float ze = gemrhLP_inSegmRef.z();
            LocalPoint extrPoint(xe,ye,ze); // in segment rest frame
          
            auto extSegm = gemrhRef->toLocal(cscChamber->toGlobal(extrPoint)); // in layer restframe

            std::cout <<" GEM Layer Id "<<rh->gemId()<<"  error on the local point "<<  gemrhLEP
            <<"\n-> Ensamble Rest Frame  RH local  position "<<gemrhLP_inSegmRef<<"  Segment extrapolation "<<extrPoint
            <<"\n-> Layer Rest Frame  RH local  position "<<gemrhLP<<"  Segment extrapolation "<<extSegm
            <<std::endl;
            GEMCSC_Residuals_x->Fill(gemrhLP.x()-extSegm.x());
            GEMCSC_Residuals_y->Fill(gemrhLP.y()-extSegm.y());
            GEMCSC_Pool_x->Fill((gemrhLP.x()-extSegm.x())/sqrt(gemrhLEP.xx()));
            GEMCSC_Pool_y->Fill((gemrhLP.y()-extSegm.y())/sqrt(gemrhLEP.yy()));
            GEMCSC_Residuals_gem_x->Fill(gemrhLP.x()-extSegm.x());
            GEMCSC_Residuals_gem_y->Fill(gemrhLP.y()-extSegm.y());
            GEMCSC_Pool_gem_x->Fill((gemrhLP.x()-extSegm.x())/sqrt(gemrhLEP.xx()));
            GEMCSC_Pool_gem_y->Fill((gemrhLP.y()-extSegm.y())/sqrt(gemrhLEP.yy()));
            
            GEMDetId id((*rh).geographicalId());
            int chamber = id.chamber();
            if(chamber%2!=0){
              //std::cout<<"camera dispari"<<chamber<<std::endl;
              GEMCSC_Residuals_gem_odd_x->Fill(gemrhLP.x()-extSegm.x());
              GEMCSC_Residuals_gem_odd_y->Fill(gemrhLP.y()-extSegm.y());
              GEMCSC_Pool_gem_odd_x->Fill((gemrhLP.x()-extSegm.x())/sqrt(gemrhLEP.xx()));
              GEMCSC_Pool_gem_odd_y->Fill((gemrhLP.y()-extSegm.y())/sqrt(gemrhLEP.yy()));
            }
            else {
              GEMCSC_Residuals_gem_even_x->Fill(gemrhLP.x()-extSegm.x());
              GEMCSC_Residuals_gem_even_y->Fill(gemrhLP.y()-extSegm.y());
              GEMCSC_Pool_gem_even_x->Fill((gemrhLP.x()-extSegm.x())/sqrt(gemrhLEP.xx()));
              GEMCSC_Pool_gem_even_y->Fill((gemrhLP.y()-extSegm.y())/sqrt(gemrhLEP.yy()));
              
            }

            switch (gemrhId.layer()){
              case 1:
                  std::cout<<"l1"<<std::endl;
                  GEMCSC_Residuals_geml1_x->Fill(gemrhLP.x()-extSegm.x());
                  GEMCSC_Residuals_geml1_y->Fill(gemrhLP.y()-extSegm.y());
                  GEMCSC_Pool_geml1_x->Fill((gemrhLP.x()-extSegm.x())/sqrt(gemrhLEP.xx()));
                  GEMCSC_Pool_geml1_y->Fill((gemrhLP.x()-extSegm.y())/sqrt(gemrhLEP.yy()));
                  break;
              case 2:
                  std::cout<<"l2"<<std::endl;
                  GEMCSC_Residuals_geml2_x->Fill(gemrhLP.x()-extSegm.x());
                  GEMCSC_Residuals_geml2_y->Fill(gemrhLP.y()-extSegm.y());
                  GEMCSC_Pool_geml2_x->Fill((gemrhLP.x()-extSegm.x())/sqrt(gemrhLEP.xx()));
                  GEMCSC_Pool_geml2_y->Fill((gemrhLP.x()-extSegm.y())/sqrt(gemrhLEP.yy()));
                  break;
              default:
                  std::cout <<" Unphysical GEMCSC layer "<<gemrhId<<std::endl;
            }
        }
      
    
        /////////////////////////////////////////////
        auto gemcscrhs = gemcscs->specificRecHits();
        for (auto rh = gemcscrhs.begin(); rh!= gemcscrhs.end(); rh++){

            if (rh->geographicalId().subdetId() == MuonSubdetId::CSC){

                std::cout<<"CSC found"<<std::endl;

                CSCDetId id(rh->geographicalId());

                //int layer = id.layer();
                int station = id.station();
                int ring = id.ring();
                //int chamber = id.chamber();
                //int roll = id.roll();

                std::cout<<"CSC Region"<<" Station "<<station<<" ring "<<ring<<std::endl;

            }
            else if (rh->geographicalId().subdetId() == MuonSubdetId::GEM){

                //std::cout<<"GEM found"<<std::endl;

                GEMDetId id(rh->geographicalId());

                int region = id.region();
                int layer = id.layer();
                int station = id.station();
                int ring = id.ring();
                int chamber = id.chamber();
                int roll = id.roll();

                std::cout<<"GEM Region"<<region<<" Station "<<station<<" ring "<<ring<<" layer "<<layer<<" chamber "<<chamber<<" roll "<<roll<<std::endl;
            }


        }

      
        auto gemcscsegGD = cscChamber->toGlobal(gemcscs->localPosition());
        for (simTrack = simTracks->begin(); simTrack != simTracks->end(); ++simTrack){
            double simEta = (*simTrack).momentum().eta();
            double simPhi = (*simTrack).momentum().phi();
            double dR = sqrt(pow((simEta-gemcscsegGD.eta()),2) + pow((simPhi-gemcscsegGD.phi()),2));
            if(dR > 0.1) continue;
            std::cout <<"Sim Track: Eta "<<simEta<<" Phi "<<simPhi<<" segment GD eta "<<gemcscsegGD.eta()<<" segment GD phi "<<gemcscsegGD.phi()<<std::endl;
        
        
        }
  
  
  }// loop gemcsc segments

  std::cout <<"------------------------------------------------------------------------------"<<std::endl;
  std::cout <<"------------------------------------------------------------------------------"<<std::endl;
   

    
}
// ------------ method called once each job just before starting event loop  ------------
void 
TestGEMCSCSegmentAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TestGEMCSCSegmentAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TestGEMCSCSegmentAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  //iSetup.get<MuonGeometryRecord>().get(gemGeom);
  //iSetup.get<MuonGeometryRecord>().get(cscGeom);



}

//define this as a plug-in
DEFINE_FWK_MODULE(TestGEMCSCSegmentAnalyzer);
