
/** \file GEMCSCSegmentBuilder.cc
 *
 *
 */

#include <RecoLocalMuon/GEMRecHit/src/GEMCSCSegmentBuilder.h>
#include "FWCore/Framework/interface/ESHandle.h"
#include <Geometry/Records/interface/MuonGeometryRecord.h>

#include <Geometry/CSCGeometry/interface/CSCChamberSpecs.h>//?
#include <Geometry/CSCGeometry/interface/CSCLayer.h>//?
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include <Geometry/GEMGeometry/interface/GEMGeometry.h>
#include <Geometry/GEMGeometry/interface/GEMEtaPartition.h>

#include <DataFormats/MuonDetId/interface/CSCDetId.h>
#include <DataFormats/MuonDetId/interface/GEMDetId.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2D.h>//??
#include <DataFormats/CSCRecHit/interface/CSCRangeMapAccessor.h>//??
#include <DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h>//??
#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>//??
#include <DataFormats/CSCRecHit/interface/CSCSegment.h>//??
#include <DataFormats/GEMRecHit/interface/GEMRecHit.h>
#include <DataFormats/GEMRecHit/interface/GEMRecHitCollection.h>


#include <RecoLocalMuon/GEMRecHit/src/GEMCSCSegmentAlgorithm.h>
#include <RecoLocalMuon/GEMRecHit/src/GEMCSCSegmentBuilderPluginFactory.h>
#include <RecoLocalMuon/GEMRecHit/interface/CSCSegtoGEM.h>

#include <FWCore/Utilities/interface/Exception.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 



GEMCSCSegmentBuilder::GEMCSCSegmentBuilder(const edm::ParameterSet& ps) : gemgeom_(0), cscgeom_(0) {

    // Algo name
    std::string algoName = ps.getParameter<std::string>("algo_name");
        
    LogDebug("GEMCSCSegment|GEMCSC")<< "GEMCSCSegmentBuilder algorithm name: " << algoName;

    // SegAlgo parameter set 	  
    edm::ParameterSet segAlgoPSet = ps.getParameter<edm::ParameterSet>("algo_psets");
 
    // Ask factory to build this algorithm, giving it appropriate ParameterSet
    algo = GEMCSCSegmentBuilderPluginFactory::get()->create(algoName, segAlgoPSet);

}

GEMCSCSegmentBuilder::~GEMCSCSegmentBuilder() {
    delete algo;
}

void GEMCSCSegmentBuilder::build(const GEMRecHitCollection* recHits, const CSCSegmentCollection* cscsegments, const edm::EventSetup& setup, GEMCSCSegmentCollection& oc) {
  	
  LogDebug("GEMCSCSegment|GEMCSC")<< "Total number of rechits in this event: " << recHits->size()<< "Total number of csc segments" << cscsegments->size();

    std::vector<CSCDetId> chambers;//?
    std::vector<CSCDetId>::const_iterator chIt;//?
    edm::ESHandle<GEMGeometry> gemGeo;
    edm::ESHandle<CSCGeometry> cscGeo;
  
    setup.get<MuonGeometryRecord>().get(gemGeo);
    setup.get<MuonGeometryRecord>().get(cscGeo);
    
    //std::map<uint32_t, std::vector<GEMRecHit*> > ensembleRH;
    std::map<uint32_t, std::vector<RecHit2DLocalPos*> > ensembleRH;
    //std::map<uint32_t, std::vector<TrackingRecHit*> > ensembleRH;

for(CSCSegmentCollection::const_iterator segm = cscsegments->begin(); segm != cscsegments->end(); segm++) {   
    CSCDetId CSCId = segm->cscDetId();
    if(CSCId.station()==1 && CSCId.ring()==1 ){
    CSCDetId id(segm->cscDetId().endcap(),segm->cscDetId().station(),segm->cscDetId().ring(),segm->cscDetId().chamber(),1);
    
    auto cscSrhs = segm->specificRecHits();
    for (auto rh = cscSrhs.begin(); rh!= cscSrhs.end(); rh++){
    std::vector<RecHit2DLocalPos* > pp = ensembleRH[id.rawId()];
    pp.push_back(rh->clone());
    ensembleRH[id.rawId()]=pp;
    }
    
    int cscEndCap = CSCId.endcap();
    int cscStation = CSCId.station();
    int cscRing = CSCId.ring();
    int gemRegion = 1; if(cscEndCap==2) gemRegion= -1;//Relacion entre las endcaps
    int gemRing = cscRing;
    //if(cscRing==4)gemRing =1; ////csc ring=4 è me1A, mentre ring 1 è me1B, gemring=1
    int gemStation = cscStation;
    int gemChamber = CSCId.chamber();
    ObjectMapCSC* TheObjectCSC = ObjectMapCSC::GetInstance(setup);
    CSCStationIndex theindex(gemRegion,gemStation,gemRing,gemChamber);
    std::set<GEMDetId> rollsForThisCSC = TheObjectCSC->GetInstance(setup)->GetRolls(theindex);
    
   /* //TrackingRecHit* lp = dynamic_cast< TrackingRecHit* >( *it2 );
    std::vector<RecHit2DLocalPos* > pp = ensembleRH[id.rawId()];
    pp.push_back(it2->clone());
    ensembleRH[id.rawId()]=pp;*/
    
   /* for (std::set<GEMDetId>::iterator iteraRoll = rollsForThisCSC.begin();iteraRoll != rollsForThisCSC.end(); iteraRoll++){
const GEMEtaPartition* rollasociated = gemGeo->etaPartitions(*iteraRoll);
GEMDetId gemId = rollasociated->id();*/

  // Loop on the GEM rechit and select the different GEM Ensemble
  for(GEMRecHitCollection::const_iterator it2 = recHits->begin(); it2 != recHits->end(); it2++) {        
    //GEMDetId id(it2->gemId().region(),it2->gemId().ring(),it2->gemId().station(),1,it2->gemId().chamber(),it2->gemId().roll());

    for (std::set<GEMDetId>::iterator iteraRoll = rollsForThisCSC.begin();iteraRoll != rollsForThisCSC.end(); iteraRoll++){
    	const GEMEtaPartition* rollasociated = gemGeo->etaPartition(*iteraRoll);
    	GEMDetId gemIdfromSegm = rollasociated->id();
    	if(it2->gemId().roll()== gemIdfromSegm.roll()){
    		std::vector<RecHit2DLocalPos* > pp = ensembleRH[id.rawId()];
    		pp.push_back(it2->clone());
    		ensembleRH[id.rawId()]=pp;
		
    	}// if roll of gem rec hit is on of the csc segment projection roll	
    }// for roll of te segment
  }// for rec hits
  }// if ME1/1b
}// for csc segments

    //for(CSCRecHit2DCollection::const_iterator it2 = recHits->begin(); it2 != recHits->end(); it2++) {
    for(auto enIt=ensembleRH.begin(); enIt != ensembleRH.end(); ++enIt) {    
        //std::vector<const GEMRecHit*> gemRecHits;
	std::vector<const RecHit2DLocalPos*> gemcscRecHits;
        std::map<uint32_t,const GEMEtaPartition* > ens;
        //std::map<uint32_t,const CSCChamber* > ens;
    
        const GEMEtaPartition* firstlayer = gemgeom_->etaPartition(enIt->first);
	//const CSCChamber* firstlayer = cscgeom_->chamber(enIt->first);
        for(auto rechit = enIt->second.begin(); rechit != enIt->second.end(); rechit++) {
          gemcscRecHits.push_back(*rechit);
          //ens[(*rechit)->gemId()]=gemgeom_->etaPartition((*rechit)->gemId());
	  ens[(*rechit)->rawId()]=gemgeom_->etaPartition((*rechit)->rawId());
	  //ens[(*rechit)->rawId()]=cscgeom_->chamber((*rechit)->rawId());
        }    
    GEMCSCSegmentAlgorithm::GEMCSCEnsamble ensamble(std::pair<const GEMEtaPartition*, std::map<uint32_t,const GEMEtaPartition *> >(firstlayer,ens));
    
    LogDebug("GEMCSCSegment|GEMCSC") << "found " << gemcscRecHits.size() << " rechits in chamber " << *enIt;

        // given the chamber select the appropriate algo... and run it
    std::vector<GEMCSCSegment> segv = algo->run(ensamble, gemcscRecHits);
    CSCDetId mid(enIt->first);
    LogDebug("GEMSegment|GEM") << "found " << segv.size() << " segments in chamber " << mid;
    
    // Add the segments to master collection
    oc.put(mid, segv.begin(), segv.end());
    }
}


void GEMCSCSegmentBuilder::setGeometry(const GEMGeometry* gemgeom, const CSCGeometry* cscgeom) {
	gemgeom_ = gemgeom;
	cscgeom_ = cscgeom;
}
/*void GEMCSCSegmentBuilder::setSetup(const edm::EventSetup& setup) {
	
	setup_ = setup;
}*/

