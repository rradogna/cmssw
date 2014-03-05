
/** \file GEMCSCSegmentBuilder.cc
 *
 *
 */

#include <RecoLocalMuon/GEMRecHit/src/GEMCSCSegmentBuilder.h>

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

void GEMCSCSegmentBuilder::build(const GEMRecHitCollection* recHits, const CSCSegmentCollection* cscsegments, GEMCSCSegmentCollection& oc) {
  	
  LogDebug("GEMCSCSegment|GEMCSC")<< "Total number of rechits in this event: " << recHits->size()<< "Total number of csc segments" << cscsegments->size();

    std::vector<CSCDetId> chambers;//?
    std::vector<CSCDetId>::const_iterator chIt;//?
    
    //std::map<uint32_t, std::vector<GEMRecHit*> > ensembleRH;
    //std::map<uint32_t, std::vector<RecHit2DLocalPos*> > ensembleRH;
    std::map<uint32_t, std::vector<TrackingRecHit*> > ensembleRH;

  // Loop on the GEM rechit and select the different GEM Ensemble
  for(GEMRecHitCollection::const_iterator it2 = recHits->begin(); it2 != recHits->end(); it2++) {        
    GEMDetId id(it2->gemId().region(),it2->gemId().ring(),it2->gemId().station(),1,it2->gemId().chamber(),it2->gemId().roll());
    TrackingRecHit* lp = dynamic_cast< TrackingRecHit* >( *it2 );
    std::vector<TrackingRecHit* > pp = ensembleRH[id.rawId()];
    pp.push_back(lp->clone());
    ensembleRH[id.rawId()]=pp;
  }

    //for(CSCRecHit2DCollection::const_iterator it2 = recHits->begin(); it2 != recHits->end(); it2++) {
    for(auto enIt=ensembleRH.begin(); enIt != ensembleRH.end(); ++enIt) {    
        std::vector<const GEMRecHit*> gemRecHits;
        std::map<uint32_t,const GEMEtaPartition* > ens;
    
        const GEMEtaPartition* firstlayer = gemgeom_->etaPartition(enIt->first);
        for(auto rechit = enIt->second.begin(); rechit != enIt->second.end(); rechit++) {
          gemRecHits.push_back(*rechit);
          ens[(*rechit)->gemId()]=gemgeom_->etaPartition((*rechit)->gemId());
        }    
    GEMCSCSegmentAlgorithm::GEMCSCEnsamble ensamble(std::pair<const GEMEtaPartition*, std::map<uint32_t,const GEMEtaPartition *> >(firstlayer,ens));
    
    LogDebug("GEMCSCSegment|GEMCSC") << "found " << gemRecHits.size() << " rechits in chamber " << *enIt;

        // given the chamber select the appropriate algo... and run it
    std::vector<GEMCSCSegment> segv = algo->run(ensamble, gemRecHits);
    GEMDetId mid(enIt->first);
    LogDebug("GEMSegment|GEM") << "found " << segv.size() << " segments in chamber " << mid;
    
    // Add the segments to master collection
    oc.put(mid, segv.begin(), segv.end());
    }
}

void GEMCSCSegmentBuilder::setGeometry(const GEMGeometry* gemgeom, const CSCGeometry* cscgeom) {
	gemgeom_ = gemgeom;
	cscgeom_ = cscgeom;
}

