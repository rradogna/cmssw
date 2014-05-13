#ifndef GEMCSCSegment_GEMCSCSegmentProducer_h
#define GEMCSCSegment_GEMCSCSegmentProducer_h

/** \class GEMCSCSegmentProducer 
 * Produces a collection of GEM-CSCSegment's in endcap muon. 
 *
 * $Date: 2014/02/06 12:19:20 $
 *
 * \author Raffaella Radogna
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "RecoLocalMuon/GEMRecHit/interface/CSCSegtoGEM.h"

class GEMCSCSegmentBuilder; 

class GEMCSCSegmentProducer : public edm::EDProducer {
public:
    /// Constructor
    explicit GEMCSCSegmentProducer(const edm::ParameterSet&);
    /// Destructor
    ~GEMCSCSegmentProducer();
    /// Produce the GEM-CSCSegment collection
    virtual void produce(edm::Event&, const edm::EventSetup&);

private:
    int iev; // events through
    edm::InputTag inputObjectsTagCSC;
    edm::InputTag inputObjectsTagGEM;// input tag labelling rechits for input
    GEMCSCSegmentBuilder* segmentBuilder_;
    ObjectMapCSC* TheObjectCSC;
};

#endif
