#ifndef GEMCSCSegment_GEMCSCSegmentBuilder_h
#define GEMCSCSegment_GEMCSCSegmentBuilder_h

/** \class GEMCSCSegmentBuilder 
 * Algorithm to build GEMCSCSegment's from GEMRecHit collection and CSCSegment collection
 * by implementing a 'build' function required by GEMCSCSegmentProducer.
 *
 * Implementation notes: <BR>
 * Configured via the Producer's ParameterSet. <BR>
 * Presume this might become an abstract base class one day. <BR>
 *
 * $Date:  $
 * $Revision: 1.3 $
 * \author Raffaella Radogna
 *
 *
 */
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>      //?
#include <Geometry/GEMGeometry/interface/GEMGeometry.h>       //?
#include <Geometry/GEMGeometry/interface/GEMEtaPartition.h>   //?

#include <DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h>  //?
#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>   //?
#include <DataFormats/GEMRecHit/interface/GEMRecHit.h>           //????
#include <DataFormats/GEMRecHit/interface/GEMRecHitCollection.h>

#include <DataFormats/GEMRecHit/interface/GEMCSCSegmentCollection.h>

#include <FWCore/ParameterSet/interface/ParameterSet.h>

class GEMCSCSegmentAlgorithm;

class GEMCSCSegmentBuilder {
public:
   
    /** Configure the algorithm via ctor.
     * Receives ParameterSet percolated down from EDProducer
     * which owns this Builder.
     */
    explicit GEMCSCSegmentBuilder(const edm::ParameterSet&);
    /// Destructor
    ~GEMCSCSegmentBuilder();

    /** Find rechits in each GEMChamber and Segment in CSCChamber, build GEMCSCSegment,
     *  and fill into output collection. 
     */
    void build(const GEMRecHitCollection* rechits, GEMCSCSegmentCollection& oc); //???
    /** Cache pointer to geometry _for current event_
     */
    void setGeometry(const GEMGeometry* geom); //2 geom?

     private:
    GEMCSCSegmentAlgorithm* algo;
    const GEMGeometry* geom_;

};

#endif
