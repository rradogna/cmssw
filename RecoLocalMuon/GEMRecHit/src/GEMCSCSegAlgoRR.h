#ifndef GEMRecHit_GEMCSCSegAlgoRR_h
#define GEMRecHit_GEMCSCSegAlgoRR_h

/**
 * \class GEMCSCSegAlgoRR
 *
 * This algorithm is very basic no attemp to deal with ambiguities , noise etc.
 * The GEMCSC track segments is built out of the rechit's in a the 2 GEM Layer
 and the CSC segment denoted
 * as the GEMCSC Ensabmle .<BR>
 *
 *  \authors Raffaella Radogna 
 *
 */

#include <RecoLocalMuon/GEMRecHit/src/GEMCSCSegmentAlgorithm.h>
#include <DataFormats/GEMRecHit/interface/GEMRecHit.h>
#include <DataFormats/CSCRecHit/interface/CSCSegment.h>

#include <deque>
#include <vector>

class GEMCSCSegAlgoRR : public GEMCSCSegmentAlgorithm {


public:

  /// Typedefs

  typedef std::vector<const GEMRecHit*> EnsambleGEMHitContainer;
  typedef std::vector<const CSCSegment*> EnsambleCSCSegContainer;
  typedef std::vector<const RecHit2DLocalPos*> EnsambleHitContainer;
  typedef std::vector<EnsambleHitContainer> ProtoSegments;
  typedef std::deque<bool> BoolContainer;

  /// Constructor
  explicit GEMCSCSegAlgoRR(const edm::ParameterSet& ps);
  /// Destructor
  virtual ~GEMCSCSegAlgoRR();

  /**
   * Build segments for all desired groups of hits
   */
  std::vector<GEMCSCSegment> run(GEMCSCEnsamble ensamble, const EnsambleCSCSegContainer& cscsegments, const EnsambleGEMHitContainer& rechits); 

private:
  /// Utility functions 

  //  Build groups of rechits that are separated in x and y to save time on the segment finding
  //ProtoSegments clusterHits(const EnsambleHitContainer & rechits);

  // Build groups of rechits that are separated in strip numbers and Z to save time on the segment finding
  //ProtoSegments chainHitsToSegm(const EnsambleCSCSegContainer cscsegments, const EnsambleGEMHitContainer & rechits);
  EnsambleHitContainer chainHitsToSegm(const CSCSegment* cscsegment, const EnsambleGEMHitContainer & rechits);
  //bool isGoodToMerge(EnsambleHitContainer & newChain, EnsambleHitContainer & oldChain);

  // Build track segments in this chamber (this is where the actual segment-building algorithm hides.)
  std::vector<GEMCSCSegment> buildSegments(const CSCSegment* cscsegment, const EnsambleHitContainer& rechits);

  void doSlopesAndChi2();
  void fitSlopes();
  void fillChiSquared();
  void fillLocalDirection();
  CLHEP::HepMatrix derivativeMatrix(void);
  AlgebraicSymMatrix weightMatrix(void);
  AlgebraicSymMatrix calculateError(void);
  void flipErrors(AlgebraicSymMatrix& protoErrors);

  // Member variables
 private:
  const std::string myName; 

  // input from .cfi file
 private:
  bool    debug;
  unsigned int     minHitsPerSegment;
  bool    preClustering;
  double  dXclusBoxMax;
  double  dYclusBoxMax;
  bool    preClustering_useChaining;
  double  dPhiChainBoxMax;
  double  dThEtaChainBoxMax;
  int     maxRecHitsInCluster;
  
 private:
  EnsambleGEMHitContainer gemrhs;
  EnsambleCSCSegContainer csc_segment;
  EnsambleHitContainer proto_segment;
  GEMCSCEnsamble theEnsamble;
  LocalPoint protoIntercept;
  float protoSlope_u;
  float protoSlope_v;
  double protoChi2;
  double protoNDF;
  LocalVector protoDirection;
  double protoChiUCorrection;

  std::vector<double> e_Cxx;
};

#endif
