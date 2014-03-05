import FWCore.ParameterSet.Config as cms

gemcscSegments = cms.EDProducer("GEMCSCSegmentProducer",
 #define input
 inputObjects = cms.InputTag("gemRecHits"),
 algo_name = cms.string("GEMCSCSegAlgoRR"),                             
 algo_pset = cms.PSet(
     ME0Debug = cms.untracked.bool(True),
     minHitsPerSegment = cms.uint32(2),
     preClustering = cms.bool(True),
     dXclusBoxMax = cms.double(1.),
     dYclusBoxMax = cms.double(5.),
     preClusteringUseChaining = cms.bool(True),
     dPhiChainBoxMax = cms.double(1.),
     dEtaChainBoxMax = cms.double(1.),
     maxRecHitsInCluster = cms.int32(6)
     )
)
