import FWCore.ParameterSet.Config as cms

gemcscSegments = cms.EDProducer("GEMCSCSegmentProducer",
 #define input
 inputObjectsGEM = cms.InputTag("gemRecHits","","NEWRECO"),
 inputObjectsCSC = cms.InputTag("cscSegments"),#2 input gem rec hits e csc segment
 algo_name = cms.string("GEMCSCSegAlgoRR"),                             
 algo_psets = cms.PSet(
     GEMCSCDebug = cms.untracked.bool(True),
     minHitsPerSegment = cms.uint32(2),
     preClustering = cms.bool(True),
     dXclusBoxMax = cms.double(1.),
     dYclusBoxMax = cms.double(5.),
     preClusteringUseChaining = cms.bool(True),
     dPhiChainBoxMax = cms.double(1.),
     dThetaChainBoxMax = cms.double(1.),
     dRChainBoxMax = cms.double(0.5),
     #dThEtaChainBoxMax = cms.double(1.),
     maxRecHitsInCluster = cms.int32(6)
     )
)
