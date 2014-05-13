import FWCore.ParameterSet.Config as cms

process = cms.Process("TestGEMCSCSegment")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration.Geometry.GeometryExtended2019Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2019_cff') #?
process.load('Configuration.EventContent.EventContent_cff') #?
#process.load("Geometry.GEMGeometry.me0Geometry_cfi")
process.load("Geometry.GEMGeometry.gemGeometry_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/lustre/cms/store/user/radogna/GEMCSCSegment/out_rec_gemcsc.test500.root'
    )
)

process.gemcscs = cms.EDAnalyzer('TestGEMCSCSegmentAnalyzer',
                              RootFileName = cms.untracked.string("TestGEMCSCSegmentHistograms500.root"),

)

process.p = cms.Path(process.gemcscs)
