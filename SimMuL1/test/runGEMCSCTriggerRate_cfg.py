import FWCore.ParameterSet.Config as cms
import os

## initialization
process = cms.Process('GEMCSCTRGANA')

## CMSSW RELEASE
cmssw = os.getenv( "CMSSW_VERSION" )

## steering
deltaMatch = 2
minSimTrDR = 0.0
pileup = 100
events = 1000000
#sample = 'dimu'
sample = 'minbias'
globalTag = 'upgrade2019'

## input
from GEMCode.SimMuL1.GEMCSCTriggerSamplesLib import files
suffix = '_gem_dphi_pat2_PU0'
inputDir = files[suffix]
inputDir = ['tempDir/']
theInputFiles = []
import os
for d in range(len(inputDir)):
  my_dir = inputDir[d]
  if not os.path.isdir(my_dir):
    print "ERROR: This is not a valid directory: ", my_dir
    if d==len(inputDir)-1:
      print "ERROR: No input files were selected"
      exit()
    continue
  print "Proceed to next directory"
  ls = os.listdir(my_dir)
  theInputFiles.extend([my_dir[16:] + x for x in ls if x.endswith('root')])
    
print "inputFiles:", theInputFiles

## readout windows
w = 3
if w==3:
    readout_windows = [ [5,7],[5,7],[5,7],[5,7] ]
if w==11:
    readout_windows = [ [1,11],[1,11],[1,11],[1,11] ]
if w==7:
    readout_windows = [ [5,11],[5,11],[5,11],[5,11] ]
if w==61:
    readout_windows = [ [5,10],[1,11],[1,11],[1,11] ]

## output
outputFileName = 'hp_' + sample + "_" + cmssw + "_" + globalTag + "_pu%d"%(pileup) + '_w%d'%(w) + suffix + '.root'
print "outputFile:", outputFileName

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2019Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2019_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.GlobalTag.globaltag = 'auto:upgrade2019'
process.load('L1TriggerConfig.L1ScalesProducers.L1MuTriggerScalesConfig_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.L1Emulator_cff")
process.load("Configuration.StandardSequences.L1Extra_cff")
process.load("SimMuon.CSCDigitizer.muonCSCDigis_cfi")
process.load('L1Trigger.CSCTrackFinder.csctfTrackDigisUngangedME1a_cfi')
process.simCsctfTrackDigis = process.csctfTrackDigisUngangedME1a.clone()
process.simCsctfTrackDigis.DTproducer = cms.untracked.InputTag("simDtTriggerPrimitiveDigis")
process.simCsctfTrackDigis.SectorReceiverInput = cms.untracked.InputTag("simCscTriggerPrimitiveDigis","MPCSORTED")
process.simCsctfTrackDigis.SectorProcessor.isCoreVerbose = cms.bool(True)

process.options = cms.untracked.PSet(
#     Rethrow = cms.untracked.vstring('ProductNotFound'),
#     FailPath = cms.untracked.vstring('ProductNotFound'),
    makeTriggerResults = cms.untracked.bool(False),
    wantSummary = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
#    inputCommands = cms.untracked.vstring(
#      'keep  *_*_*_*',
#      'drop *_simDtTriggerPrimitiveDigis_*_MUTRG'
#    ),
    fileNames = cms.untracked.vstring(
      *theInputFiles
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(events)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(outputFileName)
)

from GEMCode.GEMValidation.simTrackMatching_cfi import SimTrackMatching
#SimTrackMatching.verboseSimHit = 1
#SimTrackMatching.verboseGEMDigi = 1
#SimTrackMatching.verboseCSCDigi = 1
#SimTrackMatching.verboseCSCStub = 1
#SimTrackMatching.simMuOnlyGEM = False
#SimTrackMatching.simMuOnlyCSC = False
#SimTrackMatching.discardEleHitsCSC = False
#SimTrackMatching.discardEleHitsGEM = False

process.load('GEMCode.SimMuL1.GEMCSCTriggerRate_cfi')
process.GEMCSCTriggerRate.minDeltaWire = -1*deltaMatch
process.GEMCSCTriggerRate.maxDeltaWire = deltaMatch
process.GEMCSCTriggerRate.minDeltaStrip = -1*deltaMatch
process.GEMCSCTriggerRate.simTrackGEMMatching = process.SimTrackMatching
process.GEMCSCTriggerRate.minBxALCT = readout_windows[0][0]
process.GEMCSCTriggerRate.maxBxALCT = readout_windows[0][1]
process.GEMCSCTriggerRate.minBxCLCT = readout_windows[1][0]
process.GEMCSCTriggerRate.maxBxCLCT = readout_windows[1][1]
process.GEMCSCTriggerRate.minBxLCT = readout_windows[2][0]
process.GEMCSCTriggerRate.maxBxLCT = readout_windows[2][1]
process.GEMCSCTriggerRate.minBxMPLCT = readout_windows[3][0]
process.GEMCSCTriggerRate.maxBxMPLCT = readout_windows[3][1]
process.GEMCSCTriggerRate.minSimTrDR = minSimTrDR
process.GEMCSCTriggerRate.minSimTrEta = theMinSimTrEta
process.GEMCSCTriggerRate.maxSimTrEta = theMaxSimTrEta
process.GEMCSCTriggerRate.sectorProcessor = process.simCsctfTrackDigis.SectorProcessor
process.GEMCSCTriggerRate.strips = process.simMuonCSCDigis.strips
print process.GEMCSCTriggerRate.strips


## customization 
from SLHCUpgradeSimulations.Configuration.muonCustoms import *
process = customise_csc_Geometry(process)
process = digitizer_timing_pre3_median(process)
process = customise_csc_L1Extra_allsim(process)

## Sequence and schedule
process.ana_seq = cms.Sequence(process.GEMCSCTriggerRate)
process.l1extra_step = cms.Path(process.L1Extra)
process.ana_step = cms.Path(process.ana_seq)

process.schedule = cms.Schedule(
#    process.l1extra_step,
    process.ana_step
)
