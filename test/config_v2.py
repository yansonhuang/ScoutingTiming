import FWCore.ParameterSet.Config as cms

process = cms.Process("ScoutingTiming")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = '150X_dataRun3_HLT_v1'
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

with open('validation_6files.txt') as f:
    files = ['root://cms-xrd-global.cern.ch/' + l.strip() for l in f.readlines() if l.strip()]


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        files
    )
)

process.testScoutingTiming = cms.EDAnalyzer("ScoutingTimingAnalyzer",
        pfJetsTag = cms.InputTag( "hltScoutingPFPacker","","HLT"),
          ebRecHitsTag = cms.InputTag("hltScoutingRecHitPacker", "EB", "HLT"),
          ebCleanedRecHitsTag = cms.InputTag("hltScoutingRecHitPacker", "EBCleaned", "HLT"),
            triggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
            triggerResultsRerunTag = cms.InputTag("TriggerResults","","HLT") 
          )

# process.out = cms.OutputModule("PoolOutputModule",
#     fileName = cms.untracked.string("test.root"),
#     # fastCloning = cms.untracked.bool(False)
#     )
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("test_6files_v2.root")
        )
process.path = cms.Path(process.testScoutingTiming)
