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
process.GlobalTag.globaltag = '142X_mcRun3_2025_realistic_v7'
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 'file:/afs/cern.ch/work/m/mcitron/scouting/CMSSW_15_0_12/src/outputScoutingPF0.root',
        'file:/home/users/mcitron/scouting/CMSSW_15_0_12/src/outputScoutingPF0.root'
        # 'file:/afs/cern.ch/work/m/mcitron/scouting/CMSSW_15_0_12/src/outputScoutingPF1.root'
        # 'file:/vols/build/cms/mc3909/timingCMSSW/CMSSW_8_0_27/src/ecalTiming/T1qqqqLL10000/T1qqqqLL1000_AOD_all_10000.root',
        # 'file:/vols/build/cms/mc3909/timingCMSSW/CMSSW_8_0_27/src/ecalTiming/T1qqqqLL1000/T1qqqqLL1000_AOD_500To999.root',
        # 'file:/vols/build/cms/mc3909/timingCMSSW/CMSSW_8_0_27/src/ecalTiming/T1qqqqLL1000/T1qqqqLL1000_AOD_200.root',
        # 'file:/vols/build/cms/mc3909/timingCMSSW/CMSSW_8_0_27/src/ecalTiming/T1qqqqLL1000/T1qqqqLL1000_AOD_200To500.root'
        # 'file:/vols/build/cms/mc3909/timingCMSSW/CMSSW_8_0_27/src/ecalTiming/T1qqqqLL0p001/T1qqqqLL1000_AOD_all_0p001.root'
        # 'file:/vols/build/cms/mc3909/timingCMSSW/CMSSW_8_0_27/src/ecalTiming/tempT1qqqqLL1000_AOD_all_10000.root'
        # 'file:/vols/build/cms/mc3909/timingCMSSW/CMSSW_8_0_25/src/dataFiles/2017_04_13_16_39_PRnewco_80X_dataRun2_2016LegacyRepro_v3-v1.root'
        #'file:/vols/build/cms/mc3909/qcdAODSIM.root'
    )
)

process.testScoutingTiming = cms.EDAnalyzer("ScoutingTimingAnalyzer",
        pfJetsTag = cms.InputTag( "hltScoutingPFPacker","","HLTX"),
          ebRecHitsTag = cms.InputTag("hltScoutingRecHitPacker", "EB", "HLTX"),
            triggerResultsTag = cms.InputTag("TriggerResults","","SIM"),
            triggerResultsRerunTag = cms.InputTag("TriggerResults","","HLTX") 
          )

# process.out = cms.OutputModule("PoolOutputModule",
#     fileName = cms.untracked.string("test.root"),
#     # fastCloning = cms.untracked.bool(False)
#     )
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("test1.root")
        )
process.path = cms.Path(process.testScoutingTiming)
