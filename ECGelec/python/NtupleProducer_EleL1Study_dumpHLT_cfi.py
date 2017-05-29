import FWCore.ParameterSet.Config as cms

produceNtupleL1Study = cms.EDAnalyzer(
    "ElectronL1Study_dumpHLT",

    # What to get/printout
    GetL1  = cms.untracked.bool(False),
    GetL1M = cms.untracked.bool(False),
    GetTP  = cms.untracked.bool(False),
    GetTPmodif = cms.untracked.bool(False),
    GetTPemul  = cms.untracked.bool(False),
    GetHcalTP  = cms.untracked.bool(False),
    GetStripMask = cms.untracked.bool(False),
    GetXtalMask  = cms.untracked.bool(False),
    GetVertices  = cms.untracked.bool(False),
    PrintDebug = cms.untracked.bool(False),
    PrintDebug_HLT = cms.untracked.bool(False),
    
    # RecHits tag (get spikes)
    EcalRecHitCollectionEB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    
    # ECAL TP tags
    TPCollectionNormal = cms.InputTag("ecalDigis","EcalTriggerPrimitives"),
    TPCollectionModif  = cms.InputTag("zeroedEcalTrigPrimDigis"),
    TPEmulatorCollection  = cms.InputTag("simEcalTriggerPrimitiveDigis"),
    
    # HCAL TP tag
    #TP_HCAL = cms.InputTag("hcalDigis","HcalTriggerPrimitives"),
    TP_HCAL = cms.InputTag("hcalDigis"),
    hcalTowers = cms.InputTag("towerMaker"),
    
    # GT digis (to get the pre/post firing)
    GTRecordCollection = cms.string('gtDigis'),
    #GTRecordCollectionM = cms.string('simGtDigis'),
    
    EleTag = cms.InputTag("gsfElectrons"),
    VerticesTag = cms.InputTag("offlinePrimaryVerticesWithBS"),
    
    # Trigger tags/paths
    HLTTag      = cms.InputTag("TriggerResults","","HLT"),
    HLT_paths   = cms.vstring("HLT_Ele"),
    HLTElePaths = cms.vstring("HLT_Photon15_L1R"),
    
    # Data to analyze 
    AOD = cms.untracked.bool(False),
    DoFillEle = cms.untracked.bool(False),
    DoFillTrigger = cms.untracked.bool(False),
    DoFillSC = cms.untracked.bool(False),
    DoFillSpikes = cms.untracked.bool(False),
    
    # Configuration
    dcsTag = cms.untracked.InputTag("scalersRawToDigi"),
    type = cms.string("DATA"),
    
    # Compute H/E
    hOverEPtMin = cms.double(0.), 
    hOverEConeSize = cms.double(0.15),
    
    )

