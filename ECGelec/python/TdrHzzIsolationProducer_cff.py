import FWCore.ParameterSet.Config as cms

TdrHzzIsolationProducer = cms.EDProducer("TdrHzzIsolationProducer",
    isolationConeVeto = cms.double(0.015),
    TracksLabel = cms.InputTag("generalTracks"),
    ElectronsVetoLabel = cms.InputTag("gedGsfElectrons"),
    isoVarCut = cms.double(500.0),
    MuonsLabel = cms.InputTag("muons"),
    radiusConeExtHad = cms.double(0.2),
    eTMinHad = cms.double(0.5),
    radiusConeIntHad = cms.double(0.0),
    hcalrhits = cms.string('hbhereco'),
    ElectronsLabel = cms.InputTag("gedGsfElectrons"),
    isolationCone = cms.double(0.25)
)


