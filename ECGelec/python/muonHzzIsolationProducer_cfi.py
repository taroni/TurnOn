import FWCore.ParameterSet.Config as cms

MuonHzzIsolationProducer = cms.EDProducer("MuonHzzIsolationProducer",

    # primary vertex
    PVLabel = cms.InputTag("offlinePrimaryVertices"),

    # deposits
    ECALIsoDepositLabel = cms.InputTag("muIsoDepositCalByAssociatorTowersNew","ecal"),
    HCALIsoDepositLabel = cms.InputTag("muIsoDepositCalByAssociatorTowersNew","hcal"),
    HOCALIsoDepositLabel = cms.InputTag("muIsoDepositCalByAssociatorTowersNew","ho"),
    TrackerIsoDepositLabel = cms.InputTag("muIsoDepositTkNew"),

    # objects
    MuonsLabel = cms.InputTag("muons"),
    TracksLabel = cms.InputTag("generalTracks"),
    ElectronsLabel = cms.InputTag("gedGsfElectrons"),

    # algo parameters
    isolationCone = cms.double(0.3),
    isolationConeVeto = cms.double(0.015),
    trkIsoWeight = cms.double(1.),
    ecalWeight = cms.double(1.),
    hcalWeight = cms.double(1.),

    # isolation cut
    isolationcut = cms.double(60.0)
)


