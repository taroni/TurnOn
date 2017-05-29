import FWCore.ParameterSet.Config as cms

from RecoMuon.MuonIsolationProducers.muIsoDepositJets_cfi  import *
import RecoMuon.MuonIsolationProducers.muIsoDepositJets_cfi 
from RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi  import *
from RecoMuon.MuonIsolationProducers.muIsoDepositTk_cfi import *
import RecoMuon.MuonIsolationProducers.muIsoDepositTk_cfi
muIsoDepositTkNew=RecoMuon.MuonIsolationProducers.muIsoDepositTk_cfi.muIsoDepositTk.clone()
muIsoDepositTkNew.IOPSet.inputMuonCollection = cms.InputTag("muons")

#import RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi 
muIsoDepositCalByAssociatorTowersNew=RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi.muIsoDepositCalByAssociatorTowers.clone()
muIsoDepositCalByAssociatorTowersNew.IOPSet.inputMuonCollection = cms.InputTag("muons")

muIsoDepositJetsNew=RecoMuon.MuonIsolationProducers.muIsoDepositJets_cfi.muIsoDepositJets.clone()
muIsoDepositJetsNew.IOPSet.inputMuonCollection = cms.InputTag("muons")

muIsoDeposits_muonsNew = cms.Sequence(muIsoDepositTkNew+muIsoDepositCalByAssociatorTowersNew+muIsoDepositJetsNew)
muIsolation_muonsNew = cms.Sequence(muIsoDeposits_muonsNew)
muIsolationNew = cms.Sequence(muIsolation_muonsNew)


from EGamma.ECGelec.muonHzzIsolationProducer_cfi import *

MuonHZZIsolationSequence = cms.Sequence( muIsolationNew + MuonHzzIsolationProducer )
