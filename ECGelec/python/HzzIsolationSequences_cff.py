import FWCore.ParameterSet.Config as cms

#from Configuration.Geometry.GeometryIdeal_cff import *

## original "TDR" HZZ iso
from EGamma.ECGelec.TdrHzzIsolationProducer_cff import *
import EGamma.ECGelec.TdrHzzIsolationProducer_cff 

## Egamma
from EGamma.ECGelec.Eg4HzzIsolationSequences_cff import *
import EGamma.ECGelec.Eg4HzzIsolationSequences_cff

#HzzIsolationSequence = cms.Sequence(TdrHzzIsolationProducer + Eg4HzzElectronIsolationDepositSequence)
HzzIsolationSequence =  cms.Sequence(Eg4HzzElectronIsolationDepositSequence)







