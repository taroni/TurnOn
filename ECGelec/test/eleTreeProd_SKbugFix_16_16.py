import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
from Configuration.StandardSequences.Eras import eras
import os
import sys
import commands

options = VarParsing.VarParsing ('analysis')



from Configuration.StandardSequences.Eras import eras

process = cms.Process('RAW2DIGI',eras.Run2_2016)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

options.parseArguments()

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring( 
'/store/data/Run2016H/DoubleEG/RAW-RECO/ZElectron-PromptReco-v2/000/281/207/00000/401A3BC8-6482-E611-B39C-FA163E51B535.root',
'/store/data/Run2016H/DoubleEG/RAW-RECO/ZElectron-PromptReco-v2/000/281/209/00000/34CCF961-6682-E611-8827-02163E01414B.root',
'/store/data/Run2016H/DoubleEG/RAW-RECO/ZElectron-PromptReco-v2/000/281/222/00000/2012CF8A-6882-E611-8DCD-02163E01220E.root'
'/store/data/Run2016H/DoubleEG/RAW-RECO/ZElectron-PromptReco-v2/000/281/231/00000/9219CC0A-6882-E611-84AE-FA163EA38458.root'

       ),
)


import FWCore.PythonUtilities.LumiList as LumiList
process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt').getVLuminosityBlockRange()


process.options = cms.untracked.PSet(

)



# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('l1Ntuple nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag.toGet = cms.VPSet( #we are applying the masking.
    cms.PSet(record = cms.string("EcalTPGFineGrainStripEERcd"),             
             tag = cms.string("EcalTPGFineGrainStrip_16"),
             connect =cms.string('frontier://FrontierPrep/CMS_CONDITIONS')
             ),
    cms.PSet(record = cms.string("EcalTPGSpikeRcd"),
             tag = cms.string("EcalTPGSpike_16"),
             connect =cms.string('frontier://FrontierPrep/CMS_CONDITIONS')
             )
    )

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.endjob_step)

# customisation of the process.

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseReEmul
from L1Trigger.Configuration.customiseReEmul import L1TReEmulFromRAW

#call to customisation function L1TReEmulFromRAW imported from L1Trigger.Configuration.customiseReEmul
process = L1TReEmulFromRAW(process)

# Automatic addition of the customisation function from L1Trigger.L1TNtuples.customiseL1Ntuple
from L1Trigger.L1TNtuples.customiseL1Ntuple import L1NtupleAODRAWEMU 

#call to customisation function L1NtupleAODRAWEMU imported from L1Trigger.L1TNtuples.customiseL1Ntuple
process = L1NtupleAODRAWEMU(process)

# Automatic addition of the customisation function from L1Trigger.Configuration.customiseSettings
from L1Trigger.Configuration.customiseSettings import L1TSettingsToCaloStage2Params_v3_3


#call to customisation function L1TSettingsToCaloStage2Params_2017_v1_0_inconsistent imported from L1Trigger.Configuration.customiseSettings
process = L1TSettingsToCaloStage2Params_v3_3(process)

from EventFilter.L1TRawToDigi.caloStage2Digis_cfi import caloStage2Digis
process.rawCaloStage2Digis = caloStage2Digis.clone()
process.rawCaloStage2Digis.InputLabel = cms.InputTag('rawDataCollector')


# ---------------------------------------------------------------------
# JETS
# ---------------------------------------------------------------------
# JPT
process.load('RecoJets.Configuration.RecoPFJets_cff')


# ---------------------------------------------------------------------
# Fast Jet Rho Correction
# ---------------------------------------------------------------------
process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)


# ---------------------------------------------------------------------
# Vertexing DA
# ---------------------------------------------------------------------
from RecoVertex.Configuration.RecoVertex_cff import *



# ---------------------------------------------------------------------
# Set up electron ID (VID framework)
# ---------------------------------------------------------------------
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.AOD
switchOnVIDElectronIdProducer(process, dataFormat)
process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']
#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

# ---------------------------------------------------------------------
# Unpack Ecal Digis
# ---------------------------------------------------------------------
process.load("EventFilter.EcalRawToDigi.EcalUnpackerMapping_cfi");
process.load("EventFilter.EcalRawToDigi.EcalUnpackerData_cfi");

# ---------------------------------------------------------------------
# Simulate Ecal Trigger Primitives
# ---------------------------------------------------------------------

process.simEcalTriggerPrimitiveDigis.Label = 'ecalDigis'
process.simEcalTriggerPrimitiveDigis.InstanceEB =  'ebDigis'
process.simEcalTriggerPrimitiveDigis.InstanceEE =  'eeDigis'
process.simEcalTriggerPrimitiveDigis.BarrelOnly =False 

process.load('L1Trigger.Configuration.L1Trigger_EventContent_cff')

# L1 extra for the online candidates
process.l1extraParticlesOnline = cms.EDProducer("L1ExtraParticlesProd",
                                                muonSource = cms.InputTag("gtDigis"),
                                                etTotalSource = cms.InputTag("caloStage1LegacyFormatDigis"),
                                                nonIsolatedEmSource = cms.InputTag("caloStage1LegacyFormatDigis","nonIsoEm"),
                                                etMissSource = cms.InputTag("caloStage1LegacyFormatDigis"),
                                                htMissSource = cms.InputTag("caloStage1LegacyFormatDigis"),
                                                produceMuonParticles = cms.bool(False),
                                                forwardJetSource = cms.InputTag("caloStage1LegacyFormatDigis","forJets"),
                                                centralJetSource = cms.InputTag("caloStage1LegacyFormatDigis","cenJets"),
                                                produceCaloParticles = cms.bool(True),
                                                tauJetSource = cms.InputTag("caloStage1LegacyFormatDigis","tauJets"),
                                                isoTauJetSource = cms.InputTag("caloStage1LegacyFormatDigis","isoTauJets"),
                                                isolatedEmSource = cms.InputTag("caloStage1LegacyFormatDigis","isoEm"),
                                                etHadSource = cms.InputTag("caloStage1LegacyFormatDigis"),
                                                hfRingEtSumsSource = cms.InputTag("caloStage1LegacyFormatDigis"),
                                                hfRingBitCountsSource = cms.InputTag("caloStage1LegacyFormatDigis"),
                                                centralBxOnly = cms.bool(True),
                                                ignoreHtMiss = cms.bool(False)
                                                )
# ---------------------------------------------------------------------
# Produce Ntuple Module
# ---------------------------------------------------------------------

HLT_name = 'HLT'


process.load("EGamma.ECGelec.NtupleProducer_custom_cfi")
from EGamma.ECGelec.NtupleProducer_custom_cfi import *
process.produceNtuple = produceNtupleCustom.clone()


process.produceNtuple.NadL1M = cms.untracked.bool(True)
process.produceNtuple.NadTP = cms.untracked.bool(True)
process.produceNtuple.NadTPemul = cms.untracked.bool(True) # Need to put True when running Emulator !!
process.produceNtuple.NadTPmodif = cms.untracked.bool(False)
process.produceNtuple.PrintDebug = cms.untracked.bool(False)

process.produceNtuple.type = 'DATA'
process.produceNtuple.AOD = cms.untracked.bool(False)
process.produceNtuple.FillSC = cms.untracked.bool(True)
process.produceNtuple.functionName = cms.string("EcalClusterEnergyUncertainty")
# Trigger Stuff
process.produceNtuple.HLTTag          = 'TriggerResults::' + HLT_name
process.produceNtuple.TriggerEventTag = 'hltTriggerSummaryAOD::' + HLT_name
process.produceNtuple.HLTElePaths = cms.vstring(
'HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v1',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2',
'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_v2'
        )
process.produceNtuple.HLTMuonPaths    = cms.vstring('HLT_Mu9')
process.produceNtuple.HLTFilters = cms.VInputTag('hltL1NonIsoHLTNonIsoSingleElectronEt17TighterEleIdIsolTrackIsolFilter::'+HLT_name,
                                                      'hltL1NonIsoHLTNonIsoDoubleElectronEt17PixelMatchFilter::'+HLT_name,
                                                      #'hltL1NonIsoHLTNonIsoSingleElectronEt17TightCaloEleIdEle8HEPixelMatchFilter::'+HLT_name,
                                                      'hltL1NonIsoHLTNonIsoSingleElectronEt17TighterEleIdIsolPixelMatchFilter::'+HLT_name,
                                                      'hltL1NonIsoHLTNonIsoSingleElectronEt17TightCaloEleIdEle8HEDoublePixelMatchFilter::'+HLT_name,
                                                      # Muon Trigger
                                                      'hltSingleMu9L3Filtered9')                                                   

process.produceNtuple.eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto")
process.produceNtuple.eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose")
process.produceNtuple.eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium")
process.produceNtuple.eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight")

#process.reco = cms.Sequence(process.reconstruction)
process.L1TSeq = cms.Sequence(   process.RawToDigi #+ process.reconstruction        
#                                   + process.L1TReEmulateFromRAW
)

# ---------------------------------------------------------------------
# Sequence PATH
# ---------------------------------------------------------------------
process.p = cms.Path (
    process.L1TSeq
    + process.kt6PFJets
    + process.egmGsfElectronIDSequence
    + process.ecalEBunpacker
    + process.simEcalTriggerPrimitiveDigis
    #+ process.simRctDigis
    #+ process.simGctDigis
    #+ process.simGtDigis
    + process.caloStage1FinalDigis
    + process.caloStage1LegacyFormatDigis
    #+ process.l1extraParticles
    + process.l1extraParticlesOnline
    + process.rawCaloStage2Digis
    + process.produceNtuple

    )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.schedule.append(process.p)



