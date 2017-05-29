
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'turnOn_16_16'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'eleTreeProd_16_16.py'

config.Data.inputDataset =  '/DoubleEG/Run2016G-ZElectron-PromptReco-v1/RAW-RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 25
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-282092_13TeV_PromptReco_Collisions16_JSON.txt'
config.Data.runRange = '279993,279931,279071'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'turnOn_16_16'

config.Site.storageSite = 'T2_CH_CERN'
