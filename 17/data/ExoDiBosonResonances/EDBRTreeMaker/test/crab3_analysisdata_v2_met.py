from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName   = 'MET-17B-v1_new'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
#config.JobType.inputFiles = ['Fall17_17Nov2017B_V32_DATA_L1FastJet_AK4PFchs.txt','Fall17_17Nov2017B_V32_DATA_L2Relative_AK4PFchs.txt','Fall17_17Nov2017B_V32_DATA_L3Absolute_AK4PFchs.txt','Fall17_17Nov2017B_V32_DATA_L2L3Residual_AK4PFchs.txt','Fall17_17Nov2017B_V32_DATA_L1FastJet_AK8PFchs.txt','Fall17_17Nov2017B_V32_DATA_L2Relative_AK8PFchs.txt','Fall17_17Nov2017B_V32_DATA_L3Absolute_AK8PFchs.txt','Fall17_17Nov2017B_V32_DATA_L2L3Residual_AK8PFchs.txt','Fall17_17Nov2017B_V32_DATA_L1FastJet_AK8PFPuppi.txt','Fall17_17Nov2017B_V32_DATA_L2Relative_AK8PFPuppi.txt','Fall17_17Nov2017B_V32_DATA_L3Absolute_AK8PFPuppi.txt','Fall17_17Nov2017B_V32_DATA_L2L3Residual_AK8PFPuppi.txt','L1PrefiringMaps_new.root']
config.JobType.inputFiles = ['Fall17_17Nov2017B_V32_DATA_L1FastJet_AK4PFPuppi.txt','Fall17_17Nov2017B_V32_DATA_L2Relative_AK4PFPuppi.txt','Fall17_17Nov2017B_V32_DATA_L3Absolute_AK4PFPuppi.txt','Fall17_17Nov2017B_V32_DATA_L2L3Residual_AK4PFPuppi.txt','Fall17_17Nov2017B_V32_DATA_L1FastJet_AK8PFchs.txt','Fall17_17Nov2017B_V32_DATA_L2Relative_AK8PFchs.txt','Fall17_17Nov2017B_V32_DATA_L3Absolute_AK8PFchs.txt','Fall17_17Nov2017B_V32_DATA_L2L3Residual_AK8PFchs.txt','Fall17_17Nov2017B_V32_DATA_L1FastJet_AK8PFPuppi.txt','Fall17_17Nov2017B_V32_DATA_L2Relative_AK8PFPuppi.txt','Fall17_17Nov2017B_V32_DATA_L3Absolute_AK8PFPuppi.txt','Fall17_17Nov2017B_V32_DATA_L2L3Residual_AK8PFPuppi.txt','L1PrefiringMaps_new.root']
# Name of the CMSSW configuration file
#config.JobType.psetName    = 'bkg_ana.py'
config.JobType.psetName    = 'analysis.py'
#config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 3000
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 3000

config.section_("Data")
config.Data.inputDataset = '/MET/Run2017B-31Mar2018-v1/MINIAOD'
config.Data.inputDBS = 'global'

config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 200
config.Data.lumiMask = 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
config.Data.runRange = ''#'250843-250932' # '193093-194075'
config.Data.publication = False
name = 'WWW'
steam_dir = 'xulyu'
##config.Data.outLFNDirBase = '/store/group/dpg_trigger/comm_trigger/TriggerStudiesGroup/STEAM/' + steam_dir + '/' + name + '/'
config.Data.outputDatasetTag = 'MET-17B-v1_new'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T3_US_FNALLPC'


##config.Data.inputDataset = '/WJetsToLNu_13TeV-madgraph-pythia8-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
#config.Data.inputDBS = 'global'

##config.Data.inputDBS = 'phys03'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob =10 
#config.Data.totalUnits = 279
#config.Data.publication = False
#
## This string is used to construct the output dataset name
#config.Data.outputDatasetTag = 'WJets100To200_weight'
#
#config.section_("Site")
## Where the output files will be transmitted to
