import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'MC_42_V12::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

# this inputs the input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:PythiaMonoJet_ENERGYIN.root'
    ),
	duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

process.demo = cms.EDAnalyzer('HepMCAnalyzer',
    fileName = cms.string("PythiaMonoJet_ENERGYIN.dat"),
	nparts = cms.int32(PNUM)
)

process.p = cms.Path(process.demo)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#SOURCE process.load('#SRC')
