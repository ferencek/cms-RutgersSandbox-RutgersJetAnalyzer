###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')

options.register('runOnData',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)   
options.register('globalTag',
    'START52_V9D::All',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Global tag to be used"
)
options.register('outfilename',
    "outfile.root",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "output file name"
)
options.register('reportEvery',
    10,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=10)"
)
options.register('wantSummary',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Print out trigger and timing summary"
)
options.register('usePFchs',
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,   
    "Use PFchs"
)
## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 100)

options.parseArguments()
    
print "Running on data: %s"%('True' if options.runOnData else 'False')
print "Using PFchs: %s"%('True' if options.usePFchs else 'False')
    
## Jet energy corrections
inputJetCorrLabelAK5 = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
inputJetCorrLabelAK7 = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
if not options.usePFchs:
    inputJetCorrLabelAK5 = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute'])
    inputJetCorrLabelAK7 = ('AK7PF', ['L1FastJet', 'L2Relative', 'L3Absolute'])

import FWCore.ParameterSet.Config as cms
    
process = cms.Process("PAT")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.MessageLogger.cerr.default.limit = 10
#################################################################

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
## Make sure a correct global tag is used (please refer to https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Valid_Global_Tags_by_Release)
process.GlobalTag.globaltag = options.globalTag

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.wantSummary) )

## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/Summer12/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/PU_S7_START52_V9-v1/0000/FC07A5F9-8495-E111-9E8E-00304867903E.root'
    )
)

## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('patTuple_PF2PAT.root'),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    # save PAT Layer 1 output; you need a '*' to
    # unpack the list of commands 'patEventContent'
    outputCommands = cms.untracked.vstring('drop *', *patEventContent )
)

## Define Endpath
process.outpath = cms.EndPath(process.out)

## Configure PAT to use PF2PAT instead of AOD sources
## this function will modify the PAT sequences.
from PhysicsTools.PatAlgos.tools.pfTools import *

postfix = "PFlow"
jetAlgo="AK7"
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=not options.runOnData, postfix=postfix,
          jetCorrections=inputJetCorrLabelAK7, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'))

## Top projections in PF2PAT
getattr(process,"pfNoPileUp"+postfix).enable = options.usePFchs
getattr(process,"pfNoMuon"+postfix).enable = True
getattr(process,"pfNoElectron"+postfix).enable = True
getattr(process,"pfNoTau"+postfix).enable = False
getattr(process,"pfNoJet"+postfix).enable = True

## Define AK6 jets (GEN and RECO)
process.ak6GenJetsNoNu = process.ak5GenJetsNoNu.clone( rParam = 0.6 )
process.ak6PFlow = process.pfJetsPFlow.clone( rParam = 0.6 )

## Define AK8 jets (GEN and RECO)
process.ak8GenJetsNoNu = process.ak5GenJetsNoNu.clone( rParam = 0.8 )
process.ak8PFlow = process.pfJetsPFlow.clone( rParam = 0.8 )

## Define AK10 jets (GEN and RECO)
process.ak10GenJetsNoNu = process.ak5GenJetsNoNu.clone( rParam = 1.0 )
process.ak10PFlow = process.pfJetsPFlow.clone( rParam = 1.0 )

##AK6 jets groomed

from RecoJets.JetProducers.ak5PFJetsTrimmed_cfi import ak5PFJetsTrimmed
process.ak6TrimmedPFlow = ak5PFJetsTrimmed.clone(
    src = process.pfJetsPFlow.src,
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(0.6)
    )

from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsFiltered
process.ak6FilteredPFlow = ak5PFJetsFiltered.clone(
    src = process.pfJetsPFlow.src,
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(0.6)
    )

from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.ak6PrunedPFlow = ak5PFJetsPruned.clone(
    src = process.pfJetsPFlow.src,
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(0.6)
    )

##AK8 jets groomed

process.ak8TrimmedPFlow = ak5PFJetsTrimmed.clone(
    src = process.pfJetsPFlow.src,
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(0.8)
    )

process.ak8FilteredPFlow = ak5PFJetsFiltered.clone(
    src = process.pfJetsPFlow.src,
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(0.8)
    )

process.ak8PrunedPFlow = ak5PFJetsPruned.clone(
    src = process.pfJetsPFlow.src,
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(0.8)
    )

##AK10 jets groomed
    
process.ak10TrimmedPFlow = ak5PFJetsTrimmed.clone(
    src = process.pfJetsPFlow.src,
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(1.0)
    )
    
process.ak10FilteredPFlow = ak5PFJetsFiltered.clone(
    src = process.pfJetsPFlow.src,
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(1.0)
    )

process.ak10PrunedPFlow = ak5PFJetsPruned.clone(
    src = process.pfJetsPFlow.src,
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(1.0)
    )                          

from PhysicsTools.PatAlgos.tools.jetTools import *
## Add AK6 jets to PAT
addJetCollection(
    process,
    cms.InputTag('ak6PFlow'),
    'AK6', 'PF',
    doJTA=True,
    doBTagging=True,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID = False,
    genJetCollection = cms.InputTag("ak6GenJetsNoNu")
)

## Add AK6 Trimmed jets to PAT
addJetCollection(
    process,
    cms.InputTag('ak6TrimmedPFlow'),
    'AK6Trimmed','PF',
    doJTA=True,
    doBTagging=True,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID = False,
    genJetCollection = cms.InputTag("ak6GenJetsNoNu")
)

## Add AK6 Filtered jets to PAT
addJetCollection(
    process,
    cms.InputTag('ak6FilteredPFlow'),
    'AK6Filtered','PF',
    doJTA=True,
    doBTagging=True,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID = False,
    genJetCollection = cms.InputTag("ak6GenJetsNoNu")
)

## Add AK6 Pruned jets to PAT
addJetCollection(
    process,
    cms.InputTag('ak6PrunedPFlow'),
    'AK6Pruned','PF',
    doJTA=True,
    doBTagging=True,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID = False,
   genJetCollection = cms.InputTag("ak6GenJetsNoNu")
)

## Define a sequence for RECO jets and append it to the PF2PAT sequence
process.jetSeq = cms.Sequence(
    process.ak6GenJetsNoNu+
    process.ak6PFlow+
    process.ak6TrimmedPFlow+
    process.ak6FilteredPFlow+
    process.ak6PrunedPFlow
)

## Explicitly specify JEC levels for the default PAT jets (otherwise, the code will crash since it looks for L5Flavor and L7Parton which are not available in the latest JECs)
getattr(process,"patJetCorrFactors").levels = cms.vstring('L1Offset','L2Relative','L3Absolute')

## Produce a collection of good primary vertices
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone(
        minNdof = cms.double(4.0), # this is >= 4
        maxZ = cms.double(24.0),
        maxRho = cms.double(2.0)
    ),
    src = cms.InputTag('offlinePrimaryVertices')
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outfilename)
)

process.rutgersJetAnalyzerak6 = cms.EDAnalyzer('RutgersJetAnalyzer',
    GenJetsTag = cms.InputTag('ak12GenJets'),
    GenParticleTag = cms.InputTag('genParticles'),
    InputsTag = cms.InputTag('genParticlesForJets'),
    JetsTag = cms.InputTag('selectedPatJetsAK6PF'),
    pvtag = cms.InputTag('goodOfflinePrimaryVertices'),
    InputPtMin = cms.double(0.0),
    JetPtMin = cms.double(0.0)
)

process.rutgersJetAnalyzerak6Trimmed = cms.EDAnalyzer('RutgersJetAnalyzer',
    GenJetsTag = cms.InputTag('ak12GenJets'),
    GenParticleTag = cms.InputTag('genParticles'),
    InputsTag = cms.InputTag('genParticlesForJets'),
    JetsTag = cms.InputTag('selectedPatJetsAK6TrimmedPF'),
    pvtag = cms.InputTag('goodOfflinePrimaryVertices'),
    InputPtMin = cms.double(0.0),
    JetPtMin = cms.double(0.0)
)

process.rutgersJetAnalyzerak6Filtered = cms.EDAnalyzer('RutgersJetAnalyzer',
    GenJetsTag = cms.InputTag('ak12GenJets'),
    GenParticleTag = cms.InputTag('genParticles'),
    InputsTag = cms.InputTag('genParticlesForJets'),
    JetsTag = cms.InputTag('selectedPatJetsAK6FilteredPF'),
    pvtag = cms.InputTag('goodOfflinePrimaryVertices'),
    InputPtMin = cms.double(0.0),
    JetPtMin = cms.double(0.0)
)

process.rutgersJetAnalyzerak6Pruned = cms.EDAnalyzer('RutgersJetAnalyzer',
    GenJetsTag = cms.InputTag('ak12GenJets'),
    GenParticleTag = cms.InputTag('genParticles'),
    InputsTag = cms.InputTag('genParticlesForJets'),
    JetsTag = cms.InputTag('selectedPatJetsAK6PrunedPF'),
    pvtag = cms.InputTag('goodOfflinePrimaryVertices'),
    InputPtMin = cms.double(0.0),
    JetPtMin = cms.double(0.0)
)

## Path definition
process.p = cms.Path(
    process.goodOfflinePrimaryVertices*
    getattr(process,"patPF2PATSequence"+postfix)*
    process.jetSeq*
    process.patDefaultSequence*
    process.rutgersJetAnalyzerak6*
    process.rutgersJetAnalyzerak6Trimmed*
    process.rutgersJetAnalyzerak6Filtered*
    process.rutgersJetAnalyzerak6Pruned
)

### Add PF2PAT output to the created file
#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
##process.load("CommonTools.ParticleFlow.PF2PAT_EventContent_cff")
##process.out.outputCommands =  cms.untracked.vstring('drop *')
#process.out.outputCommands = cms.untracked.vstring('drop *',
                                                   #'keep recoPFCandidates_particleFlow_*_*',
                                                   #*patEventContentNoCleaning )

## Delete predefined Endpath (needed for running with CRAB)
del process.out
del process.outpath

## Schedule definition
process.schedule = cms.Schedule(process.p)
#process.schedule.append(process.outpath)


