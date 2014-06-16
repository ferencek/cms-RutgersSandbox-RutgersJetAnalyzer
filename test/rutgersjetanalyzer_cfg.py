###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing
import string

options = VarParsing ('python')

options.register('runOnData', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
## Make sure correct global tags are used (please refer to https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions)
options.register('mcGlobalTag', 'START53_V27',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "MC global tag"
)
options.register('dataGlobalTag', 'FT53_V21A_AN6',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Data global tag"
)
options.register('outFilename', 'outfile',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
)
options.register('useExternalInput', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use external input"
)
options.register('externalInput', '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Path of an external list of input files"
)
options.register('dumpPythonCfg', '',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Name of the rewritten cfg file"
)
options.register('reportEvery', 10,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=10)"
)
options.register('wantSummary', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Print out trigger and timing summary"
)
options.register('usePFchs', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use PFchs"
)
options.register('doJTA', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run jet-track association"
)
options.register('useExplicitJTA', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use explicit jet-track association"
)
options.register('doBTagging', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run b tagging"
)
options.register('jetRadius', 0.8,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Distance parameter R for jet clustering (default is 0.8)"
)
options.register('doBosonMatching', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Do boson matching"
)
options.register('applyBosonIsolation', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Apply boson isolation"
)
options.register('useEventWeight', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use event weight"
)
options.register('useAK5Jets', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use AK5 jets"
)
options.register('runQCDFlavorExtra', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run QCD flavor extra"
)
options.register('runOnWBkg', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run on W background"
)
options.register('runOnZBkg', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run on Z background"
)
options.register('runOnTopBkg', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run on top background"
)
options.register('useSVClustering', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use SV clustering"
)
options.register('useSVMomentum', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use SV momentum"
)
options.register('useRadionCuts', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use cuts optimized for Radion samples"
)
## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 100)

options.parseArguments()

print "Running on data: %s"%('True' if options.runOnData else 'False')
print "Using PFchs: %s"%('True' if options.usePFchs else 'False')

## Global tag
globalTag = options.mcGlobalTag
if options.runOnData:
    globalTag = options.dataGlobalTag

## Jet energy corrections
inputJetCorrLabelAK5 = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
inputJetCorrLabelAK7 = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])

if not options.usePFchs:
    inputJetCorrLabelAK5 = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute'])
    inputJetCorrLabelAK7 = ('AK7PF', ['L1FastJet', 'L2Relative', 'L3Absolute'])

if options.runOnData:
    inputJetCorrLabelAK5[1].append('L2L3Residual')
    inputJetCorrLabelAK7[1].append('L2L3Residual')

## b tagging
bTagInfos = ['impactParameterTagInfos','secondaryVertexTagInfos','inclusiveSecondaryVertexFinderTagInfos']
             #,'inclusiveSecondaryVertexFinderFilteredTagInfos','softMuonTagInfos','secondaryVertexNegativeTagInfos']
bTagDiscriminators = ['jetProbabilityBJetTags','jetBProbabilityBJetTags','combinedSecondaryVertexBJetTags'
                      #,'trackCountingHighPurBJetTags','trackCountingHighEffBJetTags'
                      #,'simpleSecondaryVertexHighPurBJetTags','simpleSecondaryVertexHighEffBJetTags'
                      ,'combinedInclusiveSecondaryVertexBJetTags']
                      #,'simpleInclusiveSecondaryVertexHighEffBJetTags','simpleInclusiveSecondaryVertexHighPurBJetTags'
                      #,'doubleSecondaryVertexHighEffBJetTags']

import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.MessageLogger.cerr.default.limit = 10
#################################################################

## Geometry and Detector Conditions
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = globalTag + '::All'

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.wantSummary) )

#-------------------------------------
## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # /RadionToHH_4b_M-600_TuneZ2star_8TeV-Madgraph_pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM
        'root://cmsxrootd-site.fnal.gov//store/user/ferencek/noreplica/RadionToHH_4b_M-600_TuneZ2star_8TeV-Madgraph_pythia6/Summer12_DR53X-PU_S10_START53_V19-v1_PATTuple_v3/patTuple_PF2PAT_v3_1_1_jHR.root'
        #'file:/cms/ferencek/store/ferencek/RadionToHH_4b_M-600_TuneZ2star_8TeV-Madgraph_pythia6/Summer12_DR53X-PU_S10_START53_V19-v1_PATTuple_v3/patTuple_PF2PAT_v3_1_1_jHR.root'
    )
)
# If using external input files
if options.useExternalInput:
    process.source.fileNames = cms.untracked.vstring( open(options.externalInput,"r").read().splitlines() )

#-------------------------------------
outFilename = string.replace(options.outFilename,'.root','') + '_mc.root'
if options.runOnData :
    outFilename = string.replace(options.outFilename,'.root','') + '_data.root'

## Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string(outFilename)
)
#-------------------------------------
## Standard PAT Configuration File
process.load("PhysicsTools.PatAlgos.patSequences_cff")

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("test.root"),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    # save PAT Layer 1 output; you need a '*' to
    # unpack the list of commands 'patEventContent'
    outputCommands = cms.untracked.vstring('drop *', *patEventContent)
)

#-------------------------------------
## CA jets (Gen and Reco)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.caGenJetsNoNu = ca4GenJets.clone(
    rParam = cms.double(options.jetRadius),
    src = cms.InputTag("genParticlesForJetsNoNu")
)
from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
process.caPFJetsCHS = ca4PFJets.clone(
    rParam = cms.double(options.jetRadius),
    src = cms.InputTag("pfNoElectronPFlow"),
    srcPVs = cms.InputTag("goodOfflinePrimaryVertices"),
    doAreaFastjet = cms.bool(True),
    jetPtMin = cms.double(20.)
)
## CA filtered jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.caGenJetsNoNuFiltered = ca4GenJets.clone(
    rParam = cms.double(options.jetRadius),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    useFiltering = cms.bool(True),
    nFilt = cms.int32(3),
    rFilt = cms.double(0.3),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsFiltered
process.caPFJetsCHSFiltered = ak5PFJetsFiltered.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(options.jetRadius),
    src = process.caPFJetsCHS.src,
    srcPVs = process.caPFJetsCHS.srcPVs,
    doAreaFastjet = process.caPFJetsCHS.doAreaFastjet,
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20.)
)
## CA MassDrop-BDRS filtered jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
## Compared to the above filtered jets, here dynamic filtering radius is used (as in arXiv:0802.2470)
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.caGenJetsNoNuMDBDRSFiltered = ca4GenJets.clone(
    rParam = cms.double(options.jetRadius),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    useMassDropTagger = cms.bool(True),
    muCut = cms.double(0.667),
    yCut = cms.double(0.08),
    useFiltering = cms.bool(True),
    useDynamicFiltering = cms.bool(True),
    nFilt = cms.int32(3),
    rFilt = cms.double(0.3),
    rFiltFactor = cms.double(0.5),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsMassDropFiltered
process.caPFJetsCHSMDBDRSFiltered = ak5PFJetsMassDropFiltered.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(options.jetRadius),
    src = process.caPFJetsCHS.src,
    srcPVs = process.caPFJetsCHS.srcPVs,
    doAreaFastjet = process.caPFJetsCHS.doAreaFastjet,
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20.),
    useDynamicFiltering = cms.bool(True),
    rFiltFactor = cms.double(0.5)
)
## CA Kt-BDRS filtered jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
## Compared to the above filtered jets, here dynamic filtering radius is used (as in arXiv:0802.2470)
## However, here the mass drop is replaced by finding two Kt subjets which then set the size of the filtering radius
from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
process.caGenJetsNoNuKtBDRSFiltered = ca4GenJets.clone(
    rParam = cms.double(options.jetRadius),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    usePruning = cms.bool(True),
    useKtPruning = cms.bool(True),
    zcut = cms.double(0.),
    rcut_factor = cms.double(9999.),
    useFiltering = cms.bool(True),
    useDynamicFiltering = cms.bool(True),
    nFilt = cms.int32(3),
    rFilt = cms.double(0.3),
    rFiltFactor = cms.double(0.5),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsFiltered
process.caPFJetsCHSKtBDRSFiltered = ak5PFJetsFiltered.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(options.jetRadius),
    src = process.caPFJetsCHS.src,
    srcPVs = process.caPFJetsCHS.srcPVs,
    doAreaFastjet = process.caPFJetsCHS.doAreaFastjet,
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20.),
    usePruning = cms.bool(True),
    useKtPruning = cms.bool(True),
    zcut = cms.double(0.),
    rcut_factor = cms.double(9999.),
    useDynamicFiltering = cms.bool(True),
    rFiltFactor = cms.double(0.5)
)
## CA pruned jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.caGenJetsNoNuPruned = ca4GenJets.clone(
    SubJetParameters,
    rParam = cms.double(options.jetRadius),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.caPFJetsCHSPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(options.jetRadius),
    src = process.caPFJetsCHS.src,
    srcPVs = process.caPFJetsCHS.srcPVs,
    doAreaFastjet = process.caPFJetsCHS.doAreaFastjet,
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20.)
)
## CA jets with Kt subjets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
## Kt subjets produced using Kt-based pruning with very loose pruning cuts (pruning is effectively disabled)
process.caGenJetsNoNuKtPruned = ca4GenJets.clone(
    SubJetParameters.clone(
        zcut = cms.double(0.),
        rcut_factor = cms.double(9999.)
    ),
    rParam = cms.double(options.jetRadius),
    src = cms.InputTag("genParticlesForJetsNoNu"),
    usePruning = cms.bool(True),
    useKtPruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.caPFJetsCHSKtPruned = ak5PFJetsPruned.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(options.jetRadius),
    src = process.caPFJetsCHS.src,
    srcPVs = process.caPFJetsCHS.srcPVs,
    doAreaFastjet = process.caPFJetsCHS.doAreaFastjet,
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(20.),
    useKtPruning = cms.bool(True),
    zcut = cms.double(0.),
    rcut_factor = cms.double(9999.)
)
## CA trimmed jets (Reco only)
from RecoJets.JetProducers.ak5PFJetsTrimmed_cfi import ak5PFJetsTrimmed
process.caPFJetsCHSTrimmed = ak5PFJetsTrimmed.clone(
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(options.jetRadius),
    src = process.caPFJetsCHS.src,
    srcPVs = process.caPFJetsCHS.srcPVs,
    doAreaFastjet = process.caPFJetsCHS.doAreaFastjet,
    jetPtMin = cms.double(20.)
)

#-------------------------------------
## PATify the above jets
from PhysicsTools.PatAlgos.tools.jetTools import *
## CA jets
switchJetCollection(process,
    jetCollection=cms.InputTag('caPFJetsCHS'),
    jetIdLabel='ca',
    rParam = options.jetRadius,
    useLegacyFlavour=False,
    doJTA=options.doJTA,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel = inputJetCorrLabelAK5,
    doType1MET   = False,
    genJetCollection = cms.InputTag("caGenJetsNoNu"),
    doJetID      = False,
)
## Filtered CA jets
addJetCollection(
    process,
    jetCollection=cms.InputTag('caPFJetsCHSFiltered'),
    algoLabel='CA',
    typeLabel='FilteredPFCHS',
    getJetMCFlavour=False,
    doJTA=False,
    doBTagging=False,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag("caGenJetsNoNu")
)
## Filtered subjets of CA jets
addJetCollection(
    process,
    jetCollection=cms.InputTag('caPFJetsCHSFiltered','SubJets'),
    algoLabel='CA',
    typeLabel='FilteredSubjetsPFCHS',
    rParam = options.jetRadius,
    useLegacyFlavour=False,
    doJTA=options.doJTA,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag('caGenJetsNoNuFiltered','SubJets')
)
## MassDrop-BDRS filtered CA jets
addJetCollection(
    process,
    jetCollection=cms.InputTag('caPFJetsCHSMDBDRSFiltered'),
    algoLabel='CA',
    typeLabel='MDBDRSFilteredPFCHS',
    getJetMCFlavour=False,
    doJTA=False,
    doBTagging=False,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag("caGenJetsNoNu")
)
## MassDrop-BDRS filtered subjets of CA jets
addJetCollection(
    process,
    jetCollection=cms.InputTag('caPFJetsCHSMDBDRSFiltered','SubJets'),
    algoLabel='CA',
    typeLabel='MDBDRSFilteredSubjetsPFCHS',
    rParam = options.jetRadius,
    useLegacyFlavour=False,
    doJTA=options.doJTA,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag('caGenJetsNoNuMDBDRSFiltered','SubJets')
)
## Kt-BDRS filtered CA jets
addJetCollection(
    process,
    jetCollection=cms.InputTag('caPFJetsCHSKtBDRSFiltered'),
    algoLabel='CA',
    typeLabel='KtBDRSFilteredPFCHS',
    getJetMCFlavour=False,
    doJTA=False,
    doBTagging=False,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag("caGenJetsNoNu")
)
## Kt-BDRS filtered subjets of CA jets
addJetCollection(
    process,
    jetCollection=cms.InputTag('caPFJetsCHSKtBDRSFiltered','SubJets'),
    algoLabel='CA',
    typeLabel='KtBDRSFilteredSubjetsPFCHS',
    rParam = options.jetRadius,
    useLegacyFlavour=False,
    doJTA=options.doJTA,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag('caGenJetsNoNuKtBDRSFiltered','SubJets')
)
## Pruned CA jets
addJetCollection(
    process,
    jetCollection=cms.InputTag('caPFJetsCHSPruned'),
    algoLabel='CA',
    typeLabel='PrunedPFCHS',
    getJetMCFlavour=False,
    doJTA=False,
    doBTagging=False,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag("caGenJetsNoNu")
)
## Pruned subjets of CA jets
addJetCollection(
    process,
    jetCollection=cms.InputTag('caPFJetsCHSPruned','SubJets'),
    algoLabel='CA',
    typeLabel='PrunedSubjetsPFCHS',
    rParam = options.jetRadius,
    useLegacyFlavour=False,
    doJTA=options.doJTA,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag('caGenJetsNoNuPruned','SubJets')
)
## Kt pruned CA jets
addJetCollection(
    process,
    jetCollection=cms.InputTag('caPFJetsCHSKtPruned'),
    algoLabel='CA',
    typeLabel='KtPrunedPFCHS',
    getJetMCFlavour=False,
    doJTA=False,
    doBTagging=False,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK7,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag("caGenJetsNoNu")
)
## Kt subjets of CA jets
addJetCollection(
    process,
    jetCollection=cms.InputTag('caPFJetsCHSKtPruned','SubJets'),
    algoLabel='CA',
    typeLabel='KtSubjetsPFCHS',
    rParam = options.jetRadius,
    useLegacyFlavour=False,
    doJTA=options.doJTA,
    doBTagging=options.doBTagging,
    btagInfo=bTagInfos,
    btagdiscriminators=bTagDiscriminators,
    jetCorrLabel=inputJetCorrLabelAK5,
    doType1MET=False,
    doL1Cleaning=False,
    doL1Counters=False,
    doJetID=False,
    genJetCollection=cms.InputTag('caGenJetsNoNuKtPruned','SubJets')
)

#-------------------------------------
## N-subjettiness
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

process.NjettinessCA = Njettiness.clone(
    src = cms.InputTag("caPFJetsCHS"),
    cone = cms.double(options.jetRadius)
)

process.patJets.userData.userFloats.src += ['NjettinessCA:tau1','NjettinessCA:tau2','NjettinessCA:tau3']

#-------------------------------------
## Grooming ValueMaps
from RecoJets.JetProducers.ca8PFJetsCHS_groomingValueMaps_cfi import ca8PFJetsCHSPrunedLinks

process.caPFJetsCHSPrunedMass = ca8PFJetsCHSPrunedLinks.clone(
    src = cms.InputTag("caPFJetsCHS"),
    matched = cms.InputTag("caPFJetsCHSPruned"),
    distMax = cms.double(options.jetRadius),
    value = cms.string('mass')
)

process.caPFJetsCHSFilteredMass = ca8PFJetsCHSPrunedLinks.clone(
    src = cms.InputTag("caPFJetsCHS"),
    matched = cms.InputTag("caPFJetsCHSFiltered"),
    distMax = cms.double(options.jetRadius),
    value = cms.string('mass')
)

process.caPFJetsCHSTrimmedMass = ca8PFJetsCHSPrunedLinks.clone(
    src = cms.InputTag("caPFJetsCHS"),
    matched = cms.InputTag("caPFJetsCHSTrimmed"),
    distMax = cms.double(options.jetRadius),
    value = cms.string('mass')
)

process.patJets.userData.userFloats.src += ['caPFJetsCHSPrunedMass','caPFJetsCHSFilteredMass','caPFJetsCHSTrimmedMass']

#-------------------------------------
if options.useSVClustering:
    ## Enable clustering-based jet-SV association for IVF vertices and AK5 jets
    #process.inclusiveSecondaryVertexFinderTagInfosAODPFlow = process.inclusiveSecondaryVertexFinderTagInfosAODPFlow.clone(
        #useSVClustering = cms.bool(True),
        ##useSVMomentum   = cms.bool(True), # otherwise using SV flight direction
        #jetAlgorithm    = cms.string("AntiKt"),
        #rParam          = cms.double(0.5),
        #ghostRescaling  = cms.double(1e-18)
    #)
    ## Enable clustering-based jet-SV association for IVF vertices and subjets of CA jets
    process.inclusiveSecondaryVertexFinderTagInfosCAFilteredSubjetsPFCHS = process.inclusiveSecondaryVertexFinderTagInfosCAFilteredSubjetsPFCHS.clone(
        useSVClustering = cms.bool(True),
        useSVMomentum   = cms.bool(options.useSVMomentum), # otherwise using SV flight direction
        jetAlgorithm    = cms.string("CambridgeAachen"),
        rParam          = cms.double(options.jetRadius),
        ghostRescaling  = cms.double(1e-18),
        fatJets         = cms.InputTag("caPFJetsCHS"),
        groomedFatJets   = cms.InputTag("caPFJetsCHSFiltered")
    )
    process.inclusiveSecondaryVertexFinderTagInfosCAMDBDRSFilteredSubjetsPFCHS = process.inclusiveSecondaryVertexFinderTagInfosCAMDBDRSFilteredSubjetsPFCHS.clone(
        useSVClustering = cms.bool(True),
        useSVMomentum   = cms.bool(options.useSVMomentum), # otherwise using SV flight direction
        jetAlgorithm    = cms.string("CambridgeAachen"),
        rParam          = cms.double(options.jetRadius),
        ghostRescaling  = cms.double(1e-18),
        fatJets         = cms.InputTag("caPFJetsCHS"),
        groomedFatJets   = cms.InputTag("caPFJetsCHSMDBDRSFiltered")
    )
    process.inclusiveSecondaryVertexFinderTagInfosCAKtBDRSFilteredSubjetsPFCHS = process.inclusiveSecondaryVertexFinderTagInfosCAKtBDRSFilteredSubjetsPFCHS.clone(
        useSVClustering = cms.bool(True),
        useSVMomentum   = cms.bool(options.useSVMomentum), # otherwise using SV flight direction
        jetAlgorithm    = cms.string("CambridgeAachen"),
        rParam          = cms.double(options.jetRadius),
        ghostRescaling  = cms.double(1e-18),
        fatJets         = cms.InputTag("caPFJetsCHS"),
        groomedFatJets   = cms.InputTag("caPFJetsCHSKtBDRSFiltered")
    )
    process.inclusiveSecondaryVertexFinderTagInfosCAPrunedSubjetsPFCHS = process.inclusiveSecondaryVertexFinderTagInfosCAPrunedSubjetsPFCHS.clone(
        useSVClustering = cms.bool(True),
        useSVMomentum   = cms.bool(options.useSVMomentum), # otherwise using SV flight direction
        jetAlgorithm    = cms.string("CambridgeAachen"),
        rParam          = cms.double(options.jetRadius),
        ghostRescaling  = cms.double(1e-18),
        fatJets         = cms.InputTag("caPFJetsCHS"),
        groomedFatJets   = cms.InputTag("caPFJetsCHSPruned")
    )
    process.inclusiveSecondaryVertexFinderTagInfosCAKtSubjetsPFCHS = process.inclusiveSecondaryVertexFinderTagInfosCAKtSubjetsPFCHS.clone(
        useSVClustering = cms.bool(True),
        useSVMomentum   = cms.bool(options.useSVMomentum), # otherwise using SV flight direction
        jetAlgorithm    = cms.string("CambridgeAachen"),
        rParam          = cms.double(options.jetRadius),
        ghostRescaling  = cms.double(1e-18),
        fatJets         = cms.InputTag("caPFJetsCHS"),
        groomedFatJets   = cms.InputTag("caPFJetsCHSKtPruned")
    )
#-------------------------------------
## New jet flavor still requires some cfg-level adjustments for subjets until it is better integrated into PAT
## Adjust the jet flavor for CA filtered subjets
process.patJetFlavourAssociationCAFilteredSubjetsPFCHS = process.patJetFlavourAssociation.clone(
    groomedJets = cms.InputTag("caPFJetsCHSFiltered"),
    subjets = cms.InputTag("caPFJetsCHSFiltered", "SubJets")
)
process.patJetsCAFilteredSubjetsPFCHS.JetFlavourInfoSource = cms.InputTag("patJetFlavourAssociationCAFilteredSubjetsPFCHS","SubJets")
## Adjust the jet flavor for CA MassDrop-BDRS filtered subjets
process.patJetFlavourAssociationCAMDBDRSFilteredSubjetsPFCHS = process.patJetFlavourAssociation.clone(
    groomedJets = cms.InputTag("caPFJetsCHSMDBDRSFiltered"),
    subjets = cms.InputTag("caPFJetsCHSMDBDRSFiltered", "SubJets")
)
process.patJetsCAMDBDRSFilteredSubjetsPFCHS.JetFlavourInfoSource = cms.InputTag("patJetFlavourAssociationCAMDBDRSFilteredSubjetsPFCHS","SubJets")
## Adjust the jet flavor for CA Kt-BDRS filtered subjets
process.patJetFlavourAssociationCAKtBDRSFilteredSubjetsPFCHS = process.patJetFlavourAssociation.clone(
    groomedJets = cms.InputTag("caPFJetsCHSKtBDRSFiltered"),
    subjets = cms.InputTag("caPFJetsCHSKtBDRSFiltered", "SubJets")
)
process.patJetsCAKtBDRSFilteredSubjetsPFCHS.JetFlavourInfoSource = cms.InputTag("patJetFlavourAssociationCAKtBDRSFilteredSubjetsPFCHS","SubJets")
## Adjust the jet flavor for CA pruned subjets
process.patJetFlavourAssociationCAPrunedSubjetsPFCHS = process.patJetFlavourAssociation.clone(
    groomedJets = cms.InputTag("caPFJetsCHSPruned"),
    subjets = cms.InputTag("caPFJetsCHSPruned", "SubJets")
)
process.patJetsCAPrunedSubjetsPFCHS.JetFlavourInfoSource = cms.InputTag("patJetFlavourAssociationCAPrunedSubjetsPFCHS","SubJets")
## Adjust the jet flavor for CA Kt subjets
process.patJetFlavourAssociationCAKtSubjetsPFCHS = process.patJetFlavourAssociation.clone(
    groomedJets = cms.InputTag("caPFJetsCHSKtPruned"),
    subjets = cms.InputTag("caPFJetsCHSKtPruned", "SubJets")
)
process.patJetsCAKtSubjetsPFCHS.JetFlavourInfoSource = cms.InputTag("patJetFlavourAssociationCAKtSubjetsPFCHS","SubJets")

#-------------------------------------
## Establish references between PATified fat jets and subjets using the BoostedJetMerger
process.selectedPatJetsCAFilteredPFCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsCAFilteredPFCHS"),
    subjetSrc=cms.InputTag("selectedPatJetsCAFilteredSubjetsPFCHS")
)

process.selectedPatJetsCAMDBDRSFilteredPFCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsCAMDBDRSFilteredPFCHS"),
    subjetSrc=cms.InputTag("selectedPatJetsCAMDBDRSFilteredSubjetsPFCHS")
)

process.selectedPatJetsCAKtBDRSFilteredPFCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsCAKtBDRSFilteredPFCHS"),
    subjetSrc=cms.InputTag("selectedPatJetsCAKtBDRSFilteredSubjetsPFCHS")
)

process.selectedPatJetsCAPrunedPFCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsCAPrunedPFCHS"),
    subjetSrc=cms.InputTag("selectedPatJetsCAPrunedSubjetsPFCHS")
)

process.selectedPatJetsCAKtPrunedPFCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsCAKtPrunedPFCHS"),
    subjetSrc=cms.InputTag("selectedPatJetsCAKtSubjetsPFCHS")
)

## Define BoostedJetMerger sequence
process.jetMergerSeq = cms.Sequence(
    process.selectedPatJetsCAFilteredPFCHSPacked
    + process.selectedPatJetsCAMDBDRSFilteredPFCHSPacked
    + process.selectedPatJetsCAKtBDRSFilteredPFCHSPacked
    + process.selectedPatJetsCAPrunedPFCHSPacked
    + process.selectedPatJetsCAKtPrunedPFCHSPacked
)

#-------------------------------------
from PhysicsTools.PatAlgos.tools.coreTools import *
## Keep only jets in the default sequence
removeAllPATObjectsBut(process, ['Jets'])

## Remove MC matching when running over data
if options.runOnData:
    removeMCMatching( process, ['All'] )

#-------------------------------------
## GenParticles for GenJets
from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu
process.genParticlesForJetsNoNu = genParticlesForJetsNoNu

#-------------------------------------
## Define instances of the RutgersJetAnalyzer
process.jetAnalyzerCAFatJets_PrunedSubjets = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('selectedPatJetsCAPrunedPFCHSPacked'),
    SubJetMode                = cms.string('Pruned'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    CalculateMassDrop         = cms.bool(True),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgIds           = cms.vint32(5)
)
process.jetAnalyzerCAFatJets_FilteredSubjets = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('selectedPatJetsCAFilteredPFCHSPacked'),
    SubJetMode                = cms.string('Filtered'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    CalculateMassDrop         = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgIds           = cms.vint32(5)
)
process.jetAnalyzerCAFatJets_MDBDRSFilteredSubjets = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('selectedPatJetsCAMDBDRSFilteredPFCHSPacked'),
    SubJetMode                = cms.string('Filtered'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    CalculateMassDrop         = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgIds           = cms.vint32(5)
)
process.jetAnalyzerCAFatJets_KtBDRSFilteredSubjets = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('selectedPatJetsCAKtBDRSFilteredPFCHSPacked'),
    SubJetMode                = cms.string('Filtered'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    CalculateMassDrop         = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgIds           = cms.vint32(5)
)
process.jetAnalyzerCAFatJets_KtSubjets = cms.EDAnalyzer('RutgersJetAnalyzer',
    UseEventWeight            = cms.bool(False),
    GenParticleTag            = cms.InputTag('genParticles'),
    JetsTag                   = cms.InputTag('selectedPatJets'),
    UseSubJets                = cms.bool(True),
    GroomedBasicJetsTag       = cms.InputTag('selectedPatJetsCAKtPrunedPFCHSPacked'),
    SubJetMode                = cms.string('Kt'),
    PvTag                     = cms.InputTag('goodOfflinePrimaryVertices'),
    JetRadius                 = cms.double(options.jetRadius),
    DoBosonMatching           = cms.bool(True),
    BosonMatchingRadius       = cms.double(0.5),
    BosonPdgId                = cms.int32(25),
    ApplyBosonIsolation       = cms.bool(True),
    DoBosonDecayProdSelection = cms.bool(True),
    BosonDecayProdPdgIds      = cms.vint32(5),
    CalculateMassDrop         = cms.bool(False),
    JetPtMin                  = cms.double(300.),
    JetPtBins                 = cms.uint32(3),
    JetPtBinWidth             = cms.double(200.),
    JetAbsEtaMax              = cms.double(1.5),
    JetMassMin                = cms.double(75.),
    JetMassMax                = cms.double(135.),
    DoJetFlavor               = cms.bool(False),
    JetFlavorPdgIds           = cms.vint32(5)
)

## Define jet analyzer sequence
process.jetAnalyzerSequence = cms.Sequence(
    process.jetAnalyzerCAFatJets_PrunedSubjets
    + process.jetAnalyzerCAFatJets_FilteredSubjets
    + process.jetAnalyzerCAFatJets_MDBDRSFilteredSubjets
    + process.jetAnalyzerCAFatJets_KtBDRSFilteredSubjets
    + process.jetAnalyzerCAFatJets_KtSubjets
)

#-------------------------------------
## If using explicit jet-track association
if options.useExplicitJTA:
    from RecoJets.JetAssociationProducers.ak5JTA_cff import ak5JetTracksAssociatorExplicit
    for m in getattr(process,"patDefaultSequence").moduleNames():
        if m.startswith('jetTracksAssociatorAtVertex'):
            print 'Switching ' + m + ' to explicit jet-track association'
            setattr( process, m, ak5JetTracksAssociatorExplicit.clone(jets = getattr(getattr(process,m),'jets')) )

#-------------------------------------
## Remove tau stuff that really shouldn't be there (probably a bug in PAT)
process.patDefaultSequence.remove(process.kt6PFJetsForRhoComputationVoronoi)
for m in getattr(process,"patDefaultSequence").moduleNames():
    if m.startswith('hpsPFTau'):
        getattr(process,"patDefaultSequence").remove(getattr(process,m))

#-------------------------------------
## Define jet sequences
process.genJetSeq = cms.Sequence(
    process.caGenJetsNoNu
    + process.caGenJetsNoNuFiltered
    + process.caGenJetsNoNuMDBDRSFiltered
    + process.caGenJetsNoNuKtBDRSFiltered
    + process.caGenJetsNoNuPruned
    + process.caGenJetsNoNuKtPruned
)
process.jetSeq = cms.Sequence(
    (
    process.caPFJetsCHS
    + process.caPFJetsCHSFiltered
    + process.caPFJetsCHSMDBDRSFiltered
    + process.caPFJetsCHSKtBDRSFiltered
    + process.caPFJetsCHSPruned
    + process.caPFJetsCHSKtPruned
    + process.caPFJetsCHSTrimmed
    )
    * (
    process.NjettinessCA
    + process.caPFJetsCHSFilteredMass
    + process.caPFJetsCHSPrunedMass
    + process.caPFJetsCHSTrimmedMass
    )
)

if not options.runOnData:
    process.jetSeq = cms.Sequence( process.genParticlesForJetsNoNu * process.genJetSeq + process.jetSeq )

#-------------------------------------
## Adapt primary vertex collection
from PhysicsTools.PatAlgos.tools.pfTools import *
adaptPVs(process, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), postfix='', sequence='jetSeq')
adaptPVs(process, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), postfix='', sequence='patDefaultSequence')

#-------------------------------------
## Add TagInfos to PAT jets
for m in ['patJets', 'patJetsCAFilteredSubjetsPFCHS', 'patJetsCAMDBDRSFilteredSubjetsPFCHS', 'patJetsCAKtBDRSFilteredSubjetsPFCHS',
          'patJetsCAPrunedSubjetsPFCHS', 'patJetsCAKtSubjetsPFCHS']:
    if hasattr(process,m) and getattr( getattr(process,m), 'addBTagInfo' ):
        print "Switching 'addTagInfos' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addTagInfos', cms.bool(True) )

#-------------------------------------
## Adapt fat jet b tagging
if options.doBTagging:
    # Set the cone size for the jet-track association to the jet radius
    process.jetTracksAssociatorAtVertex.coneSize = cms.double(options.jetRadius) # default is 0.5
    process.secondaryVertexTagInfosAOD.trackSelection.jetDeltaRMax = cms.double(options.jetRadius)   # default is 0.3
    process.secondaryVertexTagInfosAOD.vertexCuts.maxDeltaRToJetAxis = cms.double(options.jetRadius) # default is 0.5
    # Set the jet-SV dR to the jet radius
    process.inclusiveSecondaryVertexFinderTagInfosAOD.vertexCuts.maxDeltaRToJetAxis = cms.double(options.jetRadius) # default is 0.5
    process.inclusiveSecondaryVertexFinderTagInfosAOD.extSVDeltaRToJet = cms.double(options.jetRadius) # default is 0.3
    # Set the JP track dR cut to the jet radius
    process.jetProbabilityCA = process.jetProbability.clone( deltaR = cms.double(options.jetRadius) ) # default is 0.3
    process.jetProbabilityBJetTagsAOD.jetTagComputer = cms.string('jetProbabilityCA')
    # Set the JBP track dR cut to the jet radius
    process.jetBProbabilityCA = process.jetBProbability.clone( deltaR = cms.double(options.jetRadius) ) # default is 0.5
    process.jetBProbabilityBJetTagsAOD.jetTagComputer = cms.string('jetBProbabilityCA')
    # Set the CSV track dR cut to the jet radius
    process.combinedSecondaryVertexCA = process.combinedSecondaryVertex.clone()
    process.combinedSecondaryVertexCA.trackSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
    process.combinedSecondaryVertexCA.trackPseudoSelection.jetDeltaRMax = cms.double(options.jetRadius) # default is 0.3
    process.combinedSecondaryVertexBJetTagsAOD.jetTagComputer = cms.string('combinedSecondaryVertexCA')
    process.combinedInclusiveSecondaryVertexBJetTagsAOD.jetTagComputer = cms.string('combinedSecondaryVertexCA')

#-------------------------------------
## Various additional options
for m in getattr(process,'jetAnalyzerSequence').moduleNames():
    if not options.doBosonMatching:
        setattr( getattr(process,m), 'DoBosonMatching', cms.bool(False) )
    if not options.applyBosonIsolation:
        setattr( getattr(process,m), 'ApplyBosonIsolation', cms.bool(False) )
    if options.useEventWeight:
        setattr( getattr(process,m), 'UseEventWeight', cms.bool(True) )
    if options.useRadionCuts:
        setattr( getattr(process,m), 'JetPtMin',   cms.double(200.) )
        setattr( getattr(process,m), 'JetMassMin', cms.double(90.) )
        setattr( getattr(process,m), 'JetMassMax', cms.double(150.) )

#-------------------------------------
## Path definition
process.p = cms.Path(
    ( process.jetSeq * process.patDefaultSequence )
    * process.jetMergerSeq
    * process.jetAnalyzerSequence
)

## Delete output module
del process.out

## Schedule definition
process.schedule = cms.Schedule(process.p)

## Rewrite the cfg file
if options.dumpPythonCfg != '':
    open(options.dumpPythonCfg,'w').write(process.dumpPython())
