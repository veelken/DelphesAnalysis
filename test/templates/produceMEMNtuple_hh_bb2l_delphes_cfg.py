import FWCore.ParameterSet.Config as cms
import os

#--------------------------------------------------------------------------------
# CV: imports needed by analyzeConfig.py base-class
from tthAnalysis.HiggsToTauTau.configs.recommendedMEtFilters_cfi import *
from tthAnalysis.HiggsToTauTau.configs.EvtYieldHistManager_cfi import *
from tthAnalysis.HiggsToTauTau.configs.hhWeight_cfi import hhWeight
#--------------------------------------------------------------------------------

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring(),
    maxEvents = cms.int32(-1),
    outputEvery = cms.uint32(10000)
)

process.fwliteOutput = cms.PSet(
    fileName = cms.string('')
)

process.analyze_hh_bbwwMEM_dilepton = cms.PSet(
    treeName = cms.string('Events'),

    skipSelEvents = cms.int32(0),
    maxSelEvents = cms.int32(1000),

    process = cms.string(''),
    histogramDir = cms.string(''),

    apply_genWeight = cms.bool(True),

    branchName_electrons = cms.string('Electron'),
    branchName_muons = cms.string('Muon'),
    branchName_jets = cms.string('Jet'),
    branchName_met = cms.string('MET'),

    branchName_genParticles = cms.string('GenPart'),
    branchName_genJets = cms.string('GenJet'),
    branchName_genMEt = cms.string('GenMET'),

    isDEBUG = cms.bool(False)
)

#--------------------------------------------------------------------------------
# CV: parameter sets needed by analyzeConfig.py base-class
process.analyze_hh_bbwwMEM_dilepton.leptonFakeRateWeight = cms.PSet()
process.analyze_hh_bbwwMEM_dilepton.hhWeight_cfg = cms.PSet()
#--------------------------------------------------------------------------------
