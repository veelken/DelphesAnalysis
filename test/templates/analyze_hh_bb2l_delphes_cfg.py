import FWCore.ParameterSet.Config as cms
import os

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames = cms.vstring(),
    maxEvents = cms.int32(-1),
    outputEvery = cms.uint32(10000)
)

process.fwliteOutput = cms.PSet(
    fileName = cms.string('')
)

process.analyze_hh_bb2l_delphes = cms.PSet(
    treeName = cms.string('Events'),

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
