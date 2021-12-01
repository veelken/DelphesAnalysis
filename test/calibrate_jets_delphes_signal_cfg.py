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

process.calibrate_jets_delphes = cms.PSet(
    treeName = cms.string('Events'),

    process = cms.string(''),
    histogramDir = cms.string(''),

    apply_genWeight = cms.bool(True),

    branchName_electrons = cms.string('Electron'),
    branchName_muons = cms.string('Muon'),
    branchName_jets = cms.string('Jet'),

    branchName_genParticles = cms.string('GenPart'),
    branchName_genJets = cms.string('GenJet'),

    isDEBUG = cms.bool(False)
)

process.fwliteOutput.fileName = cms.string('calibrate_jets_delphes_signal.root')

process.calibrate_jets_delphes.process = cms.string('signal')
process.calibrate_jets_delphes.histogramDir = cms.string('signal')

inputFilePath = "/hdfs/local/karl/DelphesNtuples/2016/HH_DL_LO_PU40/0000/"

import re
inputFile_regex = r"tree_[a-zA-Z0-9-_]+.root"
inputFile_matcher = re.compile(inputFile_regex)

from hhAnalysis.DelphesAnalysis.tools.jobTools import getInputFileNames
inputFileNames = getInputFileNames(inputFilePath, inputFile_matcher)
process.fwliteInput.fileNames = cms.vstring(inputFileNames)
