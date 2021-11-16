#include "FWCore/ParameterSet/interface/ParameterSet.h" // edm::ParameterSet
#include "FWCore/Utilities/interface/Exception.h"       // cms::Exception
#include "PhysicsTools/FWLite/interface/TFileService.h" // fwlite::TFileService
#include "DataFormats/FWLite/interface/InputSource.h"   // fwlite::InputSource
#include "DataFormats/FWLite/interface/OutputFiles.h"   // fwlite::OutputFiles
#include "DataFormats/Math/interface/LorentzVector.h"   // math::PtEtaPhiMLorentzVector
#include "DataFormats/Math/interface/deltaR.h"          // deltaR
#include "DataFormats/Math/interface/deltaPhi.h"        // deltaPhi

#if __has_include (<FWCore/ParameterSetReader/interface/ParameterSetReader.h>)
#  include <FWCore/ParameterSetReader/interface/ParameterSetReader.h> // edm::readPSetsFrom()
#else
#  include <FWCore/PythonParameterSet/interface/MakeParameterSets.h>  // edm::readPSetsFrom()
#endif

#include <TBenchmark.h>     // TBenchmark
#include <TString.h>        // TString, Form
#include <TError.h>         // gErrorAbortLevel, kError
#include <TRandom3.h>       // TRandom3
#include <TLorentzVector.h> // TLorentzVector 
#include <TMatrixD.h>       // TMatrixD

#include "hhAnalysis/DelphesAnalysis/interface/DelphesLepton.h"                   // DelphesLepton
#include "hhAnalysis/DelphesAnalysis/interface/DelphesLeptonReader.h"             // DelphesLeptonReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesLeptonCollectionSelector.h" // DelphesLeptonCollectionSelector
#include "hhAnalysis/DelphesAnalysis/interface/DelphesLeptonCollectionCleaner.h"  // DelphesLeptonCollectionCleaner
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJet.h"                      // DelphesJet
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetReader.h"                // DelphesJetReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetCollectionSelector.h"    // DelphesJetCollectionSelector
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetCollectionCleaner.h"     // DelphesJetCollectionCleaner
#include "hhAnalysis/DelphesAnalysis/interface/DelphesMEt.h"                      // DelphesMEt
#include "hhAnalysis/DelphesAnalysis/interface/DelphesMEtReader.h"                // DelphesMEtReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesEventInfo.h"                // DelphesEventInfo
#include "hhAnalysis/DelphesAnalysis/interface/DelphesEventInfoReader.h"          // DelphesEventInfoReader

#include "tthAnalysis/HiggsToTauTau/interface/TTreeWrapper.h"                     // TTreeWrapper
#include "tthAnalysis/HiggsToTauTau/interface/histogramAuxFunctions.h"            // fillWithOverFlow
#include "tthAnalysis/HiggsToTauTau/interface/cutFlowTable.h"                     // cutFlowTableType
#include "tthAnalysis/HiggsToTauTau/interface/HistManagerBase.h"                  // makeHistManager_cfg
#include "tthAnalysis/HiggsToTauTau/interface/CutFlowTableHistManager.h"          // CutFlowTableHistManager
#include "tthAnalysis/HiggsToTauTau/interface/convert_to_ptrs.h"                  // convert_to_ptrs

#include <iostream> // std::cerr, std::fixed
#include <iomanip>  // std::setprecision(), std::setw()
#include <string>   // std::string
#include <vector>   // std::vector<>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <fstream>  // std::ofstream
#include <assert.h> // assert

template <typename T>
bool
isHigherPt(const T * particle1,
           const T * particle2)
{
  return particle1->pt() > particle2->pt();
}

std::vector<const DelphesLepton *>
mergeLeptonCollections(const std::vector<const DelphesLepton *> & electrons,
                       const std::vector<const DelphesLepton *> & muons,
                       bool (*sortFunction)(const DelphesLepton *, const DelphesLepton *))
{
  std::vector<const DelphesLepton *> leptons;
  const std::size_t numLeptons = electrons.size() + muons.size();
  if ( numLeptons > 0 )
  {
    leptons.reserve(numLeptons);
    leptons.insert(leptons.end(), electrons.begin(), electrons.end());
    leptons.insert(leptons.end(), muons.begin(), muons.end());
    std::sort(leptons.begin(), leptons.end(), sortFunction);
  }
  return leptons;
}

bool 
isBetterBJet(const DelphesJet * jet1, const DelphesJet * jet2)
{
  // CV: select jets passing b-tag discriminator first,
  //     then select jets by higher pT
  if ( jet1->btag() > jet2->btag() ) return true;
  if ( jet1->btag() < jet2->btag() ) return false;
  return jet1->pt() > jet2->pt();
}

/**
 * @brief Check validity (with Delphes) of assumptions on detector performance that we make in HH->bbWW MEM paper
 */
int main(int argc, char* argv[])
{
//--- throw an exception in case ROOT encounters an error
  gErrorAbortLevel = kError;

//--- stop ROOT from keeping track of all histograms
  TH1::AddDirectory(false);

//--- parse command-line arguments
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "<analyze_hh_bb2l_delphes>:" << std::endl;

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("analyze_hh_bb2l_delphes");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") )
    throw cms::Exception("analyze_hh_bb2l_delphes")
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfg_analyze = cfg.getParameter<edm::ParameterSet>("analyze_hh_bb2l_delphes");

  std::string treeName = cfg_analyze.getParameter<std::string>("treeName");

  std::string process_string = cfg_analyze.getParameter<std::string>("process");
  //bool isSignal = ( process_string.find("signal") != std::string::npos ) ? true : false;

  std::string histogramDir = cfg_analyze.getParameter<std::string>("histogramDir");

  std::string branchName_electrons = cfg_analyze.getParameter<std::string>("branchName_electrons");
  std::string branchName_muons = cfg_analyze.getParameter<std::string>("branchName_muons");
  std::string branchName_jets = cfg_analyze.getParameter<std::string>("branchName_jets");
  std::string branchName_met = cfg_analyze.getParameter<std::string>("branchName_met");

  std::string branchName_genMEt = cfg_analyze.getParameter<std::string>("branchName_genMEt");

  bool apply_genWeight = cfg_analyze.getParameter<bool>("apply_genWeight");

  bool isDEBUG = cfg_analyze.getParameter<bool>("isDEBUG");

  fwlite::InputSource inputFiles(cfg);
  int maxEvents = inputFiles.maxEvents();
  std::cout << " maxEvents = " << maxEvents << std::endl;
  unsigned reportEvery = inputFiles.reportAfter();

  fwlite::OutputFiles outputFile(cfg);
  fwlite::TFileService fs = fwlite::TFileService(outputFile.file().data());

  TTreeWrapper* inputTree = new TTreeWrapper(treeName.data(), inputFiles.files(), maxEvents);
  std::cout << "Loaded " << inputTree->getFileCount() << " file(s)." << std::endl;

//--- declare event-level variables
  DelphesEventInfo eventInfo;
  DelphesEventInfoReader eventInfoReader(&eventInfo);
  inputTree->registerReader(&eventInfoReader);

//--- declare collections of reconstructed particles
  DelphesLeptonReader* muonReader = new DelphesLeptonReader(branchName_muons);
  inputTree->registerReader(muonReader);
  DelphesLeptonCollectionSelector muonSelector(Era::kUndefined, -1, isDEBUG);
  muonSelector.getSelector().set_min_pt(5.);
  muonSelector.getSelector().set_max_absEta(2.4);
  muonSelector.getSelector().set_max_relIsoRhoCorr(0.10);

  DelphesLeptonReader* electronReader = new DelphesLeptonReader(branchName_electrons);
  inputTree->registerReader(electronReader);
  DelphesLeptonCollectionCleaner electronCleaner(0.3, isDEBUG);
  DelphesLeptonCollectionSelector electronSelector(Era::kUndefined, -1, isDEBUG);
  electronSelector.getSelector().set_min_pt(7.);
  electronSelector.getSelector().set_max_absEta(2.5);
  electronSelector.getSelector().set_max_relIsoRhoCorr(0.10);

  DelphesJetReader* jetReader = new DelphesJetReader(branchName_jets);
  inputTree->registerReader(jetReader);
  DelphesJetCollectionCleaner jetCleaner_dR04(0.4, isDEBUG);
  DelphesJetCollectionSelector jetSelector(Era::kUndefined, -1, isDEBUG);

  DelphesMEtReader* metReader = new DelphesMEtReader(branchName_met);
  inputTree->registerReader(metReader);

  DelphesMEtReader* genMEtReader = new DelphesMEtReader(branchName_genMEt);
  inputTree->registerReader(genMEtReader);

//--- declare histograms showing multiplicity of generator-level b-jets
//    matching reconstructed b-jets,
//    for different b-tagging cuts applied on the reconstructed b-jets
  TH1* histogram_numGenBJets_geq1BJet  = fs.make<TH1D>("numGenBJets_geq1BJet",  "numGenBJets_geq1BJet",  4, -0.5, +3.5);
  TH1* histogram_numGenBJets_geq2BJets = fs.make<TH1D>("numGenBJets_geq2BJets", "numGenBJets_geq2BJets", 4, -0.5, +3.5);

  TH1* histogram_metResPx = fs.make<TH1D>("metResPx", "metResPx", 200, -100., +100.);
  TH1* histogram_metResPy = fs.make<TH1D>("metResPy", "metResPy", 200, -100., +100.);

  int analyzedEntries = 0;
  int selectedEntries = 0;
  double selectedEntries_weighted = 0.;
  TH1* histogram_analyzedEntries = fs.make<TH1D>("analyzedEntries", "analyzedEntries", 1, -0.5, +0.5);
  TH1* histogram_selectedEntries = fs.make<TH1D>("selectedEntries", "selectedEntries", 1, -0.5, +0.5);
  cutFlowTableType cutFlowTable;
  const edm::ParameterSet cutFlowTableCfg = makeHistManager_cfg(
    process_string, Form("%s/sel/cutFlow", histogramDir.data()), "", "central"
  );
  const std::vector<std::string> cuts = {
    ">= 2 selected leptons",
    "lead lepton pT > 25 GeV && sublead lepton pT > 15 GeV",
    "lepton-pair OS charge",
    ">= 2 jets from H->bb",
    ">= 1 b-jet",
    "m(ll) > 12 GeV",
  };
  CutFlowTableHistManager * cutFlowHistManager = new CutFlowTableHistManager(cutFlowTableCfg, cuts);
  cutFlowHistManager->bookHistograms(fs);
  while ( inputTree->hasNextEvent() ) {
    if ( inputTree -> canReport(reportEvery) ) {
      std::cout << "processing Entry " << inputTree -> getCurrentMaxEventIdx()
                << " or " << inputTree -> getCurrentEventIdx() << " entry in #"
                << (inputTree -> getProcessedFileCount() - 1)
                << " (" << eventInfo
                << ") file (" << selectedEntries << " Entries selected)\n";
    }
    ++analyzedEntries;
    histogram_analyzedEntries->Fill(0.);

    if ( isDEBUG ) {
      std::cout << "event #" << inputTree -> getCurrentMaxEventIdx() << ' ' << eventInfo << '\n';
    }

    double evtWeight = 1.;
    if ( apply_genWeight ) evtWeight *= eventInfo.genweight();

//--- build collections of electrons and muons
//    resolve overlaps in order of priority: muon, electron,
    const std::vector<DelphesLepton> muons = muonReader->read();
    const std::vector<const DelphesLepton*> muon_ptrs = convert_to_ptrs(muons);
    const std::vector<const DelphesLepton*> cleanedMuons = muon_ptrs; // CV: no cleaning needed for muons, as they have the highest priority in the overlap removal
    const std::vector<const DelphesLepton*> selMuons = muonSelector(cleanedMuons, isHigherPt<DelphesLepton>);

    const std::vector<DelphesLepton> electrons = electronReader->read();
    const std::vector<const DelphesLepton*> electron_ptrs = convert_to_ptrs(electrons);
    const std::vector<const DelphesLepton*> cleanedElectrons = electronCleaner(electron_ptrs, selMuons);
    const std::vector<const DelphesLepton*> selElectrons = electronSelector(cleanedElectrons, isHigherPt<DelphesLepton>);

    const std::vector<const DelphesLepton*> selLeptonsFull = mergeLeptonCollections(selElectrons, selMuons, isHigherPt<DelphesLepton>);
    const std::vector<const DelphesLepton*> selLeptons = pickFirstNobjects(selLeptonsFull, 2);

//--- build collections of jets
    const std::vector<DelphesJet> jets = jetReader->read();
    const std::vector<const DelphesJet*> jet_ptrs = convert_to_ptrs(jets);
    const std::vector<const DelphesJet*> cleanedJets = jetCleaner_dR04(jet_ptrs, selLeptons)
    ;
    const std::vector<const DelphesJet*> selJets = jetSelector(cleanedJets, isBetterBJet);

    // require at least two leptons passing fakeable selection criteria
    if ( !(selLeptons.size() >= 2) ) {
      if ( isDEBUG ) {
        std::cout << "event " << eventInfo.str() << " FAILS selLeptons selection." << std::endl;
        printCollection("selLeptons", selLeptons);
      }
      continue;
    }
    cutFlowTable.update(">= 2 selected leptons", evtWeight);
    cutFlowHistManager->fillHistograms(">= 2 selected leptons", evtWeight);

    const DelphesLepton* selLepton_lead = selLeptons[0];
    const DelphesLepton* selLepton_sublead = selLeptons[1];

    const double minPt_lead = 25.;
    const double minPt_sublead = 15.;
    if ( !(selLepton_lead->pt() > minPt_lead && selLepton_sublead->pt() > minPt_sublead) ) {
      if ( isDEBUG ) {
        std::cout << "event " << eventInfo.str() << " FAILS lepton pT selection." << std::endl;
        std::cout << " (leading selLepton pT = " << selLepton_lead->pt() << ", minPt_lead = " << minPt_lead
                  << ", subleading selLepton pT = " << selLepton_sublead->pt() << ", minPt_sublead = " << minPt_sublead << ")" << std::endl;
      }
      continue;
    }
    cutFlowTable.update("lead lepton pT > 25 GeV && sublead lepton pT > 15 GeV", evtWeight);
    cutFlowHistManager->fillHistograms("lead lepton pT > 25 GeV && sublead lepton pT > 15 GeV", evtWeight);

    if ( !(selLepton_lead->charge()*selLepton_sublead->charge() < 0) ) {
      if ( isDEBUG ) {
	std::cout << "event " << eventInfo.str() << " FAILS lepton charge selection." << std::endl;
        std::cout << " (leading selLepton charge = " << selLepton_lead->charge()
                  << ", subleading selLepton charge = " << selLepton_sublead->charge() << ", leptonChargeSelection = OS)" << std::endl;
      }
      continue;
    }
    cutFlowTable.update("lepton-pair OS charge", evtWeight);
    cutFlowHistManager->fillHistograms("lepton-pair OS charge", evtWeight);

    const DelphesJet* selJet1_Hbb = ( selJets.size() >= 1 ) ? selJets[0] : nullptr;
    const DelphesJet* selJet2_Hbb = ( selJets.size() >= 2 ) ? selJets[1] : nullptr;
    if ( !(selJet1_Hbb && selJet2_Hbb) ) {
      if ( isDEBUG ) {
        std::cout << "event " << eventInfo.str() << " FAILS >= 2 jets from H->bb selection\n";
      }
      continue;
    }
    cutFlowTable.update(">= 2 jets from H->bb", evtWeight);
    cutFlowHistManager->fillHistograms(">= 2 jets from H->bb", evtWeight);

    int numBJets = 0;
    if ( selJet1_Hbb && selJet1_Hbb->btag() >= 1 ) ++numBJets;
    if ( selJet2_Hbb && selJet2_Hbb->btag() >= 1 ) ++numBJets;
    if ( !(numBJets >= 1) )
    {
      if ( isDEBUG ) 
      {
        std::cout << "event " << eventInfo.str() << " FAILS >= 1 b-jet selection\n";
      }
      continue;
    }
    cutFlowTable.update(">= 1 b-jet", evtWeight);
    cutFlowHistManager->fillHistograms(">= 1 b-jet", evtWeight);

    bool failsLowMassVeto = false;
    for ( std::vector<const DelphesLepton*>::const_iterator selLepton1 = selLeptonsFull.begin();
          selLepton1 != selLeptonsFull.end(); ++selLepton1 ) {
      for ( std::vector<const DelphesLepton*>::const_iterator selLepton2 = selLepton1 + 1;
            selLepton2 != selLeptonsFull.end(); ++selLepton2 ) {
        if ( ((*selLepton1)->p4() + (*selLepton2)->p4()).mass() < 12. ) 
        {
          failsLowMassVeto = true;
        }
      }
    }
    if ( failsLowMassVeto ) {
      if ( isDEBUG ) {
        std::cout << "event " << eventInfo.str() << " FAILS low mass lepton pair veto." << std::endl;
      }
      continue;
    }
    cutFlowTable.update("m(ll) > 12 GeV", evtWeight);
    cutFlowHistManager->fillHistograms("m(ll) > 12 GeV", evtWeight);

    int numGenBJets = 0;
    if ( selJet1_Hbb && std::abs(selJet1_Hbb->flavor()) == 5 ) ++numGenBJets;
    if ( selJet2_Hbb && std::abs(selJet2_Hbb->flavor()) == 5 ) ++numGenBJets;

    DelphesMEt met = metReader->read();

    DelphesMEt genMEt = genMEtReader->read();

    const double evtWeightErr = 0.;
 
    if ( numBJets >= 1 ) {
      fillWithOverFlow(histogram_numGenBJets_geq1BJet, numGenBJets, evtWeight, evtWeightErr);
    }
    if ( numBJets >= 2 ) {
      fillWithOverFlow(histogram_numGenBJets_geq2BJets, numGenBJets, evtWeight, evtWeightErr);
    }

    if ( numBJets == 2 && numGenBJets == 2 ) {
      fillWithOverFlow(histogram_metResPx, met.px() - genMEt.px(), evtWeight, evtWeightErr);
      fillWithOverFlow(histogram_metResPy, met.py() - genMEt.py(), evtWeight, evtWeightErr);
    }

    ++selectedEntries;
    selectedEntries_weighted += evtWeight;
    histogram_selectedEntries->Fill(0.);
  }

  std::cout << "max num. Entries = " << inputTree -> getCumulativeMaxEventCount()
            << " (limited by " << maxEvents << ") processed in "
            << inputTree -> getProcessedFileCount() << " file(s) (out of "
            << inputTree -> getFileCount() << ")\n"
            << " analyzed = " << analyzedEntries << '\n'
            << " selected = " << selectedEntries << " (weighted = " << selectedEntries_weighted << ")\n\n"
            << "cut-flow table" << std::endl;
  cutFlowTable.print(std::cout);
  std::cout << std::endl;

  delete muonReader;
  delete electronReader;
  delete jetReader;
  delete metReader;
  delete genMEtReader;

  delete inputTree;

  clock.Show("analyze_hh_bb2l_delphes");

  return EXIT_SUCCESS;
}
