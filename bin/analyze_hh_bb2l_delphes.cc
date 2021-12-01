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

#include <TBenchmark.h> // TBenchmark
#include <TH1.h>        // TH1D
#include <TH2.h>        // TH2D

#include "hhAnalysis/DelphesAnalysis/interface/DelphesLepton.h"                   // DelphesLepton
#include "hhAnalysis/DelphesAnalysis/interface/DelphesLeptonReader.h"             // DelphesLeptonReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesLeptonCollectionSelector.h" // DelphesLeptonCollectionSelector
#include "hhAnalysis/DelphesAnalysis/interface/DelphesLeptonCollectionCleaner.h"  // DelphesLeptonCollectionCleaner
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJet.h"                      // DelphesJet
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetReader.h"                // DelphesJetReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetCollectionSelector.h"    // DelphesJetCollectionSelector
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetCollectionCleaner.h"     // DelphesJetCollectionCleaner
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetCalibration.h"           // DelphesJetCalibration
#include "hhAnalysis/DelphesAnalysis/interface/DelphesMEt.h"                      // DelphesMEt
#include "hhAnalysis/DelphesAnalysis/interface/DelphesMEtReader.h"                // DelphesMEtReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesEventInfo.h"                // DelphesEventInfo
#include "hhAnalysis/DelphesAnalysis/interface/DelphesEventInfoReader.h"          // DelphesEventInfoReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenParticle.h"              // DelphesGenParticle
#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenParticleReader.h"        // DelphesGenParticleReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenJet.h"                   // DelphesGenJet
#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenJetReader.h"             // DelphesGenJetReader
#include "hhAnalysis/DelphesAnalysis/interface/delphesAuxFunctions.h"             // mergeLeptonCollections, isBetterBJet, genBQuarkFromHiggs_or_TopDecay

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
#include <math.h>   // sqrt 

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
  bool isSignal = ( process_string.find("signal") != std::string::npos ) ? true : false;

  std::string histogramDir = cfg_analyze.getParameter<std::string>("histogramDir");

  std::string branchName_electrons = cfg_analyze.getParameter<std::string>("branchName_electrons");
  std::string branchName_muons = cfg_analyze.getParameter<std::string>("branchName_muons");
  std::string branchName_jets = cfg_analyze.getParameter<std::string>("branchName_jets");
  std::string branchName_met = cfg_analyze.getParameter<std::string>("branchName_met");

  std::string branchName_genParticles = cfg_analyze.getParameter<std::string>("branchName_genParticles");
  std::string branchName_genJets = cfg_analyze.getParameter<std::string>("branchName_genJets");
  std::string branchName_genMEt = cfg_analyze.getParameter<std::string>("branchName_genMEt");

  bool apply_genWeight = cfg_analyze.getParameter<bool>("apply_genWeight");

  bool isDEBUG = cfg_analyze.getParameter<bool>("isDEBUG");
  std::cout << " isDEBUG = " << isDEBUG << std::endl;

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
  DelphesLeptonReader* muonReader = new DelphesLeptonReader(DelphesLepton::kMuon, branchName_muons);
  inputTree->registerReader(muonReader);
  DelphesLeptonCollectionSelector muonSelector(Era::kUndefined, -1, isDEBUG);
  muonSelector.getSelector().set_min_pt(5.);
  muonSelector.getSelector().set_max_absEta(2.4);
  muonSelector.getSelector().set_max_relIsoRhoCorr(0.10);

  DelphesLeptonReader* electronReader = new DelphesLeptonReader(DelphesLepton::kElectron, branchName_electrons);
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
  DelphesJetCalibration jetCalibration;

  DelphesMEtReader* metReader = new DelphesMEtReader(branchName_met);
  inputTree->registerReader(metReader);

  DelphesGenParticleReader* genParticleReader = new DelphesGenParticleReader(branchName_genParticles);
  inputTree->registerReader(genParticleReader);

  DelphesGenJetReader* genJetReader = new DelphesGenJetReader(branchName_genJets);
  inputTree->registerReader(genJetReader);

  DelphesMEtReader* genMEtReader = new DelphesMEtReader(branchName_genMEt);
  inputTree->registerReader(genMEtReader);

//--- declare histograms showing multiplicity of generator-level b-jets
//    matching reconstructed b-jets,
//    for different b-tagging cuts applied on the reconstructed b-jets
  TFileDirectory dir = fs.mkdir(histogramDir);

  TH1* histogram_numGenBJets_geq1BJet  = dir.make<TH1D>("numGenBJets_geq1BJet",  "numGenBJets_geq1BJet",  4, -0.5, +3.5);
  TH1* histogram_numGenBJets_geq2BJets = dir.make<TH1D>("numGenBJets_geq2BJets", "numGenBJets_geq2BJets", 4, -0.5, +3.5);

  const int numPtBins = 12;
  double ptBinning[numPtBins + 1] = { 20., 25., 30., 35., 40., 50., 60., 80., 100., 120., 150., 200., 250. };
  TH2* histogram_jetEnResolution_b    = dir.make<TH2D>("jetEnResolution_b",    "jetEnResolution_b",    numPtBins, ptBinning, 200, -5., +5.);
  TH2* histogram_jetEnResolution_udsg = dir.make<TH2D>("jetEnResolution_udsg", "jetEnResolution_udsg", numPtBins, ptBinning, 200, -5., +5.);
  TH2* histogram_jetPtResolution_b    = dir.make<TH2D>("jetPtResolution_b",    "jetPtResolution_b",    numPtBins, ptBinning, 200, -5., +5.);
  TH2* histogram_jetPtResolution_udsg = dir.make<TH2D>("jetPtResolution_udsg", "jetPtResolution_udsg", numPtBins, ptBinning, 200, -5., +5.);
  
  TH1* histogram_metResolutionPx_uncalibrated = dir.make<TH1D>("metResolutionPx_uncalibrated", "metResolutionPx_uncalibrated", 200, -100., +100.);
  TH1* histogram_metResolutionPx_calibrated   = dir.make<TH1D>("metResolutionPx_calibrated",   "metResolutionPx_calibrated",   200, -100., +100.);
  TH1* histogram_metResolutionPy_uncalibrated = dir.make<TH1D>("metResolutionPy_uncalibrated", "metResolutionPy_uncalibrated", 200, -100., +100.);
  TH1* histogram_metResolutionPy_calibrated   = dir.make<TH1D>("metResolutionPy_calibrated",   "metResolutionPy_calibrated",   200, -100., +100.);

  TH1* histogram_mbb_uncalibrated_geq1BJet  = dir.make<TH1D>("mbb_uncalibrated_geq1BJet",  "mbb_uncalibrated_geq1BJet",  100, 0., 500.);
  TH1* histogram_mbb_calibrated_geq1BJet    = dir.make<TH1D>("mbb_calibrated_geq1BJet",    "mbb_calibrated_geq1BJet",    100, 0., 500.);
  TH1* histogram_mbb_uncalibrated_geq2BJets = dir.make<TH1D>("mbb_uncalibrated_geq2BJets", "mbb_uncalibrated_geq2BJets", 100, 0., 500.);
  TH1* histogram_mbb_calibrated_geq2BJets   = dir.make<TH1D>("mbb_calibrated_geq2BJets",   "mbb_calibrated_geq2BJets",   100, 0., 500.);
  TH1* histogram_mbb_uncalibrated_2genBJets = dir.make<TH1D>("mbb_uncalibrated_2genBJets", "mbb_uncalibrated_2genBJets", 100, 0., 500.);
  TH1* histogram_mbb_calibrated_2genBJets   = dir.make<TH1D>("mbb_calibrated_2genBJets",   "mbb_calibrated_2genBJets",   100, 0., 500.);
  TH1* histogram_mbb_genjets_2genBJets      = dir.make<TH1D>("mbb_genjets_2genBJets",      "mbb_genjets_2genBJets",      100, 0., 500.);
  TH1* histogram_mbb_bquarks_2genBJets      = dir.make<TH1D>("mbb_bquarks_2genBJets",      "mbb_bquarks_2genBJets",      100, 0., 500.);
  
  int analyzedEntries = 0;
  int selectedEntries = 0;
  double selectedEntries_weighted = 0.;
  TH1* histogram_analyzedEntries = dir.make<TH1D>("analyzedEntries", "analyzedEntries", 1, -0.5, +0.5);
  TH1* histogram_selectedEntries = dir.make<TH1D>("selectedEntries", "selectedEntries", 1, -0.5, +0.5);
  cutFlowTableType cutFlowTable;
  const edm::ParameterSet cutFlowTableCfg = makeHistManager_cfg(
    process_string, Form("%s/sel/cutFlow", histogramDir.data()), "", "central"
  );
  const std::vector<std::string> cuts = {
    "gen. Z-veto (signal only)",
    ">= 2 selected leptons",
    "lead lepton pT > 25 GeV && sublead lepton pT > 15 GeV",
    "lepton-pair OS charge",
    ">= 2 jets",
    ">= 1 b-jet",
    //"m(ll) > 12 GeV",
  };
  CutFlowTableHistManager * cutFlowHistManager = new CutFlowTableHistManager(cutFlowTableCfg, cuts);
  cutFlowHistManager->bookHistograms(dir);
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

//--- read generator-level information
    std::vector<DelphesGenParticle> genParticles = genParticleReader->read();
    if ( isDEBUG ) dumpCollection(genParticles);
    const std::vector<const DelphesGenParticle*> genParticles_ptrs = convert_to_ptrs(genParticles);

    std::vector<DelphesGenJet> genJets = genJetReader->read();
    if ( isDEBUG ) dumpCollection(genJets);
    const std::vector<const DelphesGenJet*> genJets_ptrs = convert_to_ptrs(genJets);

    DelphesMEt genMEt = genMEtReader->read();

//--- read collections of electrons and muons
//    resolve overlaps in order of priority: muon, electron,
    const std::vector<DelphesLepton> muons = muonReader->read();
    if ( isDEBUG ) dumpCollection(muons);
    const std::vector<const DelphesLepton*> muon_ptrs = convert_to_ptrs(muons);
    const std::vector<const DelphesLepton*> cleanedMuons = muon_ptrs; // CV: no cleaning needed for muons, as they have the highest priority in the overlap removal
    const std::vector<const DelphesLepton*> selMuons = muonSelector(cleanedMuons, isHigherPt<DelphesLepton>);

    const std::vector<DelphesLepton> electrons = electronReader->read();
    if ( isDEBUG ) dumpCollection(electrons);
    const std::vector<const DelphesLepton*> electron_ptrs = convert_to_ptrs(electrons);
    const std::vector<const DelphesLepton*> cleanedElectrons = electronCleaner(electron_ptrs, selMuons);
    const std::vector<const DelphesLepton*> selElectrons = electronSelector(cleanedElectrons, isHigherPt<DelphesLepton>);

    const std::vector<const DelphesLepton*> selLeptonsFull = mergeLeptonCollections(selElectrons, selMuons, isHigherPt<DelphesLepton>);
    const std::vector<const DelphesLepton*> selLeptons = pickFirstNobjects(selLeptonsFull, 2);

    //std::cout << "#leptons = " << (electrons.size() + muons.size()) << " (selected = " << selLeptons.size() << "):" << std::endl;
    //std::cout << " #electrons = " << electrons.size() << " (selected = " << selElectrons.size() << ")" << std::endl;
    //std::cout << " #muons = " << muons.size() << " (selected = " << selMuons.size() << ")" << std::endl;

//--- read collections of jets
    const std::vector<DelphesJet> jets = jetReader->read();
    if ( isDEBUG ) dumpCollection(jets);
    const std::vector<const DelphesJet*> jet_ptrs = convert_to_ptrs(jets);
    const std::vector<const DelphesJet*> cleanedJets = jetCleaner_dR04(jet_ptrs, selLeptons);
    const std::vector<const DelphesJet*> selJets = jetSelector(cleanedJets, isBetterBJet);

    DelphesMEt met = metReader->read();

//--- apply event selection
    // CV: remove HH->bbZZ events (keep only HH->bbWW events)    
    bool failsGenZVeto = false;
    if ( isSignal )
    {
      for ( std::vector<const DelphesGenParticle*>::const_iterator genParticle = genParticles_ptrs.begin();
            genParticle != genParticles_ptrs.end(); ++genParticle ) {
        if ( (*genParticle)->pdgId() == 23 )
        { 
          failsGenZVeto = true;
          break;
        }
      }
    }
    if ( failsGenZVeto ) {
      if ( isDEBUG ) {
        std::cout << "event " << eventInfo.str() << " FAILS generator-level Z-veto." << std::endl;
        printCollection("genParticles", genParticles_ptrs);
      }
      continue;
    }
    cutFlowTable.update("gen. Z-veto (signal only)", evtWeight);
    cutFlowHistManager->fillHistograms("gen. Z-veto (signal only)", evtWeight);

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

    const DelphesJet* selJet1_Hbb_uncalibrated = ( selJets.size() >= 1 ) ? selJets[0] : nullptr;
    const DelphesJet* selJet2_Hbb_uncalibrated = ( selJets.size() >= 2 ) ? selJets[1] : nullptr;
    if ( !(selJet1_Hbb_uncalibrated && selJet2_Hbb_uncalibrated) ) {
      if ( isDEBUG ) {
        std::cout << "event " << eventInfo.str() << " FAILS >= 2 jets selection\n";
      }
      continue;
    }
    cutFlowTable.update(">= 2 jets", evtWeight);
    cutFlowHistManager->fillHistograms(">= 2 jets", evtWeight);

    int numBJets = 0;
    if ( selJet1_Hbb_uncalibrated && selJet1_Hbb_uncalibrated->btag() >= 1 ) ++numBJets;
    if ( selJet2_Hbb_uncalibrated && selJet2_Hbb_uncalibrated->btag() >= 1 ) ++numBJets;
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

    //bool failsLowMassVeto = false;
    //for ( std::vector<const DelphesLepton*>::const_iterator selLepton1 = selLeptonsFull.begin();
    //      selLepton1 != selLeptonsFull.end(); ++selLepton1 ) {
    //  for ( std::vector<const DelphesLepton*>::const_iterator selLepton2 = selLepton1 + 1;
    //        selLepton2 != selLeptonsFull.end(); ++selLepton2 ) {
    //    if ( ((*selLepton1)->p4() + (*selLepton2)->p4()).mass() < 12. ) 
    //    {
    //      failsLowMassVeto = true;
    //    }
    //  }
    //}
    //if ( failsLowMassVeto ) {
    //  if ( isDEBUG ) {
    //    std::cout << "event " << eventInfo.str() << " FAILS low mass lepton pair veto." << std::endl;
    //  }
    //  continue;
    //}
    //cutFlowTable.update("m(ll) > 12 GeV", evtWeight);
    //cutFlowHistManager->fillHistograms("m(ll) > 12 GeV", evtWeight);

    int numGenBJets = 0;
    if ( selJet1_Hbb_uncalibrated && std::abs(selJet1_Hbb_uncalibrated->flavor()) == 5 ) ++numGenBJets;
    if ( selJet2_Hbb_uncalibrated && std::abs(selJet2_Hbb_uncalibrated->flavor()) == 5 ) ++numGenBJets;

    const double evtWeightErr = 0.;
 
    if ( numBJets >= 1 ) {
      fillWithOverFlow(histogram_numGenBJets_geq1BJet, numGenBJets, evtWeight, evtWeightErr);
    }
    if ( numBJets >= 2 ) {
      fillWithOverFlow(histogram_numGenBJets_geq2BJets, numGenBJets, evtWeight, evtWeightErr);
    }

    std::vector<const DelphesJet*> selJets_Hbb;
    if ( selJet1_Hbb_uncalibrated ) selJets_Hbb.push_back(selJet1_Hbb_uncalibrated);
    if ( selJet2_Hbb_uncalibrated ) selJets_Hbb.push_back(selJet2_Hbb_uncalibrated);
    for ( std::vector<const DelphesJet*>::const_iterator selJet = selJets_Hbb.begin();
          selJet != selJets_Hbb.end(); ++selJet ) {
      const DelphesGenJet* genJet_bestMatch = nullptr;
      double dR_bestMatch = 1.e+3;
      for ( std::vector<const DelphesGenJet*>::const_iterator genJet = genJets_ptrs.begin();
            genJet != genJets_ptrs.end(); ++genJet ) { 
        double dR = deltaR((*selJet)->p4(), (*genJet)->p4());
        // CV: require that jet directions agree within dR < 0.1 in order to match reconstructed to generator-level jets
        const double dRmax = 0.1;
        if ( dR < dRmax && dR < dR_bestMatch )
        {
          genJet_bestMatch = *genJet;
          dR_bestMatch = dR;
        }
      }
      if ( genJet_bestMatch )
      {
        bool isGenJet_b    =  std::abs((*selJet)->flavor()) == 5;
        bool isGenJet_udsg = (std::abs((*selJet)->flavor()) >= 1 && std::abs((*selJet)->flavor()) <= 3) || (*selJet)->flavor() == 21;
        TH2* histogram_jetEnResolution = nullptr;
        if      ( isGenJet_b    ) histogram_jetEnResolution = histogram_jetEnResolution_b;
        else if ( isGenJet_udsg ) histogram_jetEnResolution = histogram_jetEnResolution_udsg;
        if ( histogram_jetEnResolution )
        {
          float jetEnResolution = ((*selJet)->p4().energy() - genJet_bestMatch->p4().energy())/sqrt(std::max((double)1., genJet_bestMatch->p4().energy()));
          //std::cout << "E(rec) = " << (*selJet)->p4().energy() << ", E(gen) = " << genJet_bestMatch->p4().energy() << ":" 
          //          << " jetEnResolution = " << jetEnResolution << " (dR = " << dR_bestMatch << ")" << std::endl;
          if ( genJet_bestMatch->p4().energy() > 20. && genJet_bestMatch->p4().energy() < 250. && std::fabs(jetEnResolution) < 5. ) 
          {
            histogram_jetEnResolution->Fill(genJet_bestMatch->p4().energy(), jetEnResolution, evtWeight);
          }
        }
        TH2* histogram_jetPtResolution = nullptr;
        if      ( isGenJet_b    ) histogram_jetPtResolution = histogram_jetPtResolution_b;
        else if ( isGenJet_udsg ) histogram_jetPtResolution = histogram_jetPtResolution_udsg;
        if ( histogram_jetPtResolution )
        {
          float jetPtResolution = ((*selJet)->pt() - genJet_bestMatch->pt())/sqrt(std::max((float)1., genJet_bestMatch->pt()));
          //std::cout << "pT(rec) = " << (*selJet)->pt() << ", pT(gen) = " << genJet_bestMatch->pt() << ":" 
          //          << " jetPtResolution = " << jetPtResolution << " (dR = " << dR_bestMatch << ")" << std::endl;
          if ( genJet_bestMatch->pt() > 20. && genJet_bestMatch->pt() < 250. && std::fabs(jetPtResolution) < 5. ) 
          {
            histogram_jetPtResolution->Fill(genJet_bestMatch->pt(), jetPtResolution, evtWeight);
          }
        }
      }
    }

    DelphesJet selJet1_Hbb_calibrated = jetCalibration(*selJet1_Hbb_uncalibrated);
    DelphesJet selJet2_Hbb_calibrated = jetCalibration(*selJet2_Hbb_uncalibrated);

    if ( numBJets == 2 && numGenBJets == 2 ) {
      double metPx_uncalibrated = met.px();
      fillWithOverFlow(histogram_metResolutionPx_uncalibrated, metPx_uncalibrated - genMEt.px(), evtWeight, evtWeightErr);
      double metPx_calibrated = metPx_uncalibrated - (selJet1_Hbb_calibrated.p4().px() - selJet1_Hbb_uncalibrated->p4().px());
      fillWithOverFlow(histogram_metResolutionPx_calibrated, metPx_calibrated - genMEt.px(), evtWeight, evtWeightErr);
      double metPy_uncalibrated = met.py();
      fillWithOverFlow(histogram_metResolutionPy_uncalibrated, metPy_uncalibrated - genMEt.py(), evtWeight, evtWeightErr);
      double metPy_calibrated = metPy_uncalibrated - (selJet1_Hbb_calibrated.p4().py() - selJet1_Hbb_uncalibrated->p4().py());
      fillWithOverFlow(histogram_metResolutionPy_calibrated, metPy_calibrated - genMEt.py(), evtWeight, evtWeightErr);
    }

    double mbb_uncalibrated = (selJet1_Hbb_uncalibrated->p4() + selJet2_Hbb_uncalibrated->p4()).mass();
    double mbb_calibrated   = (selJet1_Hbb_calibrated.p4() + selJet2_Hbb_calibrated.p4()).mass();
    if ( numBJets    >= 1 ) {
      fillWithOverFlow(histogram_mbb_uncalibrated_geq1BJet,  mbb_uncalibrated, evtWeight, evtWeightErr);
      fillWithOverFlow(histogram_mbb_calibrated_geq1BJet,    mbb_calibrated,   evtWeight, evtWeightErr);
    }
    if ( numBJets    >= 2 ) {
      fillWithOverFlow(histogram_mbb_uncalibrated_geq2BJets, mbb_uncalibrated, evtWeight, evtWeightErr);
      fillWithOverFlow(histogram_mbb_calibrated_geq2BJets,   mbb_calibrated,   evtWeight, evtWeightErr);
    }
    if ( numGenBJets == 2 ) {
      fillWithOverFlow(histogram_mbb_uncalibrated_2genBJets, mbb_uncalibrated, evtWeight, evtWeightErr);
      fillWithOverFlow(histogram_mbb_calibrated_2genBJets,   mbb_calibrated,   evtWeight, evtWeightErr);
      
      const DelphesGenJet* genJet1 = get_genJet(selJet1_Hbb_uncalibrated->p4(), genJets);
      const DelphesGenJet* genJet2 = get_genJet(selJet2_Hbb_uncalibrated->p4(), genJets);
      if ( genJet1 && genJet2 ) {
        double mbb_genjets = (genJet1->p4() + genJet2->p4()).mass();
        fillWithOverFlow(histogram_mbb_genjets_2genBJets,    mbb_genjets,      evtWeight, evtWeightErr);
      }

      int mother_pdgId = ( isSignal ) ? 25 : 6; 
      const DelphesGenParticle* genBQuark1 = get_genBQuarkFromHiggs_or_TopDecay(selJet1_Hbb_uncalibrated->p4(), genParticles, mother_pdgId);
      const DelphesGenParticle* genBQuark2 = get_genBQuarkFromHiggs_or_TopDecay(selJet2_Hbb_uncalibrated->p4(), genParticles, mother_pdgId);
      if ( genBQuark1 && genBQuark2 ) {
        double mbb_bquarks = (genBQuark1->p4() + genBQuark2->p4()).mass();
        fillWithOverFlow(histogram_mbb_bquarks_2genBJets,   mbb_bquarks,       evtWeight, evtWeightErr);
      }
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
  delete genParticleReader;
  delete genJetReader;
  delete genMEtReader;

  delete inputTree;

  clock.Show("analyze_hh_bb2l_delphes");

  return EXIT_SUCCESS;
}
