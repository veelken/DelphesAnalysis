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
#include <TRandom3.h>   // TRandom3
#include <TMatrixD.h>   // TMatrixD

#include "hhAnalysis/DelphesAnalysis/interface/DelphesLepton.h"                           // DelphesLepton
#include "hhAnalysis/DelphesAnalysis/interface/DelphesLeptonReader.h"                     // DelphesLeptonReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesLeptonCollectionSelector.h"         // DelphesLeptonCollectionSelector
#include "hhAnalysis/DelphesAnalysis/interface/DelphesLeptonCollectionCleaner.h"          // DelphesLeptonCollectionCleaner
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJet.h"                              // DelphesJet
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetReader.h"                        // DelphesJetReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetCollectionSelector.h"            // DelphesJetCollectionSelector
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetCollectionCleaner.h"             // DelphesJetCollectionCleaner
#include "hhAnalysis/DelphesAnalysis/interface/DelphesMEt.h"                              // DelphesMEt
#include "hhAnalysis/DelphesAnalysis/interface/DelphesMEtReader.h"                        // DelphesMEtReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesEventInfo.h"                        // DelphesEventInfo
#include "hhAnalysis/DelphesAnalysis/interface/DelphesEventInfoReader.h"                  // DelphesEventInfoReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenParticle.h"                      // DelphesGenParticle
#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenParticleReader.h"                // DelphesGenParticleReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenJet.h"                           // DelphesGenJet
#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenJetReader.h"                     // DelphesGenJetReader

#include "tthAnalysis/HiggsToTauTau/interface/TTreeWrapper.h"                             // TTreeWrapper
#include "tthAnalysis/HiggsToTauTau/interface/histogramAuxFunctions.h"                    // fillWithOverFlow
#include "tthAnalysis/HiggsToTauTau/interface/cutFlowTable.h"                             // cutFlowTableType
#include "tthAnalysis/HiggsToTauTau/interface/HistManagerBase.h"                          // makeHistManager_cfg
#include "tthAnalysis/HiggsToTauTau/interface/CutFlowTableHistManager.h"                  // CutFlowTableHistManager
#include "tthAnalysis/HiggsToTauTau/interface/convert_to_ptrs.h"                          // convert_to_ptrs

#include "hhAnalysis/bbwwMEM/interface/MEMbbwwAlgoDilepton.h"                             // MEMbbwwAlgoDilepton
#include "hhAnalysis/bbwwMEM/interface/MeasuredParticle.h"                                // MeasuredParticle
#include "hhAnalysis/bbwwMEM/interface/memAuxFunctions.h"                                 // mem::electronMass, mem::muonMass, mem::square

#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMEvent_dilepton.h"             // MEMEvent_dilepton
#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/MEMbbwwNtupleManager_dilepton.h" // MEMbbwwNtupleManager_dilepton
#include "hhAnalysis/bbwwMEMPerformanceStudies/interface/BJetTF_toy.h"                    // BJetTF_toy

#include <iostream> // std::cerr, std::fixed
#include <iomanip>  // std::setprecision(), std::setw()
#include <string>   // std::string
#include <vector>   // std::vector<>
#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <fstream>  // std::ofstream
#include <assert.h> // assert
#include <math.h>   // sqrt 

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
 * @brief Produce MEM Ntuple for HH signal and TT background events simulated with Delphes in HH->bbWW dilepton channel
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

  std::cout << "<produceMEMNtuple_hh_bb2l_delphes>:" << std::endl;

//--- keep track of time it takes the macro to execute
  TBenchmark clock;
  clock.Start("produceMEMNtuple_hh_bb2l_delphes");

//--- read python configuration parameters
  if ( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") )
    throw cms::Exception("produceMEMNtuple_hh_bb2l_delphes")
      << "No ParameterSet 'process' found in configuration file = " << argv[1] << " !!\n";

  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");

  edm::ParameterSet cfg_analyze = cfg.getParameter<edm::ParameterSet>("analyze_hh_bbwwMEM_dilepton");

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
  //std::string branchName_genMEt = cfg_analyze.getParameter<std::string>("branchName_genMEt");

  bool apply_genWeight = cfg_analyze.getParameter<bool>("apply_genWeight");

  bool isDEBUG = cfg_analyze.getParameter<bool>("isDEBUG");
  std::cout << " isDEBUG = " << isDEBUG << std::endl;

  fwlite::InputSource inputFiles(cfg);
  int maxEvents = inputFiles.maxEvents();
  std::cout << " maxEvents = " << maxEvents << std::endl;
  int skipSelEvents = cfg_analyze.getParameter<int>("skipSelEvents");
  std::cout << " skipSelEvents = " << skipSelEvents << std::endl;
  int maxSelEvents = cfg_analyze.getParameter<int>("maxSelEvents");
  std::cout << " maxSelEvents = " << maxSelEvents << std::endl;
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

  DelphesMEtReader* metReader = new DelphesMEtReader(branchName_met);
  inputTree->registerReader(metReader);

  DelphesGenParticleReader* genParticleReader = new DelphesGenParticleReader(branchName_genParticles);
  inputTree->registerReader(genParticleReader);

  DelphesGenJetReader* genJetReader = new DelphesGenJetReader(branchName_genJets);
  inputTree->registerReader(genJetReader);

  //DelphesMEtReader* genMEtReader = new DelphesMEtReader(branchName_genMEt);
  //inputTree->registerReader(genMEtReader);

  std::string ntupleDir = Form("%s/ntuples/%s", histogramDir.data(), process_string.data());
  MEMbbwwNtupleManager_dilepton* mem_ntuple = new MEMbbwwNtupleManager_dilepton(ntupleDir, "mem");
  mem_ntuple->makeTree(fs);
  mem_ntuple->initializeBranches();
  MEMbbwwNtupleManager_dilepton* mem_ntuple_missingBJet = new MEMbbwwNtupleManager_dilepton(ntupleDir, "mem_missingBJet");
  mem_ntuple_missingBJet->makeTree(fs);
  mem_ntuple_missingBJet->initializeBranches();

  mem::BJetTF_toy bjetTF;
  bjetTF.set_coeff(1.00); // CV: assume resolution on jet pT to be 100%*sqrt(pT)

  // random number generator for choosing b-jets
  TRandom3 rnd;
  rnd.SetSeed(12345);

  int analyzedEntries = 0;
  int skippedEntries  = 0;
  int selectedEntries = 0;
  double selectedEntries_weighted = 0.;
  TH1* histogram_analyzedEntries = fs.make<TH1D>("analyzedEntries", "analyzedEntries", 1, -0.5, +0.5);
  TH1* histogram_selectedEntries = fs.make<TH1D>("selectedEntries", "selectedEntries", 1, -0.5, +0.5);
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
    "m(ll) > 12 GeV",
  };
  CutFlowTableHistManager * cutFlowHistManager = new CutFlowTableHistManager(cutFlowTableCfg, cuts);
  cutFlowHistManager->bookHistograms(fs);
  while ( inputTree->hasNextEvent() && selectedEntries < maxSelEvents ) {
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

    //std::cout << "#leptons = " << (electrons.size() + muons.size()) << " (selected = " << selLeptons.size() << "):" << std::endl;
    //std::cout << " #electrons = " << electrons.size() << " (selected = " << selElectrons.size() << ")" << std::endl;
    //std::cout << " #muons = " << muons.size() << " (selected = " << selMuons.size() << ")" << std::endl;

//--- build collections of jets
    const std::vector<DelphesJet> jets = jetReader->read();
    const std::vector<const DelphesJet*> jet_ptrs = convert_to_ptrs(jets);
    const std::vector<const DelphesJet*> cleanedJets = jetCleaner_dR04(jet_ptrs, selLeptons)
    ;
    const std::vector<const DelphesJet*> selJets = jetSelector(cleanedJets, isBetterBJet);

    // CV: remove HH->bbZZ events (keep only HH->bbWW events)
    std::vector<DelphesGenParticle> genParticles = genParticleReader->read();
    const std::vector<const DelphesGenParticle*> genParticles_ptrs = convert_to_ptrs(genParticles);
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

    const DelphesJet* selJet1 = ( selJets.size() >= 1 ) ? selJets[0] : nullptr;
    const DelphesJet* selJet2 = ( selJets.size() >= 2 ) ? selJets[1] : nullptr;
    if ( !(selJet1 && selJet2) ) {
      if ( isDEBUG ) {
        std::cout << "event " << eventInfo.str() << " FAILS >= 2 jets selection\n";
      }
      continue;
    }
    cutFlowTable.update(">= 2 jets", evtWeight);
    cutFlowHistManager->fillHistograms(">= 2 jets", evtWeight);

    int numBJets = 0;
    if ( selJet1 && selJet1->btag() >= 1 ) ++numBJets;
    if ( selJet2 && selJet2->btag() >= 1 ) ++numBJets;
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

    bool selJet1_isFake = std::abs(selJet1->flavor()) != 5;
    bool selJet2_isFake = std::abs(selJet2->flavor()) != 5;

    int numGenBJets = 0;
    if ( !selJet1_isFake ) ++numGenBJets;
    if ( !selJet2_isFake ) ++numGenBJets;

    std::vector<DelphesGenJet> genJets = genJetReader->read();
    const std::vector<const DelphesGenJet*> genJets_ptrs = convert_to_ptrs(genJets);

    DelphesMEt met = metReader->read();
    TMatrixD metCov(2,2);
    metCov[0][0] = mem::square(25.); // CV: assume resolution on px, py components of MET to be 25 GeV
    metCov[1][0] = 0.;
    metCov[0][1] = 0.;
    metCov[1][1] = mem::square(25.);

    //DelphesMEt genMEt = genMEtReader->read();

    //---------------------------------------------------------------------------
    // CV: Skip running matrix element method (MEM) computation for the first 'skipSelEvents' events.
    //     This feature allows to process the HH signal samples in chunks of 'maxSelEvents' events per job.
    ++skippedEntries;
    if ( skippedEntries < skipSelEvents ) continue;
    //---------------------------------------------------------------------------

    //---------------------------------------------------------------------------
    // CV: Compute MEM likelihood ratio of HH signal and ttbar background hypotheses

    if ( isDEBUG ) {
      std::cout << "selLepton_lead: pT = " << selLepton_lead->pt() << "," 
                << " eta = " << selLepton_lead->eta() << ", phi = " << selLepton_lead->phi() << std::endl;
      std::cout << "selLepton_sublead: pT = " << selLepton_sublead->pt() << "," 
                << " eta = " << selLepton_sublead->eta() << ", phi = " << selLepton_sublead->phi() << std::endl;
      std::cout << "selJet1:";
      if ( selJet1 ) {
        std::cout << " pT = " << selJet1->pt() << ", eta = " << selJet1->eta() << ", phi = " << selJet1->phi() 
	          << " (isFake = " << selJet1_isFake << ")";
      } else {
        std::cout << " N/A";
      }
      std::cout << std::endl;
      std::cout << "selJet2:";
      if ( selJet2 ) {
        std::cout << " pT = " << selJet2->pt() << ", eta = " << selJet2->eta() << ", phi = " << selJet2->phi() 
	          << " (isFake = " << selJet2_isFake << ")";
      } else {
        std::cout << " N/A";
      }
      std::cout << std::endl;
    }

    int memLeptonType_lead;   
    double memLeptonMass_lead;   
    if ( selLepton_lead->is_electron() ) {
      memLeptonType_lead = mem::MeasuredParticle::kElectron;
      memLeptonMass_lead = mem::electronMass;
    } else if ( selLepton_lead->is_muon() ) {
      memLeptonType_lead = mem::MeasuredParticle::kMuon;
      memLeptonMass_lead = mem::muonMass;
    } else assert(0);
    int memLeptonType_sublead;   
    double memLeptonMass_sublead;   
    if ( selLepton_sublead->is_electron() ) {
      memLeptonType_sublead = mem::MeasuredParticle::kElectron;
      memLeptonMass_sublead = mem::electronMass;
    } else if ( selLepton_sublead->is_muon() ) {
      memLeptonType_sublead = mem::MeasuredParticle::kMuon;
      memLeptonMass_sublead = mem::muonMass;
    } else assert(0);

    mem::MeasuredParticle memMeasuredLepton_lead(memLeptonType_lead, 
      selLepton_lead->pt(), selLepton_lead->eta(), selLepton_lead->phi(), 
      memLeptonMass_lead, selLepton_lead->charge());
    mem::MeasuredParticle memMeasuredLepton_sublead(memLeptonType_sublead, 
      selLepton_sublead->pt(), selLepton_sublead->eta(), selLepton_sublead->phi(), 
      memLeptonMass_sublead, selLepton_sublead->charge());
    mem::MeasuredParticle memMeasuredBJet1(mem::MeasuredParticle::kBJet,
      selJet1->pt(), selJet1->eta(), selJet1->phi(), 
      mem::bottomQuarkMass);
    mem::MeasuredParticle memMeasuredBJet2(mem::MeasuredParticle::kBJet,
      selJet2->pt(), selJet2->eta(), selJet2->phi(), 
      mem::bottomQuarkMass);

    std::vector<mem::MeasuredParticle> memMeasuredParticles;
    memMeasuredParticles.push_back(memMeasuredLepton_lead);
    memMeasuredParticles.push_back(memMeasuredLepton_sublead);
    memMeasuredParticles.push_back(memMeasuredBJet1);
    memMeasuredParticles.push_back(memMeasuredBJet2);
    
    MEMEvent_dilepton memEvent(
      { (UInt_t)eventInfo.run(), (UInt_t)eventInfo.luminosityBlock(), (ULong64_t)eventInfo.event(), eventInfo.genweight() }, isSignal, 
      &memMeasuredBJet1, &memMeasuredBJet2, 
      &memMeasuredLepton_lead, &memMeasuredLepton_sublead,
      met.px(), met.py(), metCov);
    //addGenMatches_dilepton(memEvent, genBJetsForMatching_ptrs, genLeptonsForMatching_ptrs, genMEt.px(), genMEt.py());
    memEvent.set_numGenLeptons(2);
    memEvent.set_numGenBJets(numGenBJets);
    memEvent.set_numMeasuredBJets_loose(numBJets);
    memEvent.set_numMeasuredBJets_medium(numBJets);

    std::vector<mem::MeasuredParticle> memMeasuredParticles_missingBJet;
    memMeasuredParticles_missingBJet.push_back(memMeasuredLepton_lead);
    memMeasuredParticles_missingBJet.push_back(memMeasuredLepton_sublead);
    int selJetIdx_missingBJet = -1;
    enum { kFirst, kSecond }; 
    if        ( selJet1->btag() >= 1 && selJet2->btag() <  1 ) {
      selJetIdx_missingBJet = kFirst;
    } else if ( selJet1->btag() <  1 && selJet2->btag() >= 1 ) {
      selJetIdx_missingBJet = kSecond;
    } else {
      double u = rnd.Uniform();
      assert(u >= 0. && u <= 1.);
      if ( u > 0.50 ) selJetIdx_missingBJet = kFirst;
      else            selJetIdx_missingBJet = kSecond;
    }
    const mem::MeasuredParticle* memMeasuredBJet_missingBJet = nullptr;
    bool selJet_isFake_missingBJet;
    int numBJets_missingBJet = 0;
    if        ( selJetIdx_missingBJet == kFirst  ) {
      memMeasuredParticles_missingBJet.push_back(memMeasuredBJet1);
      memMeasuredBJet_missingBJet = &memMeasuredBJet1;
      selJet_isFake_missingBJet = selJet1_isFake;
    } else if ( selJetIdx_missingBJet == kSecond ) {
      memMeasuredParticles_missingBJet.push_back(memMeasuredBJet2);
      memMeasuredBJet_missingBJet = &memMeasuredBJet2;
      selJet_isFake_missingBJet = selJet2_isFake;
    } else assert(0);
    int numGenBJets_missingBJet = 0;
    if ( !selJet_isFake_missingBJet ) ++numGenBJets_missingBJet;

    MEMEvent_dilepton memEvent_missingBJet(
      { (UInt_t)eventInfo.run(), (UInt_t)eventInfo.luminosityBlock(), (ULong64_t)eventInfo.event(), eventInfo.genweight() }, isSignal, 
      memMeasuredBJet_missingBJet, nullptr,
      &memMeasuredLepton_lead, &memMeasuredLepton_sublead,
      met.px(), met.py(), metCov);
    //addGenMatches_dilepton(memEvent_missingBJet, genBJetsForMatching_ptrs, genLeptonsForMatching_ptrs, genMEt.px(), genMEt.py());
    memEvent_missingBJet.set_numGenLeptons(2);
    memEvent_missingBJet.set_numGenBJets(numGenBJets_missingBJet);
    memEvent_missingBJet.set_numMeasuredBJets_loose(numBJets_missingBJet);
    memEvent_missingBJet.set_numMeasuredBJets_medium(numBJets_missingBJet);

    const double sqrtS = 13.e+3;
    const std::string pdfName = "MSTW2008lo68cl";
    const std::string madgraphFileName_signal     = "hhAnalysis/bbwwMEM/data/param_hh_SM.dat";
    const std::string madgraphFileName_background = "hhAnalysis/bbwwMEM/data/param_ttbar.dat";
    const bool applyOnshellWmassConstraint_signal = false;
    const int memAlgo_verbosity = 0;
    //const int maxObjFunctionCalls_signal = 2500;
    //const int maxObjFunctionCalls_background = 25000;
    const int maxObjFunctionCalls_signal = 1000;
    const int maxObjFunctionCalls_background = 10000;

    clock.Reset();
    clock.Start("memAlgo");
    MEMbbwwAlgoDilepton memAlgo(sqrtS, pdfName, findFile(madgraphFileName_signal), findFile(madgraphFileName_background), memAlgo_verbosity);
    memAlgo.setBJet1TF(&bjetTF);
    memAlgo.setBJet2TF(&bjetTF);
    memAlgo.applyOnshellWmassConstraint_signal(applyOnshellWmassConstraint_signal);
    memAlgo.setIntMode(MEMbbwwAlgoDilepton::kVAMP);
    memAlgo.setMaxObjFunctionCalls_signal(maxObjFunctionCalls_signal);
    memAlgo.setMaxObjFunctionCalls_background(maxObjFunctionCalls_background);
    memAlgo.integrate(memMeasuredParticles, met.px(), met.py(), metCov);
    MEMbbwwResultDilepton memResult = memAlgo.getResult();
    clock.Stop("memAlgo");

    double memCpuTime = clock.GetCpuTime("memAlgo");
    if ( isDEBUG ) {
      std::cout << "MEM:"
	        << " probability for signal hypothesis = " << memResult.getProb_signal() 
                << " +/- " << memResult.getProbErr_signal() << ","
	        << " probability for background hypothesis = " << memResult.getProb_background() 
                << " +/- " << memResult.getProbErr_background() << " " 
	        << "--> likelihood ratio = " << memResult.getLikelihoodRatio() 
                << " +/- " << memResult.getLikelihoodRatioErr() 
	        << " (CPU time = " << memCpuTime << ")" << std::endl;
    }

    (const_cast<MEMEvent_dilepton*>(&memEvent))->set_memResult(memResult);
    (const_cast<MEMEvent_dilepton*>(&memEvent))->set_memCpuTime(memCpuTime);

    mem_ntuple->read(memEvent);
    mem_ntuple->fill();

    clock.Reset();
    clock.Start("memAlgo_missingBJet");
    MEMbbwwAlgoDilepton memAlgo_missingBJet(sqrtS, pdfName, findFile(madgraphFileName_signal), findFile(madgraphFileName_background), memAlgo_verbosity);
    memAlgo_missingBJet.setBJet1TF(&bjetTF);
    memAlgo_missingBJet.setBJet2TF(&bjetTF);
    memAlgo_missingBJet.applyOnshellWmassConstraint_signal(applyOnshellWmassConstraint_signal);
    memAlgo_missingBJet.setIntMode(MEMbbwwAlgoDilepton::kVAMP);
    memAlgo_missingBJet.setMaxObjFunctionCalls_signal(maxObjFunctionCalls_signal);
    memAlgo_missingBJet.setMaxObjFunctionCalls_background(maxObjFunctionCalls_background);
    memAlgo_missingBJet.integrate(memMeasuredParticles_missingBJet, met.px(), met.py(), metCov);
    MEMbbwwResultDilepton memResult_missingBJet = memAlgo_missingBJet.getResult();
    clock.Stop("memAlgo_missingBJet");
    
    double memCpuTime_missingBJet = clock.GetCpuTime("memAlgo_missingBJet");
    if ( isDEBUG ) {
      std::cout << "MEM (missing b-jet case):" 
	        << " probability for signal hypothesis = " << memResult_missingBJet.getProb_signal() 
                << " +/- " << memResult_missingBJet.getProbErr_signal() << ","
	        << " probability for background hypothesis = " << memResult_missingBJet.getProb_background() 
                << " +/- " << memResult_missingBJet.getProbErr_background() << " " 
	        << "--> likelihood ratio = " << memResult_missingBJet.getLikelihoodRatio() 
                << " +/- " << memResult_missingBJet.getLikelihoodRatioErr() 
	        << " (CPU time = " << memCpuTime_missingBJet << ")" << std::endl;
    }

    (const_cast<MEMEvent_dilepton*>(&memEvent_missingBJet))->set_memResult(memResult_missingBJet);
    (const_cast<MEMEvent_dilepton*>(&memEvent_missingBJet))->set_memCpuTime(memCpuTime_missingBJet);

    mem_ntuple_missingBJet->read(memEvent_missingBJet);
    mem_ntuple_missingBJet->fill();
    //---------------------------------------------------------------------------

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
  //delete genMEtReader;

  delete inputTree;

  clock.Show("analyze_hh_bb2l_delphes");

  return EXIT_SUCCESS;
}
