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
#include <TProfile2D.h> // TProfile2D

#include "hhAnalysis/DelphesAnalysis/interface/DelphesLepton.h"                   // DelphesLepton
#include "hhAnalysis/DelphesAnalysis/interface/DelphesLeptonReader.h"             // DelphesLeptonReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesLeptonCollectionSelector.h" // DelphesLeptonCollectionSelector
#include "hhAnalysis/DelphesAnalysis/interface/DelphesLeptonCollectionCleaner.h"  // DelphesLeptonCollectionCleaner
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJet.h"                      // DelphesJet
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetReader.h"                // DelphesJetReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetCollectionSelector.h"    // DelphesJetCollectionSelector
#include "hhAnalysis/DelphesAnalysis/interface/DelphesJetCollectionCleaner.h"     // DelphesJetCollectionCleaner
#include "hhAnalysis/DelphesAnalysis/interface/DelphesEventInfo.h"                // DelphesEventInfo
#include "hhAnalysis/DelphesAnalysis/interface/DelphesEventInfoReader.h"          // DelphesEventInfoReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenParticle.h"              // DelphesGenParticle
#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenParticleReader.h"        // DelphesGenParticleReader
#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenJet.h"                   // DelphesGenJet
#include "hhAnalysis/DelphesAnalysis/interface/DelphesGenJetReader.h"             // DelphesGenJetReader
#include "hhAnalysis/DelphesAnalysis/interface/delphesAuxFunctions.h"             // mergeLeptonCollections, isBetterBJet, genBQuarkFromHiggs_or_TopDecay

#include "tthAnalysis/HiggsToTauTau/interface/TTreeWrapper.h"                     // TTreeWrapper
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

void
fillProfile2d(TProfile2D* histogram, double x, double y, double z, double evtWeight)
{
  if ( !histogram ) return;
  const TAxis* xAxis = histogram->GetXaxis();
  double xMin = xAxis->GetXmin();
  double xMax = xAxis->GetXmax();
  const TAxis* yAxis = histogram->GetYaxis();
  double yMin = yAxis->GetXmin();
  double yMax = yAxis->GetXmax();
  //std::cout << "<fillProfile2d>:" << std::endl;
  //std::cout << " x = " << x << " [" << xMin << ".." << xMax << "]" << std::endl;
  //std::cout << " y = " << y << " [" << yMin << ".." << yMax << "]" << std::endl;
  if ( x > xMin && x < xMax && y > yMin && y < yMax ) {
    histogram->Fill(x, y, z, evtWeight);
  }
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

  edm::ParameterSet cfg_analyze = cfg.getParameter<edm::ParameterSet>("calibrate_jets_delphes");

  std::string treeName = cfg_analyze.getParameter<std::string>("treeName");

  std::string process_string = cfg_analyze.getParameter<std::string>("process");
  bool isSignal = ( process_string.find("signal") != std::string::npos ) ? true : false;

  std::string histogramDir = cfg_analyze.getParameter<std::string>("histogramDir");

  std::string branchName_electrons = cfg_analyze.getParameter<std::string>("branchName_electrons");
  std::string branchName_muons = cfg_analyze.getParameter<std::string>("branchName_muons");
  std::string branchName_jets = cfg_analyze.getParameter<std::string>("branchName_jets");

  std::string branchName_genParticles = cfg_analyze.getParameter<std::string>("branchName_genParticles");
  std::string branchName_genJets = cfg_analyze.getParameter<std::string>("branchName_genJets");

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
  jetSelector.getSelector().set_min_pt(20.);
  jetSelector.getSelector().set_max_absEta(5.0);

  DelphesGenParticleReader* genParticleReader = new DelphesGenParticleReader(branchName_genParticles);
  inputTree->registerReader(genParticleReader);

  DelphesGenJetReader* genJetReader = new DelphesGenJetReader(branchName_genJets);
  inputTree->registerReader(genJetReader);

//--- declare histograms used for jet energy calibration
  TFileDirectory dir = fs.mkdir(histogramDir);
  
  const int numPtBins = 13;
  double ptBinning[numPtBins + 1] = { 20., 25., 30., 35., 40., 50., 60., 80., 100., 120., 150., 200., 250., 1.e+3 };
  const int numEtaBins = 4;
  double etaBinning[numEtaBins + 1] = { 0.0, 1.4, 2.4, 3.0, 5.0 };
  TProfile2D* histogram_jetEnResponse      = dir.make<TProfile2D>("jetEnResponse",      "jetEnResponse",      numPtBins, ptBinning, numEtaBins, etaBinning);
  TProfile2D* histogram_jetEnResponse_b    = dir.make<TProfile2D>("jetEnResponse_b",    "jetEnResponse_b",    numPtBins, ptBinning, numEtaBins, etaBinning);
  TProfile2D* histogram_jetEnResponse_udsg = dir.make<TProfile2D>("jetEnResponse_udsg", "jetEnResponse_udsg", numPtBins, ptBinning, numEtaBins, etaBinning);
  TProfile2D* histogram_jetPtResponse      = dir.make<TProfile2D>("jetPtResponse",      "jetPtResponse",      numPtBins, ptBinning, numEtaBins, etaBinning);
  TProfile2D* histogram_jetPtResponse_b    = dir.make<TProfile2D>("jetPtResponse_b",    "jetPtResponse_b",    numPtBins, ptBinning, numEtaBins, etaBinning);
  TProfile2D* histogram_jetPtResponse_udsg = dir.make<TProfile2D>("jetPtResponse_udsg", "jetPtResponse_udsg", numPtBins, ptBinning, numEtaBins, etaBinning);
  
  TProfile2D* histogram_jetEnVisFrac       = dir.make<TProfile2D>("jetEnVisFrac",       "jetEnVisFrac",       numPtBins, ptBinning, numEtaBins, etaBinning);
  TProfile2D* histogram_jetPtVisFrac       = dir.make<TProfile2D>("jetPtVisFrac",       "jetPtVisFrac",       numPtBins, ptBinning, numEtaBins, etaBinning);

  int analyzedEntries = 0;
  int selectedEntries = 0;
  double selectedEntries_weighted = 0.;
  TH1* histogram_analyzedEntries = dir.make<TH1D>("analyzedEntries", "analyzedEntries", 1, -0.5, +0.5);
  TH1* histogram_selectedEntries = dir.make<TH1D>("selectedEntries", "selectedEntries", 1, -0.5, +0.5);
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

//--- read collections of electrons and muons
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

//--- read collections of jets
    const std::vector<DelphesJet> jets = jetReader->read();
    const std::vector<const DelphesJet*> jet_ptrs = convert_to_ptrs(jets);
    const std::vector<const DelphesJet*> cleanedJets = jetCleaner_dR04(jet_ptrs, selLeptons);
    const std::vector<const DelphesJet*> selJets = jetSelector(cleanedJets, isBetterBJet);

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

    if ( isDEBUG ) {
      std::cout << "#selJets = " << selJets.size() << std::endl;
      std::cout << "#genJets = " << genJets.size() << std::endl;
    }
    for ( std::vector<const DelphesJet*>::const_iterator recJet = selJets.begin();
          recJet != selJets.end(); ++recJet ) {
      if ( isDEBUG ) {
        std::cout << "processing recJet: pT = " << (*recJet)->pt() << ", eta = " << (*recJet)->eta() << ", phi = " << (*recJet)->phi() << std::endl;
      }
      const DelphesGenJet* genJet_bestMatch = get_genJet((*recJet)->p4(), genJets);
      if ( isDEBUG ) {
        std::cout << "genJet_bestMatch:";
        if ( genJet_bestMatch ) 
        {  
          std::cout << " pT = " << genJet_bestMatch->pt() << ", eta = " << genJet_bestMatch->eta() << ", phi = " << genJet_bestMatch->phi() << std::endl;
        }
        else
        {
          std::cout << " N/A";
        }
        std::cout << std::endl;
      }
      if ( genJet_bestMatch )
      {
        double genJetEn = genJet_bestMatch->p4().energy();
        double genJetPt = genJet_bestMatch->pt();
        double recJetEn = (*recJet)->p4().energy();
        double recJetPt = (*recJet)->pt();
        // CV: exclude badly measured jets from calibration
        if ( recJetPt > 0.5*genJetPt && recJetPt < 1.5*genJetPt ) {
          // CV: parametrize jet energy corrections by reconstructed jet pT/energy,
          //     to avoid biasing jet energy corrections by pT > 20 GeV cut applied on reconstructed jets in Delphes
          double refJetEn  = (*recJet)->p4().energy();
          double refJetPt  = (*recJet)->pt();
          double refJetEta = (*recJet)->absEta();
          fillProfile2d(histogram_jetEnResponse, refJetEn, refJetEta, recJetEn/genJetEn, evtWeight);
          fillProfile2d(histogram_jetPtResponse, refJetPt, refJetEta, recJetPt/genJetPt, evtWeight);

          bool isBJet = (*recJet)->btag() >= 1;
          TProfile2D* histogram_jetEnResponse_byFlavor = nullptr;
          TProfile2D* histogram_jetPtResponse_byFlavor = nullptr;
          if ( isBJet ) 
          {
            histogram_jetEnResponse_byFlavor = histogram_jetEnResponse_b;
            histogram_jetPtResponse_byFlavor = histogram_jetPtResponse_b;
          }
          else
          {
            histogram_jetEnResponse_byFlavor = histogram_jetEnResponse_udsg;
            histogram_jetPtResponse_byFlavor = histogram_jetPtResponse_udsg;
          }
          fillProfile2d(histogram_jetEnResponse_byFlavor, refJetEn, refJetEta, recJetEn/genJetEn, evtWeight);
          fillProfile2d(histogram_jetPtResponse_byFlavor, refJetPt, refJetEta, recJetPt/genJetPt, evtWeight);

          int mother_pdgId = ( isSignal ) ? 25 : 6; 
          const DelphesGenParticle* genBQuark = get_genBQuarkFromHiggs_or_TopDecay(genJet_bestMatch->p4(), genParticles, mother_pdgId);
          if ( genBQuark )
          {            
            double genBQuarkEn = genBQuark->p4().energy();
            double genBQuarkPt = genBQuark->pt();
            if ( genJetPt < 1.5*genBQuarkPt ) {
              // CV: parametrize correction for energy carried away by neutrinos in b-jets by corrected jet pT/energy,
              //     instead of by reconstructed jet pT/energy ??
              fillProfile2d(histogram_jetEnVisFrac, genJetEn, refJetEta, genJetEn/genBQuarkEn, evtWeight);
              fillProfile2d(histogram_jetPtVisFrac, genJetPt, refJetEta, genJetPt/genBQuarkPt, evtWeight);
            }
          }
        }
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
            << " selected = " << selectedEntries << " (weighted = " << selectedEntries_weighted << ")\n\n";
  std::cout << std::endl;

  delete muonReader;
  delete electronReader;
  delete jetReader;
  delete genParticleReader;
  delete genJetReader;

  delete inputTree;

  clock.Show("calibrate_jets_delphes");

  return EXIT_SUCCESS;
}
