
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <assert.h>

TH1* 
loadHistogram1d(TFile* inputFile, const std::string& histogramName)
{
  TH1* histogram = dynamic_cast<TH1*>(inputFile->Get(histogramName.data()));
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();
  return histogram;
}

TH2* 
loadHistogram2d(TFile* inputFile, const std::string& histogramName)
{
  TH1* histogram1d = loadHistogram1d(inputFile, histogramName);
  TH2* histogram2d = dynamic_cast<TH2*>(histogram1d);
  if ( !histogram2d ) {
    std::cerr << "Histogram = " << histogramName << " is not of type TH2 !!" << std::endl;
    assert(0);
  }
  return histogram2d;
}

void dump_hh_bb2l_delphes()
{
  gROOT->SetBatch(true);

  TH1::AddDirectory(false);

  std::string inputFilePath = "/home/veelken/CMSSW_11_1_2/CMSSW_11_1_2/src/hhAnalysis/DelphesAnalysis/test/";

  std::vector<std::string> inputFileNames;
  inputFileNames.push_back("analyze_hh_bb2l_delphes_signal.root");
  inputFileNames.push_back("analyze_hh_bb2l_delphes_background.root");

  std::vector<std::string> histogramNames_genMatch;
  histogramNames_genMatch.push_back("numGenBJets_geq1BJet");
  histogramNames_genMatch.push_back("numGenBJets_geq2BJets");

  std::vector<std::string> histogramNames_jetEnRes;
  histogramNames_jetEnRes.push_back("jetEnRes_b");
  histogramNames_jetEnRes.push_back("jetEnRes_udsg");

  std::vector<std::string> histogramNames_metRes;
  histogramNames_metRes.push_back("metResPx");
  histogramNames_metRes.push_back("metResPy");

  for ( std::vector<std::string>::const_iterator inputFileName = inputFileNames.begin();
        inputFileName != inputFileNames.end(); ++inputFileName ) {
    TString inputFileName_full = inputFilePath.data();
    if ( !inputFileName_full.EndsWith("/") ) inputFileName_full.Append("/");
    inputFileName_full.Append(inputFileName->data());
    TFile* inputFile = new TFile(inputFileName_full.Data());
    if ( !inputFile ) {
      std::cerr << "Failed to open input file = " << inputFileName_full.Data() << " !!" << std::endl;
      assert(0);
    }
    std::cout << "processing input file = '" << (*inputFileName) << "'..." << std::endl;

    for ( std::vector<std::string>::const_iterator histogramName_genMatch = histogramNames_genMatch.begin();
          histogramName_genMatch != histogramNames_genMatch.end(); ++histogramName_genMatch ) {
      TH1* histogram_genMatch = loadHistogram1d(inputFile, *histogramName_genMatch);
      std::cout << " " << (*histogramName_genMatch) << ": " 
                << histogram_genMatch->GetBinContent(1)/histogram_genMatch->Integral() << " / "
                << histogram_genMatch->GetBinContent(2)/histogram_genMatch->Integral() << " / "
                << histogram_genMatch->GetBinContent(3)/histogram_genMatch->Integral() << std::endl;
    }

    for ( std::vector<std::string>::const_iterator histogramName_jetEnRes = histogramNames_jetEnRes.begin();
          histogramName_jetEnRes != histogramNames_jetEnRes.end(); ++histogramName_jetEnRes ) {
      TH2* histogram_jetEnRes = loadHistogram2d(inputFile, *histogramName_jetEnRes);
      std::cout << " " << (*histogramName_jetEnRes) << ":" << std::endl;
      TAxis* xAxis = histogram_jetEnRes->GetXaxis();
      int numBinsX = xAxis->GetNbins();
      for ( int idxBinX = 1; idxBinX <= numBinsX; ++idxBinX ) 
      {
        double xMin = xAxis->GetBinLowEdge(idxBinX);
        double xMax = xAxis->GetBinUpEdge(idxBinX);
        std::string histogramName_proj = Form("%s_proj%i", histogram_jetEnRes->GetName(), idxBinX);
        TH1* histogram_proj = histogram_jetEnRes->ProjectionY(histogramName_proj.data(), idxBinX, idxBinX);
        std::cout << "  " << xMin << " < pT < " << xMax << ":"
                  << " mean = " << histogram_proj->GetMean() << " +/- " << histogram_proj->GetMeanError() << ","
                  << " rms = " << histogram_proj->GetRMS() << " +/- " << histogram_proj->GetRMSError() << std::endl;
        delete histogram_proj;
      }
    }

    for ( std::vector<std::string>::const_iterator histogramName_metRes = histogramNames_metRes.begin();
          histogramName_metRes != histogramNames_metRes.end(); ++histogramName_metRes ) {
      TH1* histogram_metRes = loadHistogram1d(inputFile, *histogramName_metRes);
      std::cout << " " << (*histogramName_metRes) << ":" 
                << " mean = " << histogram_metRes->GetMean() << " +/- " << histogram_metRes->GetMeanError() << ","
                << " rms = " << histogram_metRes->GetRMS() << " +/- " << histogram_metRes->GetRMSError() << std::endl;
    }

    std::cout << "...done." << std::endl;
    std::cout << std::endl;

    delete inputFile;
  }
}
