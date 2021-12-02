
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TF1.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TBenchmark.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <assert.h>

enum { kDisabled, kEnabled, kInverted }; 

enum { kUndefined, kSignal_lo, kBackground_lo };

enum { kProbS, kProbB, kLR, kMbb, kMll };

std::string getHistogramName(int idxHistogram)
{
  if      ( idxHistogram == kProbS ) return "probS";
  else if ( idxHistogram == kProbB ) return "probB";
  else if ( idxHistogram == kLR    ) return "memLR";
  else if ( idxHistogram == kMbb   ) return "mbb";
  else if ( idxHistogram == kMll   ) return "mll";
  else assert(0);
  return "";
}

bool makePlots_png  = true;
bool makePlots_pdf  = true;
bool makePlots_root = true;

double square(double x)
{
  return x*x;
}

//-------------------------------------------------------------------------------
// CV: the function getLikelihoodRatio has been copied from
//       hhAnalysis/bbwwMEM/interface/MEMResult.h
//     A monotonous transformation of the likelihood ratio is applied in case the "finebin" option is used,
//     to increase the "dynamic range" of the ROC curve.
double getLikelihoodRatio(double prob_signal, double prob_background, bool finebin = false)
{
  double memLR = 0.5;
  double prob_SplusB = prob_signal + prob_background;    
  if ( prob_SplusB > 0. ) {
    memLR = prob_signal/prob_SplusB;
  } else {
    memLR = 0.;
  }

  double retVal = memLR;
  if ( finebin ) {
    if      ( memLR < 0.5 ) retVal =  TMath::Log(memLR/0.5);
    else if ( memLR > 0.5 ) retVal = -TMath::Log((1.0 - memLR)/0.5);
    else                    retVal =  0.;
  }
  return retVal;
}
//-------------------------------------------------------------------------------

//-------------------------------------------------------------------------------
// CV: the functions openFile, loadTree, getHistogram, fillWithOverFlow, fillWithOverFlow_logx, fillHistograms
//     and the defintion of the struct histogramEntryType have been copied from
//       hhAnalysis/bbwwMEMPerformanceStudies/macros/debug_bbww_dilepton.C
//     The function fillHistograms has been modified in order to select the same events from the "regular" and the "missingBJet" trees.
//
TFile* openFile(const std::string& inputFilePath, const std::string& inputFileName)
{
  TString inputFileName_full = inputFilePath.data();
  if ( !inputFileName_full.EndsWith("/") ) inputFileName_full.Append("/");
  inputFileName_full.Append(inputFileName.data());
  TFile* inputFile = new TFile(inputFileName_full.Data());
  if ( !inputFile ) {
    std::cerr << "Failed to open input file = " << inputFileName_full.Data() << " !!" << std::endl;
    assert(0);
  }
  return inputFile;
}

TTree* loadTree(TFile* inputFile, const std::string& directory, const std::string& treeName)
{  
  TString treeName_full = directory.data();
  if ( !treeName_full.EndsWith("/") ) treeName_full.Append("/");
  treeName_full.Append(treeName.data());
  TTree* tree = (TTree*)inputFile->Get(treeName_full.Data());
  if ( !tree ) {
    std::cerr << "Failed to load tree = " << treeName_full.Data() << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  std::cout << "Successfully loaded tree = '" << treeName_full.Data() << "' from file = " << inputFile->GetName() << std::endl;
  std::cout << " #entries = " << tree->GetEntries() << std::endl;
  return tree;
}

struct histogramEntryType
{
  histogramEntryType()
  {
    histogram_memLR_         = new TH1D("memLR",         "memLR",           40,   0.,    1.);
    histogram_memLR_finebin_ = new TH1D("memLR_finebin", "memLR_finebin", 1000, -50.,  +50.);
    histogram_memProbS_      = new TH1D("memProbS",      "memProbS",        40, -70.,  +10.);
    histogram_memProbB_      = new TH1D("memProbB",      "memProbB",        40, -70.,  +10.);
    histogram_drbb_          = new TH1D("drbb",          "drbb",            40,   0.,    4.);
    histogram_mbb_           = new TH1D("mbb",           "mbb",           1000,   0., 1000.);
    histogram_drll_          = new TH1D("drll",          "drll",            40,   0.,    4.);
    histogram_dphill_        = new TH1D("dphill",        "dphill",          36,   0., TMath::Pi());
    histogram_mll_           = new TH1D("mll",           "mll",             40,   0.,  100.); 
  }
  ~histogramEntryType()
  {
    delete histogram_memLR_;
    delete histogram_memLR_finebin_;
    delete histogram_memProbS_;
    delete histogram_memProbB_;
    delete histogram_drbb_;
    delete histogram_mbb_;
    delete histogram_drll_;
    delete histogram_dphill_;
    delete histogram_mll_;
  }
  TH1* histogram_memLR_;
  TH1* histogram_memLR_finebin_;
  TH1* histogram_memProbS_;
  TH1* histogram_memProbB_;
  TH1* histogram_drbb_;
  TH1* histogram_mbb_;
  TH1* histogram_drll_;
  TH1* histogram_dphill_;
  TH1* histogram_mll_;
};

TH1* getHistogram(histogramEntryType* histograms, int idxHistogram, bool finebin = false)
{
  if      ( idxHistogram == kProbS ) return histograms->histogram_memProbS_;
  else if ( idxHistogram == kProbB ) return histograms->histogram_memProbB_;
  else if ( idxHistogram == kLR    ) return ( finebin ) ? histograms->histogram_memLR_finebin_ : histograms->histogram_memLR_;  
  else if ( idxHistogram == kMbb   ) return histograms->histogram_mbb_;
  else if ( idxHistogram == kMll   ) return histograms->histogram_mll_;
  else assert(0);
  return nullptr;
}

void fillWithOverFlow(TH1* histogram, double x, double evtWeight)
{
  if ( !histogram ) return;
  const TAxis* const xAxis = histogram->GetXaxis();
  int idxBin = xAxis->FindBin(x);
  if ( idxBin < 1                 ) idxBin = 1;
  if ( idxBin > xAxis->GetNbins() ) idxBin = xAxis->GetNbins();
  double binCenter = histogram->GetBinCenter(idxBin);
  histogram->Fill(binCenter, evtWeight);
}

void fillWithOverFlow_logx(TH1* histogram, double x, double evtWeight)
{    
  const double nonzero = 1.e-30;
  fillWithOverFlow(histogram, TMath::Log(TMath::Max(nonzero, x)), evtWeight);
}

void fillHistograms(histogramEntryType* histograms, 
                    TTree* tree, double sf_memProbS, double sf_memProbB,
                    histogramEntryType* histograms_missingBJet, 
                    TTree* tree_missingBJet, double sf_memProbS_missingBJet, double sf_memProbB_missingBJet,
                    const std::string& selection)
{
  assert(histograms && tree);

  int numEntries = tree->GetEntries();

  TTreeFormula* treeFormula = nullptr;
  if ( selection != "" ) {
    treeFormula = new TTreeFormula("treeFormula", selection.data(), tree);

    int numEntries_selected = 0;
    for ( int idxEntry = 0; idxEntry < numEntries; ++idxEntry ) {
      tree->GetEntry(idxEntry);
      if ( !treeFormula->EvalInstance() ) continue;
      ++numEntries_selected;
    }

    std::cout << "Applying selection: '" << selection << "'" << std::endl;
    std::cout << " " << numEntries_selected << " out of " << numEntries << " entries selected." << std::endl;
  }

  Double_t memProbS, memProbSerr, memProbB, memProbBerr;
  Double_t memLR    = -1.;
  Double_t memLRerr = -1.;
  tree->SetBranchAddress("memProbS",    &memProbS);
  tree->SetBranchAddress("memProbSerr", &memProbSerr);
  tree->SetBranchAddress("memProbB",    &memProbB);
  tree->SetBranchAddress("memProbBerr", &memProbBerr);
  Float_t drbb, mbb;
  tree->SetBranchAddress("drbb",        &drbb);
  tree->SetBranchAddress("mbb",         &mbb);
  Float_t drll, dphill, mll;
  tree->SetBranchAddress("drll",        &drll);
  tree->SetBranchAddress("dphill",      &dphill);
  tree->SetBranchAddress("mll",         &mll);

  Double_t memProbS_missingBJet, memProbSerr_missingBJet, memProbB_missingBJet, memProbBerr_missingBJet;
  Double_t memLR_missingBJet    = -1.;
  Double_t memLRerr_missingBJet = -1.;
  if ( histograms_missingBJet && tree_missingBJet ) {
    tree_missingBJet->SetBranchAddress("memProbS",    &memProbS_missingBJet);
    tree_missingBJet->SetBranchAddress("memProbSerr", &memProbSerr_missingBJet);
    tree_missingBJet->SetBranchAddress("memProbB",    &memProbB_missingBJet);
    tree_missingBJet->SetBranchAddress("memProbBerr", &memProbBerr_missingBJet);
  }

  for ( int idxEntry = 0; idxEntry < numEntries; ++idxEntry ) {
    tree->GetEntry(idxEntry);
    if ( treeFormula && !treeFormula->EvalInstance() ) continue;
    if ( histograms_missingBJet && tree_missingBJet ) {
      tree_missingBJet->GetEntry(idxEntry);
    }
    
    const double evtWeight = 1.;

    double memLR = getLikelihoodRatio(
      sf_memProbS*memProbS, 
      sf_memProbB*memProbB);
    double memLR_finebin = getLikelihoodRatio(
      sf_memProbS*memProbS, 
      sf_memProbB*memProbB, true);
    fillWithOverFlow(histograms->histogram_memLR_,                     memLR,                          evtWeight);
    fillWithOverFlow(histograms->histogram_memLR_finebin_,             memLR_finebin,                  evtWeight);
    fillWithOverFlow_logx(histograms->histogram_memProbS_,             memProbS,                       evtWeight);
    fillWithOverFlow_logx(histograms->histogram_memProbB_,             memProbB,                       evtWeight);
    fillWithOverFlow(histograms->histogram_drbb_,                      drbb,                           evtWeight);
    fillWithOverFlow(histograms->histogram_mbb_,                       mbb,                            evtWeight);
    fillWithOverFlow(histograms->histogram_drll_,                      drll,                           evtWeight);
    fillWithOverFlow(histograms->histogram_dphill_,                    TMath::Abs(dphill),             evtWeight);
    fillWithOverFlow(histograms->histogram_mll_,                       mll,                            evtWeight);

    if ( histograms_missingBJet && tree_missingBJet ) {
      double memLR_missingBJet = getLikelihoodRatio(
        sf_memProbS_missingBJet*memProbS_missingBJet, 
        sf_memProbB_missingBJet*memProbB_missingBJet);
      double memLR_missingBJet_finebin = getLikelihoodRatio(
        sf_memProbS_missingBJet*memProbS_missingBJet, 
        sf_memProbB_missingBJet*memProbB_missingBJet, true);
      fillWithOverFlow(histograms_missingBJet->histogram_memLR_,         memLR_missingBJet,              evtWeight);
      fillWithOverFlow(histograms_missingBJet->histogram_memLR_finebin_, memLR_missingBJet_finebin,      evtWeight);
      fillWithOverFlow_logx(histograms_missingBJet->histogram_memProbS_, memProbS_missingBJet,           evtWeight);
      fillWithOverFlow_logx(histograms_missingBJet->histogram_memProbB_, memProbB_missingBJet,           evtWeight);
    }
  }

  delete treeFormula;
}
//-------------------------------------------------------------------------------

TH1* addHistograms(const std::string& histogramSumName, const TH1* histogram1, const TH1* histogram2, const TH1* histogram3 = nullptr)
{
  TH1* histogramSum = (TH1*)histogram1->Clone(histogramSumName.data());
  histogramSum->Reset();
  if ( !histogramSum->GetSumw2N() ) histogramSum->Sumw2();
  histogramSum->Add(histogram1);
  histogramSum->Add(histogram2);
  if ( histogram3 ) histogramSum->Add(histogram3);
  return histogramSum;
}

TH1* compNormalizedHistogram(TH1* histogram)
{
  std::string histogramName_normalized = Form("%s_normalized", histogram->GetName());
  TH1* histogram_normalized = (TH1*)histogram->Clone(histogramName_normalized.data());
  if ( !histogram_normalized->GetSumw2N() ) histogram_normalized->Sumw2();
  double integral = histogram_normalized->Integral(1, histogram_normalized->GetNbinsX()); // CV: exclude underflow and overflow bins
  if ( integral > 0. ) histogram_normalized->Scale(1./integral);
  return histogram_normalized;
}

void showHistograms(double canvasSizeX, double canvasSizeY,
		    TH1* histogram1, const std::string& legendEntry1,
		    TH1* histogram2, const std::string& legendEntry2,
		    TH1* histogram3, const std::string& legendEntry3,
		    TH1* histogram4, const std::string& legendEntry4,
		    int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[], const std::vector<std::string>& drawOptions,
		    double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, const std::vector<std::string>& legendOptions, 
		    const std::string& labelText, double labelTextSize,
		    double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		    const std::string& xAxisTitle, double xAxisOffset,
		    bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		    const std::string& outputFileName)
{
  const double margin_top      = ( labelText != "" ) ? 0.065 : 0.025;
  const double margin_left     = ( outputFileName.find("_prob") != std::string::npos ) ? 0.162 : 0.150;
  const double margin_bottom   = ( outputFileName.find("_prob") != std::string::npos ) ? 0.140 : 0.125;
  const double margin_right    = 0.015;

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(margin_top);
  canvas->SetLeftMargin(margin_left);
  canvas->SetBottomMargin(margin_bottom);
  canvas->SetRightMargin(margin_right); 
  canvas->SetLogx(false);
  canvas->SetLogy(useLogScale);
  canvas->Draw();
  canvas->cd();
  
  assert(histogram1);
  TH1* histogram1_normalized = compNormalizedHistogram(histogram1);
  histogram1_normalized->SetFillColor(0);
  histogram1_normalized->SetFillStyle(0);
  histogram1_normalized->SetLineColor(colors[0]);
  histogram1_normalized->SetLineStyle(lineStyles[0]);
  histogram1_normalized->SetLineWidth(lineWidths[0]);
  histogram1_normalized->SetMarkerColor(colors[0]);
  histogram1_normalized->SetMarkerStyle(markerStyles[0]);
  histogram1_normalized->SetMarkerSize(markerSizes[0]);

  assert(histogram2);
  TH1* histogram2_normalized = compNormalizedHistogram(histogram2);
  histogram2_normalized->SetFillColor(0);
  histogram2_normalized->SetFillStyle(0);
  histogram2_normalized->SetLineColor(colors[1]);
  histogram2_normalized->SetLineStyle(lineStyles[1]);
  histogram2_normalized->SetLineWidth(lineWidths[1]);
  histogram2_normalized->SetMarkerColor(colors[1]);
  histogram2_normalized->SetMarkerStyle(markerStyles[1]);
  histogram2_normalized->SetMarkerSize(markerSizes[1]);

  TH1* histogram3_normalized = nullptr;
  if ( histogram3 ) {
    histogram3_normalized = compNormalizedHistogram(histogram3);
    histogram3_normalized->SetFillColor(0);
    histogram3_normalized->SetFillStyle(0);
    histogram3_normalized->SetLineColor(colors[2]);
    histogram3_normalized->SetLineStyle(lineStyles[2]);
    histogram3_normalized->SetLineWidth(lineWidths[2]);
    histogram3_normalized->SetMarkerColor(colors[2]);
    histogram3_normalized->SetMarkerStyle(markerStyles[2]);
    histogram3_normalized->SetMarkerSize(markerSizes[2]);
  }

  TH1* histogram4_normalized = nullptr;
  if ( histogram4 ) {
    histogram4_normalized = compNormalizedHistogram(histogram4);
    histogram4_normalized->SetFillColor(0);
    histogram4_normalized->SetFillStyle(0);
    histogram4_normalized->SetLineColor(colors[3]);
    histogram4_normalized->SetLineStyle(lineStyles[3]);
    histogram4_normalized->SetLineWidth(lineWidths[3]);
    histogram4_normalized->SetMarkerColor(colors[3]);
    histogram4_normalized->SetMarkerStyle(markerStyles[3]);
    histogram4_normalized->SetMarkerSize(markerSizes[3]);
  }
  
  TAxis* xAxis = histogram1_normalized->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);
  xAxis->SetTitleSize(60);
  xAxis->SetTitleFont(43);
  //xAxis->SetLabelOffset(-0.01);
  xAxis->SetLabelSize(0.050);
  xAxis->SetLabelFont(42);
  xAxis->SetTickLength(0.040);
  xAxis->SetNdivisions(505);

  TAxis* yAxis = histogram1_normalized->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);
  yAxis->SetTitleSize(60);
  yAxis->SetTitleFont(43);
  if ( yMax > yMin ) {
    yAxis->SetRangeUser(yMin, yMax);
  }
  //yAxis->SetLabelOffset(0.010);
  yAxis->SetLabelSize(0.055);
  yAxis->SetLabelFont(42);
  yAxis->SetTickLength(0.040);  
  yAxis->SetNdivisions(505);

  histogram1_normalized->SetTitle("");
  histogram1_normalized->SetStats(false);

  histogram1_normalized->Draw(Form("%ssame", drawOptions[0].data()));
  histogram2_normalized->Draw(Form("%ssame", drawOptions[1].data()));
  if ( histogram3_normalized ) histogram3_normalized->Draw(Form("%ssame", drawOptions[2].data()));
  if ( histogram4_normalized ) histogram4_normalized->Draw(Form("%ssame", drawOptions[3].data()));
  histogram1_normalized->Draw("axissame");

  TPaveText* label = nullptr;
  if ( labelText != "" ) {
    label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "NDC");
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->AddText(labelText.data());
    label->SetTextFont(42);
    label->SetTextSize(labelTextSize);
    label->SetTextColor(1);
    label->SetTextAlign(13);
    label->Draw();
  }

  TLegend* legend = nullptr;
  if ( legendEntry1 != "" && legendEntry2 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, NULL, "brNDC");
    legend->SetFillColor(10);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(legendTextSize);
    legend->SetTextColor(1);
    legend->SetMargin(0.20);
    legend->AddEntry(histogram1_normalized, legendEntry1.data(), legendOptions[0].data());
    legend->AddEntry(histogram2_normalized, legendEntry2.data(), legendOptions[1].data());
    if ( histogram3_normalized ) legend->AddEntry(histogram3_normalized, legendEntry3.data(), legendOptions[2].data());
    if ( histogram4_normalized ) legend->AddEntry(histogram4_normalized, legendEntry4.data(), legendOptions[3].data());
    legend->Draw();
  }

  canvas->Update();

  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  //if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  if ( makePlots_png  ) canvas->Print(std::string(outputFileName_plot).append(".png").data());
  if ( makePlots_pdf  ) canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  if ( makePlots_root ) canvas->Print(std::string(outputFileName_plot).append(".root").data());

  delete histogram1_normalized;
  delete histogram2_normalized;
  delete histogram3_normalized;
  delete histogram4_normalized;
  delete label;
  delete legend;
  delete canvas;
}

TH1* compRatioHistogram(const TH1* histogram_numerator, const TH1* histogram_denominator)
{
  assert(histogram_numerator->GetNbinsX() == histogram_denominator->GetNbinsX());

  std::string histogramRatioName = Form("%s_div_%s", histogram_numerator->GetName(), histogram_denominator->GetName());
  TH1* histogramRatio = (TH1*)histogram_numerator->Clone(histogramRatioName.data());
  histogramRatio->Reset();
  if ( !histogramRatio->GetSumw2N() ) histogramRatio->Sumw2();
  histogramRatio->SetTitle("");
  histogramRatio->SetStats(false);

  int numBinsX = histogram_numerator->GetNbinsX();
  for ( int idxBin = 1; idxBin <= numBinsX; ++idxBin ) {
    double binContent_numerator   = histogram_numerator->GetBinContent(idxBin);
    double binError_numerator     = histogram_numerator->GetBinError(idxBin);
    double binContent_denominator = histogram_denominator->GetBinContent(idxBin);
    double binError_denominator   = histogram_denominator->GetBinError(idxBin);

    if ( binContent_denominator > 0. ) {
      //double binContent_ratio = binContent_numerator/binContent_denominator - 1.;
      double binContent_ratio = binContent_numerator/binContent_denominator;
      double binErr2_ratio = 0.;
      binErr2_ratio += square(binError_numerator/binContent_denominator);
      binErr2_ratio += square(binError_denominator*binContent_numerator/square(binContent_denominator));
      double binErr_ratio = TMath::Sqrt(binErr2_ratio);
      histogramRatio->SetBinContent(idxBin, binContent_ratio);
      histogramRatio->SetBinError(idxBin, binErr_ratio);
    }
  }
  
  return histogramRatio;
}

void copyHistogramStyle(const TH1* histogram_source, TH1* histogram_target)
{
  histogram_target->SetMarkerColor(histogram_source->GetMarkerColor());
  histogram_target->SetMarkerSize(histogram_source->GetMarkerSize());
  histogram_target->SetMarkerStyle(histogram_source->GetMarkerStyle());
  histogram_target->SetLineColor(histogram_source->GetLineColor());
  histogram_target->SetLineWidth(histogram_source->GetLineWidth());
  histogram_target->SetLineStyle(histogram_source->GetLineStyle());
  histogram_target->SetFillColor(histogram_source->GetFillColor());
  histogram_target->SetFillStyle(histogram_source->GetFillStyle());
}

void showHistograms_wRatio(double canvasSizeX, double canvasSizeY,
			   TH1* histogramRef, const std::string& legendEntryRef,
			   TH1* histogram2, const std::string& legendEntry2,
			   TH1* histogram3, const std::string& legendEntry3,
			   TH1* histogram4, const std::string& legendEntry4,
			   int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[], const std::vector<std::string>& drawOptions,
			   double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, const std::vector<std::string>& legendOptions, 
			   const std::string& labelText, double labelTextSize,
			   double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
			   const std::string& xAxisTitle, double xAxisOffset,
			   bool useLogScale, double yMin, double yMax, double yMin_ratio, double yMax_ratio, const std::string& yAxisTitle, double yAxisOffset,
			   const std::string& outputFileName)
{
  const double margin_top      = ( labelText != "" ) ? 0.065 : 0.025;
  const double margin_left     = 0.170;
  const double margin_bottom   = 0.155;
  const double margin_right    = 0.015;
  
  const double topPad_sizeY    = 0.69;
  const double clearance_sizeY = 0.02;
  const double bottomPad_sizeY = 0.29;
  const double sum_sizeY       = topPad_sizeY + clearance_sizeY + bottomPad_sizeY;

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(0.);
  canvas->SetLeftMargin(0.);
  canvas->SetBottomMargin(0.);
  canvas->SetRightMargin(0.); 
  canvas->Draw();
  canvas->cd();

  const double topPad_x0 = 0.;
  const double topPad_x1 = 1.;
  const double topPad_y1 = 1.;
  const double topPad_y0 = topPad_y1 - topPad_sizeY/sum_sizeY;
  assert(topPad_x1 > topPad_x0 && topPad_y1 > topPad_y0);

  TPad* topPad = new TPad("topPad", "topPad", topPad_x0, topPad_y0, topPad_x1, topPad_y1);
  topPad->SetFillColor(10);
  topPad->SetBorderSize(0);
  topPad->SetTopMargin(margin_top);
  topPad->SetLeftMargin(margin_left);
  topPad->SetBottomMargin(0.000);
  topPad->SetRightMargin(margin_right); 
  topPad->SetLogx(false);
  topPad->SetLogy(useLogScale);
  topPad->Draw();
  topPad->cd();
  
  assert(histogramRef);
  TH1* histogramRef_normalized = compNormalizedHistogram(histogramRef);
  histogramRef_normalized->SetFillColor(0);
  histogramRef_normalized->SetFillStyle(0);
  histogramRef_normalized->SetLineColor(colors[0]);
  histogramRef_normalized->SetLineStyle(lineStyles[0]);
  histogramRef_normalized->SetLineWidth(lineWidths[0]);
  histogramRef_normalized->SetMarkerColor(colors[0]);
  histogramRef_normalized->SetMarkerStyle(markerStyles[0]);
  histogramRef_normalized->SetMarkerSize(markerSizes[0]);

  assert(histogram2);
  TH1* histogram2_normalized = compNormalizedHistogram(histogram2);
  histogram2_normalized->SetFillColor(0);
  histogram2_normalized->SetFillStyle(0);
  histogram2_normalized->SetLineColor(colors[1]);
  histogram2_normalized->SetLineStyle(lineStyles[1]);
  histogram2_normalized->SetLineWidth(lineWidths[1]);
  histogram2_normalized->SetMarkerColor(colors[1]);
  histogram2_normalized->SetMarkerStyle(markerStyles[1]);
  histogram2_normalized->SetMarkerSize(markerSizes[1]);

  TH1* histogram3_normalized = nullptr;
  if ( histogram3 ) {
    histogram3_normalized = compNormalizedHistogram(histogram3);
    histogram3_normalized->SetFillColor(0);
    histogram3_normalized->SetFillStyle(0);
    histogram3_normalized->SetLineColor(colors[2]);
    histogram3_normalized->SetLineStyle(lineStyles[2]);
    histogram3_normalized->SetLineWidth(lineWidths[2]);
    histogram3_normalized->SetMarkerColor(colors[2]);
    histogram3_normalized->SetMarkerStyle(markerStyles[2]);
    histogram3_normalized->SetMarkerSize(markerSizes[2]);
  }

  TH1* histogram4_normalized = nullptr;
  if ( histogram4 ) {
    histogram4_normalized = compNormalizedHistogram(histogram4);
    histogram4->SetFillColor(0);
    histogram4->SetFillStyle(0);
    histogram4->SetLineColor(colors[3]);
    histogram4->SetLineStyle(lineStyles[3]);
    histogram4->SetLineWidth(lineWidths[3]);
    histogram4->SetMarkerColor(colors[3]);
    histogram4->SetMarkerStyle(markerStyles[3]);
    histogram4->SetMarkerSize(markerSizes[3]);
  }
  
  TAxis* xAxis_top = histogramRef_normalized->GetXaxis();
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetTitleSize(65);
  xAxis_top->SetTitleFont(43);
  //xAxis_top->SetLabelOffset(-0.01);
  xAxis_top->SetLabelSize(0.050);
  xAxis_top->SetLabelFont(42);
  xAxis_top->SetTickLength(0.040);
  xAxis_top->SetNdivisions(505);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);

  TAxis* yAxis_top = histogramRef_normalized->GetYaxis();
  yAxis_top->SetTitle(yAxisTitle.data());
  yAxis_top->SetTitleOffset(yAxisOffset);
  yAxis_top->SetTitleSize(65);
  yAxis_top->SetTitleFont(43);
  if ( yMax > yMin ) {
    yAxis_top->SetRangeUser(yMin, yMax);
  }
  //yAxis_top->SetLabelOffset(0.010);
  yAxis_top->SetLabelSize(0.055);
  yAxis_top->SetLabelFont(42);
  yAxis_top->SetTickLength(0.040);  
  yAxis_top->SetNdivisions(505);

  histogramRef_normalized->SetTitle("");
  histogramRef_normalized->SetStats(false);

  histogramRef_normalized->Draw(drawOptions[0].data());
  histogram2_normalized->Draw(Form("%ssame", drawOptions[1].data()));
  if ( histogram3_normalized ) histogram3_normalized->Draw(Form("%ssame", drawOptions[2].data()));
  if ( histogram4_normalized ) histogram4_normalized->Draw(Form("%ssame", drawOptions[3].data()));
  histogramRef_normalized->Draw("axissame");

  TPaveText* label = nullptr;
  if ( labelText != "" ) {
    label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "NDC");
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->AddText(labelText.data());
    label->SetTextFont(42);
    label->SetTextSize(labelTextSize);
    label->SetTextColor(1);
    label->SetTextAlign(13);
    label->Draw();
  }

  TLegend* legend = nullptr;
  if ( legendEntryRef != "" && legendEntry2 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, NULL, "brNDC");
    legend->SetFillColor(10);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(legendTextSize);
    legend->SetTextColor(1);
    legend->SetMargin(0.20);
    legend->AddEntry(histogramRef_normalized, legendEntryRef.data(), legendOptions[0].data());
    legend->AddEntry(histogram2_normalized, legendEntry2.data(), legendOptions[1].data());
    if ( histogram3_normalized ) legend->AddEntry(histogram3_normalized, legendEntry3.data(), legendOptions[2].data());
    if ( histogram4_normalized ) legend->AddEntry(histogram4_normalized, legendEntry4.data(), legendOptions[3].data());
    legend->Draw();
  }

  canvas->cd();

  const double bottomPad_x0 = 0.;
  const double bottomPad_x1 = 1.;
  const double bottomPad_y0 = 0.;
  const double bottomPad_y1 = bottomPad_y0 + bottomPad_sizeY/sum_sizeY;
  assert(bottomPad_x1 > bottomPad_x0 && bottomPad_y1 > bottomPad_y0);

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", bottomPad_x0, bottomPad_y0, bottomPad_x1, bottomPad_y1);
  bottomPad->SetFillColor(10);
  bottomPad->SetBorderSize(0);
  bottomPad->SetTopMargin(0.000);
  bottomPad->SetLeftMargin(margin_left);
  bottomPad->SetBottomMargin(margin_bottom*topPad_sizeY/bottomPad_sizeY);
  bottomPad->SetRightMargin(margin_right);
  bottomPad->SetLogx(false);
  bottomPad->SetLogy(false);
  bottomPad->Draw();
  bottomPad->cd();

  TH1* histogramRatio2 = compRatioHistogram(histogram2_normalized, histogramRef_normalized);
  copyHistogramStyle(histogram2_normalized, histogramRatio2);

  TH1* histogramRatio3 = nullptr;
  if ( histogram3_normalized ) {
    histogramRatio3 = compRatioHistogram(histogram3_normalized, histogramRef_normalized);
    copyHistogramStyle(histogram3_normalized, histogramRatio3);
  }

  TH1* histogramRatio4 = nullptr;
  if ( histogram4_normalized ) {
    histogramRatio4 = compRatioHistogram(histogram4_normalized, histogramRef_normalized);
    copyHistogramStyle(histogram4_normalized, histogramRatio4);
  }

  TAxis* xAxis_bottom = histogramRatio2->GetXaxis();
  xAxis_bottom->SetTitle(xAxisTitle.data());
  xAxis_bottom->SetTitleOffset(xAxisOffset);
  xAxis_bottom->SetTitleSize(65);
  xAxis_bottom->SetTitleFont(43);
  //xAxis_bottom->SetLabelOffset(-0.01);
  xAxis_bottom->SetLabelSize(0.132);
  xAxis_bottom->SetLabelFont(42);
  xAxis_bottom->SetTickLength(0.105);
  xAxis_bottom->SetNdivisions(505);

  TAxis* yAxis_bottom = histogramRatio2->GetYaxis();
  yAxis_bottom->SetTitle("Ratio");
  yAxis_bottom->SetTitleOffset(yAxisOffset);
  yAxis_bottom->SetTitleSize(65);
  yAxis_bottom->SetTitleFont(43);
  //yAxis_bottom->SetLabelOffset(0.010);
  yAxis_bottom->SetLabelSize(0.132);
  yAxis_bottom->SetLabelFont(42);
  yAxis_bottom->SetTickLength(0.060);
  yAxis_bottom->SetNdivisions(505);

  histogramRatio2->SetMinimum(yMin_ratio);
  histogramRatio2->SetMaximum(yMax_ratio);

  histogramRatio2->Draw(drawOptions[1].data());
  
  TGraph* graph_line = new TGraph(2);
  graph_line->SetPoint(0, xAxis_bottom->GetXmin(), 1.);
  graph_line->SetPoint(1, xAxis_bottom->GetXmax(), 1.);
  graph_line->SetLineColor(colors[0]);
  graph_line->SetLineStyle(lineStyles[0]);
  graph_line->SetLineWidth(lineWidths[0]);
  graph_line->Draw("L");

  histogramRatio2->Draw("same");

  if ( histogramRatio3 ) histogramRatio3->Draw(Form("%ssame", drawOptions[2].data()));
  if ( histogramRatio4 ) histogramRatio4->Draw(Form("%ssame", drawOptions[3].data()));
  histogramRatio2->Draw("axissame");

  canvas->Update();

  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  //if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  if ( makePlots_png  ) canvas->Print(std::string(outputFileName_plot).append(".png").data());
  if ( makePlots_pdf  ) canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  if ( makePlots_root ) canvas->Print(std::string(outputFileName_plot).append(".root").data());

  delete histogramRef_normalized;
  delete histogram2_normalized;
  delete histogram3_normalized;
  delete histogram4_normalized;
  delete label;
  delete legend;
  delete histogramRatio2;
  delete histogramRatio3;
  delete histogramRatio4;
  delete graph_line;
  delete topPad;
  delete bottomPad;
  delete canvas;
}

TGraph* compGraphEfficiency(const std::string& graphName, const TH1* histogram)
{
  const TAxis* xAxis = histogram->GetXaxis();
  int numPoints = xAxis->GetNbins();
  TGraph* graphEfficiency = new TGraph(numPoints + 2);
  graphEfficiency->SetName(graphName.data());
  graphEfficiency->SetPoint(0, xAxis->GetXmin(), 1.0);
  double integral = histogram->Integral(1, histogram->GetNbinsX());
  double sum = 0.;
  for ( int idxPoint = 1; idxPoint <= numPoints; ++idxPoint ) {
    int idxBin = idxPoint;
    double binCenter = xAxis->GetBinCenter(idxBin);
    double binContent = histogram->GetBinContent(idxBin);
    sum += binContent;
    graphEfficiency->SetPoint(idxPoint, binCenter, 1.0 - (sum/integral));
  }
  graphEfficiency->SetPoint(numPoints + 1, xAxis->GetXmax(), 0.0);
  return graphEfficiency;
} 

TGraph* compGraphROC(const std::string& graphName, const TGraph* graphEfficiency_signal, const TGraph* graphEfficiency_background, bool useLogScale)
{
  //std::cout << "<compGraphROC>:" << std::endl;
  //std::cout << " graphName = " << graphName << std::endl;
  assert(graphEfficiency_signal->GetN() == graphEfficiency_background->GetN());
  int numPoints = graphEfficiency_signal->GetN();
  TGraph* graphROC = new TGraph(numPoints);
  graphROC->SetName(graphName.data());
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double x, efficiency_signal;
    graphEfficiency_signal->GetPoint(idxPoint, x, efficiency_signal);
    double efficiency_background = graphEfficiency_background->Eval(x);
    double xROC = efficiency_signal;
    double yROC;
    if ( useLogScale ) yROC = efficiency_background;
    else yROC = 1.0 - efficiency_background;
    //std::cout << "point #" << idxPoint << ": x = " << xROC << ", y = " << yROC << std::endl;
    graphROC->SetPoint(idxPoint, xROC, yROC);
  }
  return graphROC;
}

TGraph* compGraphROC(const std::string& graphName, const TH1* histogram_signal, const TH1* histogram_background, bool useLogScale)
{
  //std::cout << "<compGraphROC>:" << std::endl;
  //std::cout << " graphName = " << graphName << std::endl;
  assert(histogram_signal->GetNbinsX() == histogram_background->GetNbinsX());
  TGraph* graphEfficiency_signal = compGraphEfficiency(Form("%s_signal", graphName.data()), histogram_signal);
  TGraph* graphEfficiency_background = compGraphEfficiency(Form("%s_background", graphName.data()), histogram_background);
  TGraph* graphROC = compGraphROC(graphName, graphEfficiency_signal, graphEfficiency_background, useLogScale);
  delete graphEfficiency_signal;
  delete graphEfficiency_background;
  return graphROC;
}

double compS_over_SplusB(double binContent_signal, double binContent_background)
{
  double binContent_SplusB = binContent_signal + binContent_background;
  double S_over_SplusB = 0.;
  if   ( binContent_SplusB > 0. ) S_over_SplusB = binContent_signal/binContent_SplusB;
  else                            S_over_SplusB = 0.;
  return S_over_SplusB;
}

void normalize_efficiency(std::vector<double>& efficiency, double sum)
{
  assert(sum > 0.);
  int numPoints = efficiency.size();
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) 
  {
    efficiency[idxPoint] = 1.0 - efficiency[idxPoint]/sum;
  }
}

TGraph* compGraphROC_combined(const std::string& graphName, 
                              const TH1* histogram_2mediumBJets_signal, const TH1* histogram_2mediumBJets_background, 
                              const TH1* histogram_missingBJet_1mediumBJet_signal, const TH1* histogram_missingBJet_1mediumBJet_background, 
                              bool useLogScale)
{
  //std::cout << "<compGraphROC_combined>:" << std::endl;
  //std::cout << " graphName = " << graphName << std::endl;
  assert(histogram_2mediumBJets_signal->GetNbinsX() == histogram_2mediumBJets_background->GetNbinsX());
  int numBins_2mediumBJets = histogram_2mediumBJets_signal->GetNbinsX();
  assert(histogram_missingBJet_1mediumBJet_signal->GetNbinsX() == histogram_missingBJet_1mediumBJet_background->GetNbinsX());
  int numBins_missingBJet_1mediumBJet = histogram_missingBJet_1mediumBJet_signal->GetNbinsX();
  int idxBin_2mediumBJets = 1;
  int idxBin_missingBJet_1mediumBJet = 1;
  double sum_signal = 0.;
  double sum_background = 0.;
  std::vector<double> efficiency_signal;
  efficiency_signal.push_back(0.);
  std::vector<double> efficiency_background;
  efficiency_background.push_back(0.);
  bool isDone_2mediumBJets = false;
  bool isDone_missingBJet_1mediumBJet = false;
  while ( !(isDone_2mediumBJets && isDone_missingBJet_1mediumBJet) ) 
  {
    //if ( ((idxBin_2mediumBJets + idxBin_missingBJet_1mediumBJet) % 10) == 0 ) {
    //  std::cout << "idxBin: 2mediumBJets = " << idxBin_2mediumBJets << ", missingBJet_1mediumBJet = " << idxBin_missingBJet_1mediumBJet << std::endl;
    //}
    double binContent_2mediumBJets_signal = histogram_2mediumBJets_signal->GetBinContent(idxBin_2mediumBJets);
    double binContent_2mediumBJets_background = histogram_2mediumBJets_background->GetBinContent(idxBin_2mediumBJets);
    double binContent_missingBJet_1mediumBJet_signal = histogram_missingBJet_1mediumBJet_signal->GetBinContent(idxBin_missingBJet_1mediumBJet);
    double binContent_missingBJet_1mediumBJet_background = histogram_missingBJet_1mediumBJet_background->GetBinContent(idxBin_missingBJet_1mediumBJet);
    // CV: in principle, we should take the bin from the histogram with the lower S/B (i.e. the most background-like bin),
    //     but we prefer to take the bin with the lower S/(S+B), because the latter is numerically more robust.
    //     The result should be the same, as switching from S/B to S/(S+B) is a monotonous transformation.
    double S_over_B_2mediumBJets = compS_over_SplusB(binContent_2mediumBJets_signal, binContent_2mediumBJets_background);
    double S_over_B_missingBJet_1mediumBJet = compS_over_SplusB(binContent_missingBJet_1mediumBJet_signal, binContent_missingBJet_1mediumBJet_background);
    if ( (S_over_B_2mediumBJets <= S_over_B_missingBJet_1mediumBJet || isDone_missingBJet_1mediumBJet) && !isDone_2mediumBJets )
    {
      sum_signal += binContent_2mediumBJets_signal;
      sum_background += binContent_2mediumBJets_background;
      if ( idxBin_2mediumBJets < numBins_2mediumBJets ) ++idxBin_2mediumBJets;
      else isDone_2mediumBJets = true;
    } 
    else if ( !isDone_missingBJet_1mediumBJet )
    {
      sum_signal += binContent_missingBJet_1mediumBJet_signal;
      sum_background += binContent_missingBJet_1mediumBJet_background;
      if ( idxBin_missingBJet_1mediumBJet < numBins_missingBJet_1mediumBJet ) ++idxBin_missingBJet_1mediumBJet;
      else isDone_missingBJet_1mediumBJet = true;
    }
    else assert (0);
    efficiency_signal.push_back(sum_signal);
    efficiency_background.push_back(sum_background);
  }
  //std::cout << "sum_signal = " << sum_signal << ", sum_background = " << sum_background << std::endl;
  normalize_efficiency(efficiency_signal, sum_signal);
  normalize_efficiency(efficiency_background, sum_background);
  int numPoints_expected = numBins_2mediumBJets + numBins_missingBJet_1mediumBJet + 1;
  assert(efficiency_signal == numPoints_expected);
  assert(efficiency_background == numPoints_expected);
  int numPoints = efficiency_signal.size();
  TGraph* graphROC = new TGraph(numPoints);
  graphROC->SetName(graphName.data());
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) 
  {
    double xROC = efficiency_signal[idxPoint];
    double yROC;
    if ( useLogScale ) yROC = efficiency_background[idxPoint];
    else yROC = 1.0 - efficiency_background[idxPoint];
    //std::cout << "point #" << idxPoint << ": x = " << xROC << ", y = " << yROC << std::endl;
    graphROC->SetPoint(idxPoint, xROC, yROC);
  }
  return graphROC;
}

struct graphPoint
{
  graphPoint(double x, double y)
    : x_(x)
    , y_(y)
  {}
  ~graphPoint() {}
  double x_;
  double y_;
};

TGraph* sparsifyGraph(TGraph* graph, double minDeltaX = 0.025, double maxDeltaY = 0.100)
{
  std::vector<graphPoint> graphPoints_sparsified;
  double x_last = -1.e+3;
  double y_last = -1.e+3;
  int numPoints = graph->GetN();
  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double x, y;
    graph->GetPoint(idxPoint, x, y);
    if ( x < 0.01 ) continue; // CV: prevent point @ zero signal efficiency and zero background rate from being drawn
    //if ( x > 0.99 ) continue; // CV: prevent point @ 100% signal efficiency and 100% background rate from being drawn
    if ( idxPoint == 0 || TMath::Abs(x - x_last) > minDeltaX || TMath::Abs(y - y_last) > maxDeltaY || idxPoint == (numPoints - 1) ) {
      graphPoints_sparsified.push_back(graphPoint(x, y));
      x_last = x;
      y_last = y;
    }
  }
  int numPoints_sparsified = graphPoints_sparsified.size();
  TGraph* graph_sparsified = new TGraph(numPoints_sparsified);
  graph_sparsified->SetName(Form("%s_sparsified", graph->GetName()));
  for ( int idxPoint = 0; idxPoint < numPoints_sparsified; ++idxPoint ) {
    const graphPoint& graphPoint_sparsified = graphPoints_sparsified[idxPoint];
    graph_sparsified->SetPoint(idxPoint, graphPoint_sparsified.x_, graphPoint_sparsified.y_);
  }
  return graph_sparsified;
}

void setLineStyle(TLine* line)
{
  line->SetLineColor(14);
  line->SetLineStyle(7);
  line->SetLineWidth(1);
}

void showGraphs(double canvasSizeX, double canvasSizeY,
		TGraph* graph1, const std::string& legendEntry1,
		TGraph* graph2, const std::string& legendEntry2,
		TGraph* graph3, const std::string& legendEntry3,
		TGraph* graph4, const std::string& legendEntry4,
		int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[], const std::vector<std::string>& drawOptions,
		double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, const std::vector<std::string>& legendOptions,
		const std::string& labelText, double labelTextSize,
		double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		int numBinsX, double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
		bool useLogScale, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		const std::string& outputFileName)
{
  const double margin_top      = ( labelText != "" ) ? 0.065 : 0.025;
  const double margin_left     = 0.170;
  const double margin_bottom   = 0.145;
  const double margin_right    = 0.015;

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(margin_top);
  canvas->SetLeftMargin(margin_left);
  canvas->SetBottomMargin(margin_bottom);
  canvas->SetRightMargin(margin_right); 
  canvas->SetLogx(false);
  canvas->SetLogy(useLogScale);
  canvas->Draw();
  canvas->cd();

  assert(graph1);
  TGraph* graph1_sparsified = sparsifyGraph(graph1);
  graph1_sparsified->SetLineColor(colors[0]);
  graph1_sparsified->SetLineStyle(lineStyles[0]);
  graph1_sparsified->SetLineWidth(lineWidths[0]);
  graph1_sparsified->SetMarkerColor(colors[0]);
  graph1_sparsified->SetMarkerStyle(markerStyles[0]);
  graph1_sparsified->SetMarkerSize(markerSizes[0]);

  TGraph* graph2_sparsified = nullptr;
  if ( graph2 ) {    
    graph2_sparsified = sparsifyGraph(graph2);
    graph2_sparsified->SetLineColor(colors[1]);
    graph2_sparsified->SetLineStyle(lineStyles[1]);
    graph2_sparsified->SetLineWidth(lineWidths[1]);
    graph2_sparsified->SetMarkerColor(colors[1]);
    graph2_sparsified->SetMarkerStyle(markerStyles[1]);
    graph2_sparsified->SetMarkerSize(markerSizes[1]);
  }
  
  TGraph* graph3_sparsified = nullptr;
  if ( graph3 ) {    
    graph3_sparsified = sparsifyGraph(graph3);
    graph3_sparsified->SetLineColor(colors[2]);
    graph3_sparsified->SetLineStyle(lineStyles[2]);
    graph3_sparsified->SetLineWidth(lineWidths[2]);
    graph3_sparsified->SetMarkerColor(colors[2]);
    graph3_sparsified->SetMarkerStyle(markerStyles[2]);
    graph3_sparsified->SetMarkerSize(markerSizes[2]);
  }
  
  TGraph* graph4_sparsified = nullptr;
  if ( graph4 ) {
    graph4_sparsified = sparsifyGraph(graph4);
    graph4_sparsified->SetLineColor(colors[3]);
    graph4_sparsified->SetLineStyle(lineStyles[3]);
    graph4_sparsified->SetLineWidth(lineWidths[3]);
    graph4_sparsified->SetMarkerColor(colors[3]);
    graph4_sparsified->SetMarkerStyle(markerStyles[3]);
    graph4_sparsified->SetMarkerSize(markerSizes[3]);
  }

  TH1* dummyHistogram = new TH1D("dummyHistogram", "dummyHistogram", numBinsX, xMin, xMax);
  dummyHistogram->SetTitle("");
  dummyHistogram->SetStats(false);
  assert(yMax > yMin);
  dummyHistogram->SetMinimum(yMin);
  dummyHistogram->SetMaximum(yMax);

  TAxis* xAxis = dummyHistogram->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);
  xAxis->SetTitleSize(55);
  xAxis->SetTitleFont(43);
  //xAxis->SetLabelOffset(-0.01);
  xAxis->SetLabelSize(0.050);
  xAxis->SetLabelFont(42);
  xAxis->SetTickLength(0.040);
  xAxis->SetNdivisions(505);

  TAxis* yAxis = dummyHistogram->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);
  yAxis->SetTitleSize(55);
  yAxis->SetTitleFont(43);
  //yAxis->SetLabelOffset(0.010);
  yAxis->SetLabelSize(0.050);
  yAxis->SetLabelFont(42);
  yAxis->SetTickLength(0.040);  
  yAxis->SetNdivisions(505);
  
  dummyHistogram->Draw("axis");
  graph1_sparsified->Draw(drawOptions[0].data());
  if ( graph2_sparsified ) graph2_sparsified->Draw(drawOptions[1].data());
  if ( graph3_sparsified ) graph3_sparsified->Draw(drawOptions[2].data());
  if ( graph4_sparsified ) graph4_sparsified->Draw(drawOptions[3].data());

  std::vector<TLine*> lines;
  if ( outputFileName == "hh_bbwwMEM_dilepton_effectOfFakes_2graphs_ROC.pdf" ) {
    double x = 0.35;
    assert(graph1);
    double y1 = graph1->Eval(x);
    TLine* line1_horizontal = new TLine(0., y1, x, y1);
    setLineStyle(line1_horizontal);
    line1_horizontal->Draw();
    lines.push_back(line1_horizontal);
    assert(graph2);
    double y2 = graph2->Eval(x);
    TLine* line2_horizontal = new TLine(0., y2, x, y2);
    setLineStyle(line2_horizontal);
    line2_horizontal->Draw();
    lines.push_back(line2_horizontal);
    TLine* line_vertical = new TLine(x, yMin, x, TMath::Max(y1, y2));
    setLineStyle(line_vertical);
    line_vertical->Draw();
    lines.push_back(line_vertical);
  } else if ( outputFileName == "hh_bbwwMEM_dilepton_effectOfFakes_ROC_missingBJet.pdf" ) {
    double x = 0.35;
    assert(graph1);
    double y = graph1->Eval(x);
    TLine* line_horizontal = new TLine(0., y, x, y);
    setLineStyle(line_horizontal);
    line_horizontal->Draw();
    lines.push_back(line_horizontal);
    TLine* line_vertical = new TLine(x, yMin, x, y);
    setLineStyle(line_vertical);
    line_vertical->Draw();
    lines.push_back(line_vertical);
  }

  dummyHistogram->Draw("axissame");

  TPaveText* label = nullptr;
  if ( labelText != "" ) {
    label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "NDC");
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->AddText(labelText.data());
    label->SetTextFont(42);
    label->SetTextSize(labelTextSize);
    label->SetTextColor(1);
    label->SetTextAlign(13);
    label->Draw();
  }

  TLegend* legend = nullptr;
  if ( graph2_sparsified || graph3_sparsified || graph4_sparsified ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, NULL, "brNDC");
    legend->SetFillColor(10);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(legendTextSize);
    legend->SetTextColor(1);
    legend->SetMargin(0.20);
    legend->AddEntry(graph1_sparsified, legendEntry1.data(), legendOptions[0].data());
    if ( graph2_sparsified ) legend->AddEntry(graph2_sparsified, legendEntry2.data(), legendOptions[1].data());
    if ( graph3_sparsified ) legend->AddEntry(graph3_sparsified, legendEntry3.data(), legendOptions[2].data());
    if ( graph4_sparsified ) legend->AddEntry(graph4_sparsified, legendEntry4.data(), legendOptions[3].data());
    legend->Draw();
  }

  canvas->Update();

  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  //if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  if ( makePlots_png  ) canvas->Print(std::string(outputFileName_plot).append(".png").data());
  if ( makePlots_pdf  ) canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  if ( makePlots_root ) canvas->Print(std::string(outputFileName_plot).append(".root").data());

  delete label;
  delete legend;
  delete graph1_sparsified;
  delete graph2_sparsified;
  delete graph3_sparsified;
  delete graph4_sparsified;
  for ( std::vector<TLine*>::iterator line = lines.begin();
        line != lines.end(); ++line ) {
    delete (*line);
  }
  delete canvas;
}

TGraph* compRatioGraph(const TGraph* graph_numerator, const TGraph* graph_denominator)
{
  assert(graph_numerator->GetN() == graph_denominator->GetN());

  std::string graphRatioName = Form("%s_div_%s", graph_numerator->GetName(), graph_denominator->GetName());
  int numPoints = graph_numerator->GetN();
  TGraph* graphRatio = new TGraph(numPoints);

  for ( int idxPoint = 0; idxPoint < numPoints; ++idxPoint ) {
    double x, y_numerator;
    graph_numerator->GetPoint(idxPoint, x, y_numerator);

    double y_denominator = graph_denominator->Eval(x);

    if ( y_denominator > 0. ) {
      //double y_ratio = y_numerator/y_denominator - 1.;
      double y_ratio = y_numerator/y_denominator;
      graphRatio->SetPoint(idxPoint, x, y_ratio);
    }
  }
  
  return graphRatio;
}

void copyGraphStyle(const TGraph* graph_source, TGraph* graph_target)
{
  graph_target->SetMarkerColor(graph_source->GetMarkerColor());
  graph_target->SetMarkerSize(graph_source->GetMarkerSize());
  graph_target->SetMarkerStyle(graph_source->GetMarkerStyle());
  graph_target->SetLineColor(graph_source->GetLineColor());
  graph_target->SetLineWidth(graph_source->GetLineWidth());
  graph_target->SetLineStyle(graph_source->GetLineStyle());
}

void showGraphs_wRatio(double canvasSizeX, double canvasSizeY,
		       TGraph* graphRef, const std::string& legendEntryRef,
		       TGraph* graph2, const std::string& legendEntry2,
		       TGraph* graph3, const std::string& legendEntry3,
		       TGraph* graph4, const std::string& legendEntry4,
		       int colors[], int markerStyles[], int markerSizes[], int lineStyles[], int lineWidths[], const std::vector<std::string>& drawOptions,
		       double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY, const std::vector<std::string>& legendOptions, 
		       const std::string& labelText, double labelTextSize,
		       double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
		       int numBinsX, double xMin, double xMax, const std::string& xAxisTitle, double xAxisOffset,
		       bool useLogScale, double yMin, double yMax, double yMin_ratio, double yMax_ratio, const std::string& yAxisTitle, double yAxisOffset,
		       const std::string& outputFileName)
{
  const double margin_top      = ( labelText != "" ) ? 0.065 : 0.025;
  const double margin_left     = 0.170;
  const double margin_bottom   = 0.145;
  const double margin_right    = 0.015;
  
  const double topPad_sizeY    = 0.69;
  const double clearance_sizeY = 0.02;
  const double bottomPad_sizeY = 0.29;
  const double sum_sizeY       = topPad_sizeY + clearance_sizeY + bottomPad_sizeY;

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(0.);
  canvas->SetLeftMargin(0.);
  canvas->SetBottomMargin(0.);
  canvas->SetRightMargin(0.); 
  canvas->Draw();
  canvas->cd();

  const double topPad_x0 = 0.;
  const double topPad_x1 = 1.;
  const double topPad_y1 = 1.;
  const double topPad_y0 = topPad_y1 - topPad_sizeY/sum_sizeY;
  assert(topPad_x1 > topPad_x0 && topPad_y1 > topPad_y0);

  TPad* topPad = new TPad("topPad", "topPad", topPad_x0, topPad_y0, topPad_x1, topPad_y1);
  topPad->SetFillColor(10);
  topPad->SetBorderSize(0);
  topPad->SetTopMargin(margin_top);
  topPad->SetLeftMargin(margin_left);
  topPad->SetBottomMargin(0.000);
  topPad->SetRightMargin(margin_right); 
  topPad->SetLogx(false);
  topPad->SetLogy(useLogScale);
  topPad->Draw();
  topPad->cd();

  assert(graphRef);
  TGraph* graphRef_sparsified = sparsifyGraph(graphRef);
  graphRef_sparsified->SetLineColor(colors[0]);
  graphRef_sparsified->SetLineStyle(lineStyles[0]);
  graphRef_sparsified->SetLineWidth(lineWidths[0]);
  graphRef_sparsified->SetMarkerColor(colors[0]);
  graphRef_sparsified->SetMarkerStyle(markerStyles[0]);
  graphRef_sparsified->SetMarkerSize(markerSizes[0]);

  assert(graph2);
  TGraph* graph2_sparsified = sparsifyGraph(graph2);
  graph2_sparsified->SetLineColor(colors[1]);
  graph2_sparsified->SetLineStyle(lineStyles[1]);
  graph2_sparsified->SetLineWidth(lineWidths[1]);
  graph2_sparsified->SetMarkerColor(colors[1]);
  graph2_sparsified->SetMarkerStyle(markerStyles[1]);
  graph2_sparsified->SetMarkerSize(markerSizes[1]);
  
  TGraph* graph3_sparsified = nullptr;
  if ( graph3 ) {    
    graph3_sparsified = sparsifyGraph(graph3);
    graph3_sparsified->SetLineColor(colors[2]);
    graph3_sparsified->SetLineStyle(lineStyles[2]);
    graph3_sparsified->SetLineWidth(lineWidths[2]);
    graph3_sparsified->SetMarkerColor(colors[2]);
    graph3_sparsified->SetMarkerStyle(markerStyles[2]);
    graph3_sparsified->SetMarkerSize(markerSizes[2]);
  }
  
  TGraph* graph4_sparsified = nullptr;
  if ( graph4 ) {
    graph4_sparsified = sparsifyGraph(graph4);
    graph4_sparsified->SetLineColor(colors[3]);
    graph4_sparsified->SetLineStyle(lineStyles[3]);
    graph4_sparsified->SetLineWidth(lineWidths[3]);
    graph4_sparsified->SetMarkerColor(colors[3]);
    graph4_sparsified->SetMarkerStyle(markerStyles[3]);
    graph4_sparsified->SetMarkerSize(markerSizes[3]);
  }

  TH1* dummyHistogram_top = new TH1D("dummyHistogram_top", "dummyHistogram_top", numBinsX, xMin, xMax);
  dummyHistogram_top->SetTitle("");
  dummyHistogram_top->SetStats(false);
  assert(yMax > yMin);
  dummyHistogram_top->SetMinimum(yMin);
  dummyHistogram_top->SetMaximum(yMax);

  TAxis* xAxis_top = dummyHistogram_top->GetXaxis();
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetTitle(xAxisTitle.data());
  xAxis_top->SetTitleOffset(xAxisOffset);
  xAxis_top->SetTitleSize(60);
  xAxis_top->SetTitleFont(43);
  //xAxis_top->SetLabelOffset(-0.01);
  xAxis_top->SetLabelSize(0.050);
  xAxis_top->SetLabelFont(42);
  xAxis_top->SetTickLength(0.040);
  xAxis_top->SetNdivisions(505);
  xAxis_top->SetLabelColor(10);
  xAxis_top->SetTitleColor(10);

  TAxis* yAxis_top = dummyHistogram_top->GetYaxis();
  yAxis_top->SetTitle(yAxisTitle.data());
  yAxis_top->SetTitleOffset(yAxisOffset);
  yAxis_top->SetTitleSize(60);
  yAxis_top->SetTitleFont(43);
  //yAxis_top->SetLabelOffset(0.010);
  yAxis_top->SetLabelSize(0.055);
  yAxis_top->SetLabelFont(42);
  yAxis_top->SetTickLength(0.040);  
  yAxis_top->SetNdivisions(505);
  
  dummyHistogram_top->Draw("axis");
  graphRef_sparsified->Draw(drawOptions[0].data());
  graph2_sparsified->Draw(drawOptions[1].data());
  if ( graph3_sparsified ) graph3_sparsified->Draw(drawOptions[2].data());
  if ( graph4_sparsified ) graph4_sparsified->Draw(drawOptions[3].data());
  dummyHistogram_top->Draw("axissame");

  TPaveText* label = nullptr;
  if ( labelText != "" ) {
    label = new TPaveText(labelPosX, labelPosY, labelPosX + labelSizeX, labelPosY + labelSizeY, "NDC");
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->AddText(labelText.data());
    label->SetTextFont(42);
    label->SetTextSize(labelTextSize);
    label->SetTextColor(1);
    label->SetTextAlign(13);
    label->Draw();
  }

  TLegend* legend = nullptr;
  if ( legendEntryRef != "" && legendEntry2 != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, NULL, "brNDC");
    legend->SetFillColor(10);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(legendTextSize);
    legend->SetTextColor(1);
    legend->SetMargin(0.20);
    legend->AddEntry(graphRef_sparsified, legendEntryRef.data(), legendOptions[0].data());
    legend->AddEntry(graph2_sparsified, legendEntry2.data(), legendOptions[1].data());
    if ( graph3_sparsified ) legend->AddEntry(graph3_sparsified, legendEntry3.data(), legendOptions[2].data());
    if ( graph4_sparsified ) legend->AddEntry(graph4_sparsified, legendEntry4.data(), legendOptions[3].data());
    legend->Draw();
  }

  canvas->cd();

  const double bottomPad_x0 = 0.;
  const double bottomPad_x1 = 1.;
  const double bottomPad_y0 = 0.;
  const double bottomPad_y1 = bottomPad_y0 + bottomPad_sizeY/sum_sizeY;
  assert(bottomPad_x1 > bottomPad_x0 && bottomPad_y1 > bottomPad_y0);

  TPad* bottomPad = new TPad("bottomPad", "bottomPad", bottomPad_x0, bottomPad_y0, bottomPad_x1, bottomPad_y1);
  bottomPad->SetFillColor(10);
  bottomPad->SetBorderSize(0);
  bottomPad->SetTopMargin(0.000);
  bottomPad->SetLeftMargin(margin_left);
  bottomPad->SetBottomMargin(margin_bottom*topPad_sizeY/bottomPad_sizeY);
  bottomPad->SetRightMargin(margin_right);
  bottomPad->SetLogx(false);
  bottomPad->SetLogy(false);
  bottomPad->Draw();
  bottomPad->cd();

  TH1* dummyHistogram_bottom = new TH1D("dummyHistogram_bottom", "dummyHistogram_bottom", numBinsX, xMin, xMax);
  dummyHistogram_bottom->SetTitle("");
  dummyHistogram_bottom->SetStats(false);

  TGraph* graphRatio2_sparsified = compRatioGraph(graph2_sparsified, graphRef_sparsified);
  copyGraphStyle(graph2_sparsified, graphRatio2_sparsified);

  TGraph* graphRatio3_sparsified = nullptr;
  if ( graph3_sparsified ) {
    graphRatio3_sparsified = compRatioGraph(graph3_sparsified, graphRef_sparsified);
    copyGraphStyle(graph3_sparsified, graphRatio3_sparsified);
  }

  TGraph* graphRatio4_sparsified = nullptr;
  if ( graph4_sparsified ) {
    graphRatio4_sparsified = compRatioGraph(graph4_sparsified, graphRef_sparsified);
    copyGraphStyle(graph4_sparsified, graphRatio4_sparsified);
  }

  TAxis* xAxis_bottom = dummyHistogram_bottom->GetXaxis();
  xAxis_bottom->SetTitle(xAxisTitle.data());
  xAxis_bottom->SetTitleOffset(xAxisOffset);
  xAxis_bottom->SetTitleSize(60);
  xAxis_bottom->SetTitleFont(43);
  //xAxis_bottom->SetLabelOffset(-0.01);
  xAxis_bottom->SetLabelSize(0.132);
  xAxis_bottom->SetLabelFont(42);
  xAxis_bottom->SetTickLength(0.105);
  xAxis_bottom->SetNdivisions(505);

  TAxis* yAxis_bottom = dummyHistogram_bottom->GetYaxis();
  yAxis_bottom->SetTitle("Ratio");
  yAxis_bottom->SetTitleOffset(yAxisOffset);
  yAxis_bottom->SetTitleSize(60);
  yAxis_bottom->SetTitleFont(43);
  //yAxis_bottom->SetLabelOffset(0.010);
  yAxis_bottom->SetLabelSize(0.132);
  yAxis_bottom->SetLabelFont(42);
  yAxis_bottom->SetTickLength(0.060);
  yAxis_bottom->SetNdivisions(505);

  dummyHistogram_bottom->SetMinimum(yMin_ratio);
  dummyHistogram_bottom->SetMaximum(yMax_ratio);

  dummyHistogram_bottom->Draw("axis");

  TGraph* graph_line = new TGraph(2);
  graph_line->SetPoint(0, xAxis_bottom->GetXmin(), 1.);
  graph_line->SetPoint(1, xAxis_bottom->GetXmax(), 1.);
  graph_line->SetLineColor(colors[0]);
  graph_line->SetLineStyle(lineStyles[0]);
  graph_line->SetLineWidth(lineWidths[0]);
  graph_line->Draw("L");

  graphRatio2_sparsified->Draw(drawOptions[1].data());
  if ( graphRatio3_sparsified ) graphRatio3_sparsified->Draw(drawOptions[2].data());
  if ( graphRatio4_sparsified ) graphRatio4_sparsified->Draw(drawOptions[3].data());
  dummyHistogram_bottom->Draw("axissame");

  canvas->Update();

  std::string outputFileName_plot = "plots/";
  size_t idx = outputFileName.find_last_of('.');
  outputFileName_plot.append(std::string(outputFileName, 0, idx));
  //if ( idx != std::string::npos ) canvas->Print(std::string(outputFileName_plot).append(std::string(outputFileName, idx)).data());
  if ( makePlots_png  ) canvas->Print(std::string(outputFileName_plot).append(".png").data());
  if ( makePlots_pdf  ) canvas->Print(std::string(outputFileName_plot).append(".pdf").data());
  if ( makePlots_root ) canvas->Print(std::string(outputFileName_plot).append(".root").data());

  delete label;
  delete legend;
  delete dummyHistogram_top;
  delete dummyHistogram_bottom;
  delete graphRef_sparsified;
  delete graph2_sparsified;
  delete graphRatio2_sparsified;
  delete graph3_sparsified;
  delete graphRatio3_sparsified;
  delete graph4_sparsified;
  delete graphRatio4_sparsified;
  delete graph_line;
  delete canvas;
}

void add_to_selection(std::string& selection, const std::string& to_add)
{
  if ( selection != "" ) selection += " && ";
  selection += to_add;
}

void comp_mbb_window(TH1* histogram, double& xMin_window, double& xMax_window)
{
  double xMin_fit =  50.;
  double xMax_fit = 200.;

  TAxis* xAxis = histogram->GetXaxis();
  xAxis->SetRangeUser(xMin_fit, xMax_fit);

  double mean = histogram->GetMean();
  double rms = histogram->GetRMS();

  xAxis->SetRange(); // reset x-axis range

  std::string fitFunctionFormula = "[0] + [1]*x + [2]*TMath::Gaus(x, [3], [4])";
  std::string fitFunctionName = "fitFunction";
  TF1* fitFunction = new TF1(fitFunctionName.data(), fitFunctionFormula.data(), xMin_fit, xMax_fit);
  fitFunction->SetParameter(0, 0.);
  fitFunction->SetParameter(1, 0.);
  fitFunction->SetParameter(2, histogram->Integral());
  fitFunction->SetParameter(3, mean);
  fitFunction->SetParameter(4, rms);
  histogram->Fit(fitFunction);

  std::cout << "fitted Gaussian: mean = " << fitFunction->GetParameter(3) << ", rms = " << fitFunction->GetParameter(4) << std::endl;
 
  double numSigma = 2.;
  xMin_window = mean - numSigma*rms;
  xMax_window = mean + numSigma*rms;

  delete fitFunction;  
}

void makeMEMPerformancePlotsFromNtuples_bbww_dilepton_delphes()
{
  gROOT->SetBatch(true);

  TH1::AddDirectory(false);
 
  bool makePlots_signal_vs_background = true;
  bool makePlots_effectOfFakes_by_numBJets = true;
  bool makePlots_effectOfFakes_by_mbb = true;

  std::string inputFilePath = "/hdfs/local/veelken/hhAnalysis/2016/2021Nov29_wJetCalibration/histograms/hh_bbwwMEM_dilepton_delphes/";
  std::string inputFileName = "histograms_harvested_stage2.root";
  TFile* inputFile = openFile(inputFilePath, inputFileName);

  std::string directory = "hh_bbwwMEM_dilepton/ntuples";
  std::string treeName = "mem";
  std::string treeName_missingBJet = "mem_missingBJet";

  std::map<int, std::string> subdirectories; // key = idxProcess
  subdirectories[kSignal_lo] = "signal_lo";
  subdirectories[kBackground_lo] = "background_lo";

  typedef std::map<int, TH1*>          histogramMap1;
  typedef std::map<int, histogramMap1> histogramMap2;
  typedef std::map<int, histogramMap2> histogramMap3;
  typedef std::map<int, histogramMap3> histogramMap4;
  typedef std::map<int, histogramMap4> histogramMap5;
  histogramMap5 histograms;                           // keys = numBJets_medium, numGenBJets, mbbCut, idxProcess, idxHistogram
  histogramMap4 histograms_memLR_finebin;             // keys = numBJets_medium, numGenBJets, mbbCut, idxProcess
  histogramMap5 histograms_missingBJet;               // keys = numBJets_medium, numGenBJets, mbbCut, idxProcess, idxHistogram
  histogramMap4 histograms_missingBJet_memLR_finebin; // keys = numBJets_medium, numGenBJets, mbbCut, idxProcess

  TBenchmark clock;
  clock.Start("makeMEMPerformancePlotsFromNtuples_bbww_dilepton_delphes");

  double mbb_min = -1.;
  double mbb_max = -1.;
  histogramEntryType* tmpHistograms_mbb_signal = new histogramEntryType();
  TString directory_mbb_signal = directory.data();
  if ( !directory_mbb_signal.EndsWith("/") ) directory_mbb_signal.Append("/");
  directory_mbb_signal.Append(subdirectories[kSignal_lo]);
  TTree* tree_mbb_signal = loadTree(inputFile, directory_mbb_signal.Data(), treeName);
  fillHistograms(
    tmpHistograms_mbb_signal, 
    tree_mbb_signal, 1., 1.,
    nullptr,
    nullptr, 1., 1.,
    "gen_nbjets == 2");
  TH1* histogram_mbb_signal = getHistogram(tmpHistograms_mbb_signal, kMbb);
  comp_mbb_window(histogram_mbb_signal, mbb_min, mbb_max);

  for ( int numBJets_medium = -1; numBJets_medium <= 2; ++numBJets_medium ) {
    for ( int numGenBJets = -1; numGenBJets <= 2; ++numGenBJets ) {
      for ( int mbbCut = kDisabled; mbbCut <= kInverted; ++mbbCut ) {
        for ( int idxProcess = kSignal_lo; idxProcess <= kBackground_lo; ++idxProcess ) {

          histogramEntryType* tmpHistograms = new histogramEntryType();
  
          histogramEntryType* tmpHistograms_missingBJet = new histogramEntryType();

          TString directory_full = directory.data();
          if ( !directory_full.EndsWith("/") ) directory_full.Append("/");
          directory_full.Append(subdirectories[idxProcess]);

          TTree* tree             = loadTree(inputFile, directory_full.Data(), treeName);
          TTree* tree_missingBJet = loadTree(inputFile, directory_full.Data(), treeName_missingBJet);
          
          std::string selection;
          if      ( numBJets_medium == 2         ) add_to_selection(selection, "nbjets_medium == 2");
          else if ( numBJets_medium == 1         ) add_to_selection(selection, "nbjets_medium == 1");
          else if ( numBJets_medium == 0         ) add_to_selection(selection, "nbjets_medium == 0");
          if      ( numGenBJets     == 2         ) add_to_selection(selection, "gen_nbjets == 2");
          else if ( numGenBJets     == 1         ) add_to_selection(selection, "gen_nbjets == 1");
          else if ( numGenBJets     == 0         ) add_to_selection(selection, "gen_nbjets == 0");
          if      ( mbbCut          == kEnabled  ) add_to_selection(selection, Form("(mbb > %1.1f && mbb < %1.1f)", mbb_min, mbb_max));
          if      ( mbbCut          == kInverted ) add_to_selection(selection, Form("(mbb < %1.1f || mbb > %1.1f)", mbb_min, mbb_max));
          if ( selection != "" ) {
            double sf_memProbS = 1.e+5;
            //double sf_memProbS = 1.;
            double sf_memProbB = 1.;

            double sf_memProbS_missingBJet = 1.e+5;
            //double sf_memProbS_missingBJet = 1.;
            double sf_memProbB_missingBJet = 1.;

            fillHistograms(
              tmpHistograms, 
              tree, sf_memProbS, sf_memProbB,
              tmpHistograms_missingBJet,
              tree_missingBJet, sf_memProbS_missingBJet, sf_memProbB_missingBJet,
              selection);
          }
 
          for ( int idxHistogram = kProbS; idxHistogram <= kMll; ++idxHistogram ) {
            histograms[numBJets_medium][numGenBJets][mbbCut][idxProcess][idxHistogram] = getHistogram(tmpHistograms, idxHistogram);
            histograms_missingBJet[numBJets_medium][numGenBJets][mbbCut][idxProcess][idxHistogram] = getHistogram(tmpHistograms_missingBJet, idxHistogram);
          }
          histograms_memLR_finebin[numBJets_medium][numGenBJets][mbbCut][idxProcess] = getHistogram(tmpHistograms, kLR, true);
          histograms_missingBJet_memLR_finebin[numBJets_medium][numGenBJets][mbbCut][idxProcess] = getHistogram(tmpHistograms_missingBJet, kLR, true);
        }
      }
    }
  }

  delete inputFile;

  TFile* outputFile_mbb_delphes = new TFile("histogramsForPaper_delphes_mbb.root", "RECREATE");
  outputFile_mbb_delphes->cd();
  TH1* histogram_2genuineBJets_signal_mbb = histograms[-1][2][kDisabled][kSignal_lo][kMbb];
  histogram_2genuineBJets_signal_mbb->SetName("signal_lo_mbb_delphes");
  histogram_2genuineBJets_signal_mbb->Write();
  TH1* histogram_2genuineBJets_background_mbb = histograms[-1][2][kDisabled][kBackground_lo][kMbb];
  histogram_2genuineBJets_background_mbb->SetName("background_lo_mbb_delphes");
  histogram_2genuineBJets_background_mbb->Write();
  delete outputFile_mbb_delphes;

//return;

  std::string labelText_signal = "HH #rightarrow b#bar{b} WW^{*} #rightarrow b#bar{b} l^{+}#nu l^{-}#bar{#nu}";
  std::string labelText_background = "t#bar{t} #rightarrow bW #bar{b}W #rightarrow b l^{+}#nu #bar{b} l^{-}#bar{#nu}";
  //std::string labelText_signal_vs_background = Form("%s vs %s", labelText_signal.data(), labelText_background.data());
  std::string labelText_signal_vs_background = "";

  int showHistograms_canvasSizeX = 1050;
  int showHistograms_canvasSizeY =  950;
  int showHistograms_canvasSizeY_wRatio = 1150;
  double showHistograms_xAxisOffset = 0.96;
  double showHistograms_xAxisOffset_wRatio = 2.90;
  double showHistograms_yAxisOffset = 1.21;
  double showHistograms_yAxisOffset_wRatio = 1.45;
  int showHistograms_colors[4]       = { kGreen - 6, kBlack, kBlue - 7, 28 };
  int showHistograms_markerStyles[4] = { 20, 24, 21, 25 };
  int showHistograms_markerSizes[4]  = { 2, 2, 2, 2 };
  int showHistograms_lineStyles[4]   = { 1, 1, 1, 7 };
  int showHistograms_lineWidths[4]   = { 3, 2, 2, 3 };
  std::vector<std::string> showHistograms_drawOptions = { "hist", "ep", "ep", "hist" };
  std::vector<std::string> showHistograms_legendOptions = { "l", "p", "p", "l" };

  int showHistograms_signal_vs_background_colors[2]       = { kBlack, kRed };
  int showHistograms_signal_vs_background_markerStyles[2] = { 20, 24 };
  int showHistograms_signal_vs_background_markerSizes[2]  = { 2, 2 };
  int showHistograms_signal_vs_background_lineStyles[2]   = { 1, 1 };
  int showHistograms_signal_vs_background_lineWidths[2]   = { 2, 2 };
  std::vector<std::string> showHistograms_signal_vs_background_drawOptions = { "ep", "ep" };
  std::vector<std::string> showHistograms_signal_vs_background_legendOptions = { "p", "p" };

  std::map<int, std::string> xAxisTitle;               // key = idxHistogram
  std::map<int, std::string> xAxisTitle_missingBJet;   // key = idxHistogram
  std::map<int, std::string> xAxisTitle_missingWJet;   // key = idxHistogram
  std::map<int, std::string> xAxisTitle_missingBnWJet; // key = idxHistogram
  std::map<int, std::string> yAxisTitle;               // key = idxHistogram
  std::map<int, std::string> yAxisTitle_missingBJet;   // key = idxHistogram
  std::map<int, std::string> yAxisTitle_missingWJet;   // key = idxHistogram
  std::map<int, std::string> yAxisTitle_missingBnWJet; // key = idxHistogram
  std::map<int, double>      yMin;                     // key = idxHistogram
  std::map<int, double>      yMin_wRatio;              // key = idxHistogram
  std::map<int, double>      yMax;                     // key = idxHistogram

  xAxisTitle[kLR]                  = "P";
  xAxisTitle_missingBJet[kLR]      = "P_{mB}";
  xAxisTitle_missingWJet[kLR]      = "P_{mW}";
  xAxisTitle_missingBnWJet[kLR]    = "P_{mBW}";
  yAxisTitle[kLR]                  = "dN/dP";
  yAxisTitle_missingBJet[kLR]      = "dN/dP_{mB}";
  yAxisTitle_missingWJet[kLR]      = "dN/dP_{mW}";
  yAxisTitle_missingBnWJet[kLR]    = "dN/dP_{mBW}";
  yMin[kLR]                        = 1.1e-5;
  yMin_wRatio[kLR]                 = 1.1e-5;
  yMax[kLR]                        = 1.9e0;

  xAxisTitle[kProbS]               = "log w_{0}";
  xAxisTitle_missingBJet[kProbS]   = "log w_{0}^{mB}";
  xAxisTitle_missingWJet[kProbS]   = "log w_{0}^{mW}";
  xAxisTitle_missingBnWJet[kProbS] = "log w_{0}^{mBW}";
  yAxisTitle[kProbS]               = "dN/dlog w_{0}";
  yAxisTitle_missingBJet[kProbS]   = "dN/dlog w_{0}^{mB}";
  yAxisTitle_missingWJet[kProbS]   = "dN/dlog w_{0}^{mW}";
  yAxisTitle_missingBnWJet[kProbS] = "dN/dlog w_{0}^{mBW}";
  yMin[kProbS]                     = 1.1e-5;
  yMin_wRatio[kProbS]              = 1.1e-5;
  yMax[kProbS]                     = 1.9e0;

  xAxisTitle[kProbB]               = "log w_{1}";
  xAxisTitle_missingBJet[kProbB]   = "log w_{1}^{mB}";
  xAxisTitle_missingWJet[kProbB]   = "log w_{1}^{mW}";
  xAxisTitle_missingBnWJet[kProbB] = "log w_{1}^{mBW}";
  yAxisTitle[kProbB]               = "dN/dlog w_{1}";
  yAxisTitle_missingBJet[kProbB]   = "dN/dlog w_{1}^{mB}";
  yAxisTitle_missingWJet[kProbB]   = "dN/dlog w_{1}^{mW}";
  yAxisTitle_missingBnWJet[kProbB] = "dN/dlog w_{1}^{mBW}";
  yMin[kProbB]                     = 3.1e-5;
  yMin_wRatio[kProbB]              = 3.1e-5;
  yMax[kProbB]                     = 1.9e0;

  xAxisTitle[kMbb]                 = "m_{bb} [GeV]";
  yAxisTitle[kMbb]                 = "dN/dm_{bb} [1/GeV]";
  yMin[kMbb]                       = 3.1e-5;
  yMin_wRatio[kMbb]                = 3.1e-5;
  yMax[kMbb]                       = 1.9e0;

  xAxisTitle[kMll]                 = "m_{ll} [GeV]";
  yAxisTitle[kMll]                 = "dN/dm_{ll} [1/GeV]";
  yMin[kMll]                       = 3.1e-5;
  yMin_wRatio[kMll]                = 3.1e-5;
  yMax[kMll]                       = 1.9e0;

  int showGraphs_canvasSizeX = 1050;
  int showGraphs_canvasSizeY =  950;
  int showGraphs_canvasSizeY_wRatio = 1150;
  double showGraphs_xAxisOffset = 1.18;
  double showGraphs_xAxisOffset_wRatio = 3.00;
  double showGraphs_yAxisOffset = 1.44;
  double showGraphs_yAxisOffset_wRatio = 1.65;
  int showGraphs_colors[4]       = { kGreen - 6, kBlack, kBlue - 7, 28 };
  int showGraphs_markerStyles[4] = { 20, 24, 21, 25 };
  int showGraphs_markerSizes[4]  = { 2, 2, 2, 2 };
  int showGraphs_lineStyles[4]   = { 1, 1, 1, 1 };
  int showGraphs_lineWidths[4]   = { 2, 2, 2, 2 };
  std::vector<std::string> showGraphs_drawOptions = { "Lp", "Lp", "Lp", "Lp" };
  std::vector<std::string> showGraphs_legendOptions = { "lp", "lp", "lp", "lp" };

  std::string legendEntry_mbbPassed = "within m_{bb} window";
  std::string legendEntry_mbbFailed = "outside m_{bb} window";

  if ( makePlots_signal_vs_background ) {
    for ( int idxHistogram = kProbS; idxHistogram <= kMll; ++idxHistogram ) {
      std::string histogramName = getHistogramName(idxHistogram);

      TH1* histogram_2genuineBJets_signal     = histograms[-1][2][kDisabled][kSignal_lo][idxHistogram];
      TH1* histogram_2genuineBJets_background = histograms[-1][2][kDisabled][kBackground_lo][idxHistogram];

      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_2genuineBJets_signal,     "Signal",
        histogram_2genuineBJets_background, "Background",
        nullptr, "",
        nullptr, "",
        showHistograms_signal_vs_background_colors, showHistograms_signal_vs_background_markerStyles, showHistograms_signal_vs_background_markerSizes, 
        showHistograms_signal_vs_background_lineStyles, showHistograms_signal_vs_background_lineWidths, showHistograms_drawOptions,
        0.055, 0.23, 0.78, 0.33, 0.15, showHistograms_signal_vs_background_legendOptions,
        "", 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        xAxisTitle[idxHistogram], showHistograms_xAxisOffset,
        true, yMin[idxHistogram], yMax[idxHistogram], yAxisTitle[idxHistogram], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_delphes_signal_vs_background_%s.pdf", histogramName.data()));

      if ( idxHistogram == kLR ) {
        histogram_2genuineBJets_signal     = histograms_memLR_finebin[-1][2][kDisabled][kSignal_lo];
        histogram_2genuineBJets_background = histograms_memLR_finebin[-1][2][kDisabled][kBackground_lo];
        TGraph* graph_ROC_2genuineBJets_logScale = compGraphROC(
          "graph_ROC_2genuineBJets",
          histogram_2genuineBJets_signal, 
          histogram_2genuineBJets_background, true);

        showGraphs(
          showGraphs_canvasSizeX, showGraphs_canvasSizeY,
          graph_ROC_2genuineBJets_logScale, "",
          nullptr, "",
          nullptr, "",
          nullptr, "",
          showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
          showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
          0.055, 0.23, 0.86, 0.33, 0.08, showGraphs_legendOptions,
          labelText_signal_vs_background, 0.040,
          0.1600, 0.9525, 0.2900, 0.0600,
          10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
          true, 1.e-3, 1.e0, "Background Rate", showGraphs_yAxisOffset, 
          "hh_bbwwMEM_dilepton_delphes_ROC.pdf");
      }
    }
  }

  if ( makePlots_effectOfFakes_by_numBJets ) {
    for ( int idxHistogram = kProbS; idxHistogram <= kMbb; ++idxHistogram ) {
      std::string histogramName = getHistogramName(idxHistogram);

      TH1* histogram_2mediumBJets_signal     = histograms[2][-1][kDisabled][kSignal_lo][idxHistogram];
      TH1* histogram_1mediumBJet_signal      = histograms[1][-1][kDisabled][kSignal_lo][idxHistogram];
      TH1* histogram_0mediumBJets_signal     = histograms[0][-1][kDisabled][kSignal_lo][idxHistogram];
  
      TH1* histogram_leq1mediumBJets_signal = addHistograms(
        Form("histogram_%s_leq1mediumBJets_signal", histogramName.data()),
        histogram_1mediumBJet_signal, 
        histogram_0mediumBJets_signal);

      TH1* histogram_2mediumBJets_background = histograms[2][-1][kDisabled][kBackground_lo][idxHistogram];
      TH1* histogram_1mediumBJet_background  = histograms[1][-1][kDisabled][kBackground_lo][idxHistogram];
      TH1* histogram_0mediumBJets_background = histograms[0][-1][kDisabled][kBackground_lo][idxHistogram];
  
      TH1* histogram_leq1mediumBJets_background = addHistograms(
        Form("histogram_%s_leq1mediumBJets_background", histogramName.data()),
        histogram_1mediumBJet_background, 
        histogram_0mediumBJets_background);

      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_2mediumBJets_signal,    "2 b-tagged jets",
        histogram_leq1mediumBJets_signal, "#leq 1 b-tagged jets",
        nullptr, "",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        xAxisTitle[idxHistogram], showHistograms_xAxisOffset,
        true, yMin[idxHistogram], yMax[idxHistogram], yAxisTitle[idxHistogram], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_numBJets_2histograms_%s_signal.pdf", histogramName.data()));
      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_2mediumBJets_signal, "2 b-tagged jets",
        histogram_1mediumBJet_signal,  "1 b-tagged jet",
        histogram_0mediumBJets_signal, "0 b-tagged jets",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.045, 0.23, 0.66, 0.33, 0.28, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        xAxisTitle[idxHistogram], showHistograms_xAxisOffset,
        true, yMin[idxHistogram], yMax[idxHistogram], yAxisTitle[idxHistogram], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_numBJets_3histograms_%s_signal.pdf", histogramName.data()));
      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_2mediumBJets_background,    "2 b-tagged jets",
        histogram_leq1mediumBJets_background, "#leq 1 b-tagged jets",
        nullptr, "",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        xAxisTitle[idxHistogram], showHistograms_xAxisOffset,
        true, yMin[idxHistogram], yMax[idxHistogram], yAxisTitle[idxHistogram], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_numBJets_2histograms_%s_background.pdf", histogramName.data()));
      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_2mediumBJets_background, "2 b-tagged jets",
        histogram_1mediumBJet_background,  "1 b-tagged jet",
        histogram_0mediumBJets_background, "0 b-tagged jets",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.045, 0.23, 0.66, 0.33, 0.28, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        xAxisTitle[idxHistogram], showHistograms_xAxisOffset,
        true, yMin[idxHistogram], yMax[idxHistogram], yAxisTitle[idxHistogram], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_numBJets_3histograms_%s_background.pdf", histogramName.data()));

      TH1* histogram_missingBJet_2mediumBJets_signal     = histograms_missingBJet[2][-1][kDisabled][kSignal_lo][idxHistogram];
      TH1* histogram_missingBJet_1mediumBJet_signal      = histograms_missingBJet[1][-1][kDisabled][kSignal_lo][idxHistogram];
      TH1* histogram_missingBJet_0mediumBJets_signal     = histograms_missingBJet[0][-1][kDisabled][kSignal_lo][idxHistogram];

      TH1* histogram_missingBJet_geq1mediumBJets_signal = addHistograms(
        Form("histogram_missingBJet_%s_geq1mediumBJets_signal", histogramName.data()),
        histogram_missingBJet_2mediumBJets_signal, 
        histogram_missingBJet_1mediumBJet_signal);

      TH1* histogram_missingBJet_2mediumBJets_background = histograms_missingBJet[2][-1][kDisabled][kBackground_lo][idxHistogram];
      TH1* histogram_missingBJet_1mediumBJet_background  = histograms_missingBJet[1][-1][kDisabled][kBackground_lo][idxHistogram];
      TH1* histogram_missingBJet_0mediumBJets_background = histograms_missingBJet[0][-1][kDisabled][kBackground_lo][idxHistogram];

      TH1* histogram_missingBJet_geq1mediumBJets_background = addHistograms(
        Form("histogram_missingBJet_%s_geq1mediumBJets_background", histogramName.data()),
        histogram_missingBJet_2mediumBJets_background, 
        histogram_missingBJet_1mediumBJet_background);

      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_missingBJet_geq1mediumBJets_signal, "#geq 1 b-tagged jets",
        histogram_missingBJet_0mediumBJets_signal,    "0 b-tagged jets",
        nullptr, "",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        xAxisTitle[idxHistogram], showHistograms_xAxisOffset,
        true, yMin[idxHistogram], yMax[idxHistogram], yAxisTitle[idxHistogram], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_numBJets_2histograms_%s_missingBJet_signal.pdf", histogramName.data()));
      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_missingBJet_2mediumBJets_signal, "2 b-tagged jets",
        histogram_missingBJet_1mediumBJet_signal,  "1 b-tagged jet",
        histogram_missingBJet_0mediumBJets_signal, "0 b-tagged jets",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.045, 0.23, 0.66, 0.33, 0.28, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        xAxisTitle[idxHistogram], showHistograms_xAxisOffset,
        true, yMin[idxHistogram], yMax[idxHistogram], yAxisTitle[idxHistogram], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_numBJets_3histograms_%s_missingBJet_signal.pdf", histogramName.data()));
      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_missingBJet_geq1mediumBJets_background, "#geq 1 b-tagged jets",
        histogram_missingBJet_0mediumBJets_background,    "0 b-tagged jets",
        nullptr, "",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        xAxisTitle[idxHistogram], showHistograms_xAxisOffset,
        true, yMin[idxHistogram], yMax[idxHistogram], yAxisTitle[idxHistogram], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_numBJets_2histograms_%s_missingBJet_background.pdf", histogramName.data()));
      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_missingBJet_2mediumBJets_background, "2 b-tagged jets",
        histogram_missingBJet_1mediumBJet_background,  "1 b-tagged jet",
        histogram_missingBJet_0mediumBJets_background, "0 b-tagged jets",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.045, 0.23, 0.66, 0.33, 0.28, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        xAxisTitle[idxHistogram], showHistograms_xAxisOffset,
        true, yMin[idxHistogram], yMax[idxHistogram], yAxisTitle[idxHistogram], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_numBJets_3histograms_%s_missingBJet_background.pdf", histogramName.data()));

      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_missingBJet_geq1mediumBJets_signal,     "Signal",
        histogram_missingBJet_geq1mediumBJets_background, "Background",
        nullptr, "",
        nullptr, "",
        showHistograms_signal_vs_background_colors, showHistograms_signal_vs_background_markerStyles, showHistograms_signal_vs_background_markerSizes, 
        showHistograms_signal_vs_background_lineStyles, showHistograms_signal_vs_background_lineWidths, showHistograms_drawOptions,
        0.055, 0.23, 0.78, 0.33, 0.15, showHistograms_signal_vs_background_legendOptions,
        "", 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        xAxisTitle[idxHistogram], showHistograms_xAxisOffset,
        true, yMin[idxHistogram], yMax[idxHistogram], yAxisTitle[idxHistogram], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_numBJets_%s_missingBJet.pdf", histogramName.data()));

      if ( idxHistogram == kLR ) {
        histogram_2mediumBJets_signal     = histograms_memLR_finebin[2][-1][kDisabled][kSignal_lo];
        histogram_2mediumBJets_background = histograms_memLR_finebin[2][-1][kDisabled][kBackground_lo];
        TGraph* graph_ROC_2mediumBJets_logScale = compGraphROC(
          "graph_ROC_2mediumBJets",
          histogram_2mediumBJets_signal, 
          histogram_2mediumBJets_background, true);
        histogram_1mediumBJet_signal      = histograms_memLR_finebin[1][-1][kDisabled][kSignal_lo];
        histogram_1mediumBJet_background  = histograms_memLR_finebin[1][-1][kDisabled][kBackground_lo];
        TGraph* graph_ROC_1mediumBJet_logScale = compGraphROC(
          "graph_ROC_1mediumBJet",
          histogram_1mediumBJet_signal, 
          histogram_1mediumBJet_background, true);
        histogram_0mediumBJets_signal     = histograms_memLR_finebin[0][-1][kDisabled][kSignal_lo];
        histogram_0mediumBJets_background = histograms_memLR_finebin[0][-1][kDisabled][kBackground_lo];
        TGraph* graph_ROC_0mediumBJets_logScale = compGraphROC(
          "graph_ROC_0mediumBJets",
          histogram_0mediumBJets_signal, 
          histogram_0mediumBJets_background, true);

        TH1* histogram_geq1mediumBJets_signal = addHistograms(
          Form("histogram_%s_geq1mediumBJets_signal", histogramName.data()),
          histogram_2mediumBJets_signal, 
          histogram_1mediumBJet_signal);
        TH1* histogram_geq1mediumBJets_background = addHistograms(
          Form("histogram_%s_geq1mediumBJets_background", histogramName.data()),
          histogram_2mediumBJets_background, 
          histogram_1mediumBJet_background);
        TGraph* graph_ROC_leq1mediumBJets_logScale = compGraphROC(
          "graph_ROC_leq1mediumBJets",
          histogram_leq1mediumBJets_signal, 
          histogram_leq1mediumBJets_background, true);

        showGraphs(
          showGraphs_canvasSizeX, showGraphs_canvasSizeY,
          graph_ROC_2mediumBJets_logScale,    "2 b-tagged jets",
          graph_ROC_leq1mediumBJets_logScale, "#leq 1 b-tagged jets",
          nullptr, "",
          nullptr, "",
          showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
          showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
          0.055, 0.23, 0.79, 0.33, 0.15, showGraphs_legendOptions,
          labelText_signal_vs_background, 0.040,
          0.1600, 0.9525, 0.2900, 0.0600,
          10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
          true, 2.1e-4, 9.9e0, "Background Rate", showGraphs_yAxisOffset, 
          "hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_numBJets_2graphs_ROC.pdf");
        showGraphs(
          showGraphs_canvasSizeX, showGraphs_canvasSizeY,
          graph_ROC_2mediumBJets_logScale, "2 b-tagged jets",
          graph_ROC_1mediumBJet_logScale,  "1 b-tagged jet",
          graph_ROC_0mediumBJets_logScale, "0 b-tagged jets",
          nullptr, "",
          showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
          showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
          0.045, 0.23, 0.66, 0.33, 0.28, showGraphs_legendOptions,
          labelText_signal_vs_background, 0.040,
          0.1600, 0.9525, 0.2900, 0.0600,
          10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
          true, 2.1e-4, 9.9e0, "Background Rate", showGraphs_yAxisOffset, 
          "hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_numBJets_3graphs_ROC.pdf");

        histogram_missingBJet_2mediumBJets_signal     = histograms_missingBJet_memLR_finebin[2][-1][kDisabled][kSignal_lo];
        histogram_missingBJet_2mediumBJets_background = histograms_missingBJet_memLR_finebin[2][-1][kDisabled][kBackground_lo]; 
        TGraph* graph_ROC_missingBJet_2mediumBJets_logScale = compGraphROC(
          "graph_ROC_missingBJet_2mediumBJets",
          histogram_missingBJet_2mediumBJets_signal, 
          histogram_missingBJet_2mediumBJets_background, true);
        histogram_missingBJet_1mediumBJet_signal     = histograms_missingBJet_memLR_finebin[1][-1][kDisabled][kSignal_lo];
        histogram_missingBJet_1mediumBJet_background = histograms_missingBJet_memLR_finebin[1][-1][kDisabled][kBackground_lo]; 
        TGraph* graph_ROC_missingBJet_1mediumBJet_logScale = compGraphROC(
          "graph_ROC_missingBJet_1mediumBJet",
          histogram_missingBJet_1mediumBJet_signal, 
          histogram_missingBJet_1mediumBJet_background, true);
        histogram_missingBJet_0mediumBJets_signal     = histograms_missingBJet_memLR_finebin[0][-1][kDisabled][kSignal_lo];
        histogram_missingBJet_0mediumBJets_background = histograms_missingBJet_memLR_finebin[0][-1][kDisabled][kBackground_lo]; 
        TGraph* graph_ROC_missingBJet_0mediumBJets_logScale = compGraphROC(
          "graph_ROC_missingBJet_0mediumBJets",
          histogram_missingBJet_0mediumBJets_signal, 
          histogram_missingBJet_0mediumBJets_background, true);

        TH1* histogram_missingBJet_geq1mediumBJets_signal = addHistograms(
          Form("histogram_missingBJet_%s_geq1mediumBJets_signal", histogramName.data()),
          histogram_missingBJet_2mediumBJets_signal, 
          histogram_missingBJet_1mediumBJet_signal);
        TH1* histogram_missingBJet_geq1mediumBJets_background = addHistograms(
          Form("histogram_missingBJet_%s_geq1mediumBJets_background", histogramName.data()),
          histogram_missingBJet_2mediumBJets_background, 
          histogram_missingBJet_1mediumBJet_background);
        TGraph* graph_ROC_missingBJet_geq1mediumBJets_logScale = compGraphROC(
          "graph_ROC_missingBJet_geq1mediumBJets",
          histogram_missingBJet_geq1mediumBJets_signal, 
          histogram_missingBJet_geq1mediumBJets_background, true);

        showGraphs(
          showGraphs_canvasSizeX, showGraphs_canvasSizeY,
          //graph_ROC_missingBJet_geq1mediumBJets_logScale, "#geq 1 b-tagged jets",
          //graph_ROC_missingBJet_0mediumBJets_logScale,    "0 b-tagged jets",
          graph_ROC_missingBJet_geq1mediumBJets_logScale,   "",
          nullptr, "",
          nullptr, "",
          nullptr, "",
          showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
          showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
          0.055, 0.23, 0.79, 0.33, 0.15, showGraphs_legendOptions,
          labelText_signal_vs_background, 0.040,
          0.1600, 0.9525, 0.2900, 0.0600,
          10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
          true, 2.1e-4, 9.9e0, "Background Rate", showGraphs_yAxisOffset, 
          "hh_bbwwMEM_dilepton_effectOfFakes_by_numBJets_ROC_missingBJet.pdf");

        TGraph* graph_ROC_combined_logScale = compGraphROC_combined(
          "graph_ROC_combined",
          histogram_2mediumBJets_signal, 
          histogram_2mediumBJets_background,
          histogram_missingBJet_1mediumBJet_signal, 
          histogram_missingBJet_1mediumBJet_background, true);

        showGraphs(
          showGraphs_canvasSizeX, showGraphs_canvasSizeY,
          graph_ROC_2mediumBJets_logScale,            "P: 2 b-tagged jets",
          graph_ROC_1mediumBJet_logScale,             "P: 1 b-tagged jet",
          graph_ROC_missingBJet_1mediumBJet_logScale, "P_{m}: 1 b-tagged jet",
          graph_ROC_combined_logScale,                "Combination",
          showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
          showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
          0.050, 0.23, 0.73, 0.33, 0.22, showGraphs_legendOptions,
          labelText_signal_vs_background, 0.040,
          0.1600, 0.9525, 0.2900, 0.0600,
          10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
          true, 2.1e-4, 9.9e0, "Background Rate", showGraphs_yAxisOffset, 
          "hh_bbwwMEM_dilepton_effectOfFakes_by_numBJets_ROC_combined.pdf");
      }
    }
  }

  if ( makePlots_effectOfFakes_by_mbb ) {
    for ( int idxHistogram = kProbS; idxHistogram <= kMll; ++idxHistogram ) {
      std::string histogramName = getHistogramName(idxHistogram);

      if ( !(idxHistogram == kProbS || idxHistogram == kProbB || idxHistogram == kLR) ) continue;

      TH1* histogram_mbbPassed_signal     = histograms[2][-1][kEnabled][kSignal_lo][idxHistogram];
      TH1* histogram_mbbFailed_signal     = histograms[2][-1][kInverted][kSignal_lo][idxHistogram];

      TH1* histogram_mbbPassed_background = histograms[2][-1][kEnabled][kBackground_lo][idxHistogram];
      TH1* histogram_mbbFailed_background = histograms[2][-1][kInverted][kBackground_lo][idxHistogram];
  
      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_mbbPassed_signal, legendEntry_mbbPassed,
        histogram_mbbFailed_signal, legendEntry_mbbFailed,
        nullptr, "",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        xAxisTitle[idxHistogram], showHistograms_xAxisOffset,
        true, yMin[idxHistogram], yMax[idxHistogram], yAxisTitle[idxHistogram], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_mbb_%s_signal.pdf", histogramName.data()));
      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_mbbPassed_background, legendEntry_mbbPassed,
        histogram_mbbFailed_background, legendEntry_mbbFailed,
        nullptr, "",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        xAxisTitle[idxHistogram], showHistograms_xAxisOffset,
        true, yMin[idxHistogram], yMax[idxHistogram], yAxisTitle[idxHistogram], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_mbb_%s_background.pdf", histogramName.data()));

      TH1* histogram_missingBJet_mbbPassed_signal     = histograms_missingBJet[2][-1][kEnabled][kSignal_lo][idxHistogram];
      TH1* histogram_missingBJet_mbbFailed_signal     = histograms_missingBJet[2][-1][kInverted][kSignal_lo][idxHistogram];

      TH1* histogram_missingBJet_mbbPassed_background = histograms_missingBJet[2][-1][kEnabled][kBackground_lo][idxHistogram];
      TH1* histogram_missingBJet_mbbFailed_background = histograms_missingBJet[2][-1][kInverted][kBackground_lo][idxHistogram];

      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_missingBJet_mbbPassed_signal, legendEntry_mbbPassed,
        histogram_missingBJet_mbbFailed_signal, legendEntry_mbbFailed,
        nullptr, "",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        xAxisTitle[idxHistogram], showHistograms_xAxisOffset,
        true, yMin[idxHistogram], yMax[idxHistogram], yAxisTitle[idxHistogram], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_mbb_%s_missingBJet_signal.pdf", histogramName.data()));
      showHistograms(
        showHistograms_canvasSizeX, showHistograms_canvasSizeY,
        histogram_missingBJet_mbbPassed_background, legendEntry_mbbPassed,
        histogram_missingBJet_mbbFailed_background, legendEntry_mbbFailed,
        nullptr, "",
        nullptr, "",
        showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
        showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
        0.055, 0.23, 0.77, 0.33, 0.15, showHistograms_legendOptions,
        labelText_signal, 0.055,
        0.1800, 0.9525, 0.2900, 0.0900,
        xAxisTitle[idxHistogram], showHistograms_xAxisOffset,
        true, yMin[idxHistogram], yMax[idxHistogram], yAxisTitle[idxHistogram], showHistograms_yAxisOffset, 
        Form("hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_mbb_%s_missingBJet_background.pdf", histogramName.data()));

      if ( idxHistogram == kLR ) {
        histogram_mbbPassed_signal     = histograms_memLR_finebin[2][-1][kEnabled][kSignal_lo];
        histogram_mbbPassed_background = histograms_memLR_finebin[2][-1][kEnabled][kBackground_lo];
        TGraph* graph_ROC_mbbPassed_logScale = compGraphROC(
          "graph_ROC_mbbPassed",
          histogram_mbbPassed_signal, 
          histogram_mbbPassed_background, true);
        histogram_mbbFailed_signal     = histograms_memLR_finebin[2][-1][kInverted][kSignal_lo];
        histogram_mbbFailed_background = histograms_memLR_finebin[2][-1][kInverted][kBackground_lo];
        TGraph* graph_ROC_mbbFailed_logScale = compGraphROC(
          "graph_ROC_mbbFailed",
          histogram_mbbFailed_signal, 
          histogram_mbbFailed_background, true);

        showGraphs(
          showGraphs_canvasSizeX, showGraphs_canvasSizeY,
          graph_ROC_mbbPassed_logScale, legendEntry_mbbPassed,
          graph_ROC_mbbFailed_logScale, legendEntry_mbbFailed,
          nullptr, "",
          nullptr, "",
          showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
          showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
          0.055, 0.23, 0.79, 0.33, 0.15, showGraphs_legendOptions,
          labelText_signal_vs_background, 0.040,
          0.1600, 0.9525, 0.2900, 0.0600,
          10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
          true, 2.1e-4, 9.9e0, "Background Rate", showGraphs_yAxisOffset, 
          "hh_bbwwMEM_dilepton_delphes_effectOfFakes_by_mbb_ROC.pdf");

        histogram_missingBJet_mbbPassed_signal     = histograms_missingBJet_memLR_finebin[2][-1][kEnabled][kSignal_lo];
        histogram_missingBJet_mbbPassed_background = histograms_missingBJet_memLR_finebin[2][-1][kEnabled][kBackground_lo];
        TGraph* graph_ROC_missingBJet_mbbPassed_logScale = compGraphROC(
          "graph_ROC_missingBJet_mbbPassed",
          histogram_missingBJet_mbbPassed_signal, 
          histogram_missingBJet_mbbPassed_background, true);
        histogram_missingBJet_mbbFailed_signal     = histograms_missingBJet_memLR_finebin[2][-1][kInverted][kSignal_lo];
        histogram_missingBJet_mbbFailed_background = histograms_missingBJet_memLR_finebin[2][-1][kInverted][kBackground_lo];
        TGraph* graph_ROC_missingBJet_mbbFailed_logScale = compGraphROC(
          "graph_ROC_missingBJet_mbbFailed",
          histogram_missingBJet_mbbFailed_signal, 
          histogram_missingBJet_mbbFailed_background, true);

        showGraphs(
          showGraphs_canvasSizeX, showGraphs_canvasSizeY,
          graph_ROC_missingBJet_mbbPassed_logScale, legendEntry_mbbPassed,
          graph_ROC_missingBJet_mbbFailed_logScale, legendEntry_mbbFailed,
          nullptr, "",
          nullptr, "",
          showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
          showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
          0.055, 0.23, 0.79, 0.33, 0.15, showGraphs_legendOptions,
          labelText_signal_vs_background, 0.040,
          0.1600, 0.9525, 0.2900, 0.0600,
          10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
          true, 2.1e-4, 9.9e0, "Background Rate", showGraphs_yAxisOffset, 
          "hh_bbwwMEM_dilepton_effectOfFakes_by_mbb_ROC_missingBJet.pdf");

        TH1* histogram_signal     = histograms_memLR_finebin[2][-1][kDisabled][kSignal_lo];
        TH1* histogram_background = histograms_memLR_finebin[2][-1][kDisabled][kBackground_lo];
        TGraph* graph_ROC_logScale = compGraphROC(
          "graph_ROC",
          histogram_signal, 
          histogram_background, true);

        TFile* outputFile_ROC_combined = new TFile("graphs_ROC_combined_DEBUG.root", "RECREATE");
        outputFile_ROC_combined->cd();        
        graph_ROC_mbbPassed_logScale->SetName("graph_ROC_mbbPassed");
        graph_ROC_mbbPassed_logScale->Write();
        graph_ROC_mbbFailed_logScale->SetName("graph_ROC_mbbFailed");
        graph_ROC_mbbFailed_logScale->Write();
        graph_ROC_missingBJet_mbbFailed_logScale->SetName("graph_ROC_missingBJet_mbbFailed");
        graph_ROC_missingBJet_mbbFailed_logScale->Write();
        graph_ROC_logScale->SetName("graph_ROC");
        graph_ROC_logScale->Write();
        histogram_signal->SetName("histogram_signal");
        histogram_signal->Write();
        histogram_background->SetName("histogram_background");
        histogram_background->Write();
        histogram_mbbPassed_signal->SetName("histogram_mbbPassed_signal");
        histogram_mbbPassed_signal->Write();
        histogram_mbbPassed_background->SetName("histogram_mbbPassed_background");
        histogram_mbbPassed_background->Write();
        histogram_missingBJet_mbbFailed_signal->SetName("histogram_missingBJet_mbbFailed_signal");
        histogram_missingBJet_mbbFailed_signal->Write();
        histogram_missingBJet_mbbFailed_background->SetName("histogram_missingBJet_mbbFailed_background");
        histogram_missingBJet_mbbFailed_background->Write();
        delete outputFile_ROC_combined;

        TGraph* graph_ROC_combined_logScale = compGraphROC_combined(
          "graph_ROC_combined",
          histogram_mbbPassed_signal, 
          histogram_mbbPassed_background,
          histogram_missingBJet_mbbFailed_signal, 
          histogram_missingBJet_mbbFailed_background, true);

        showGraphs(
          showGraphs_canvasSizeX, showGraphs_canvasSizeY,
          graph_ROC_logScale,          "P",
          graph_ROC_combined_logScale, "P and P_{m} combined",
          nullptr, "",
          nullptr, "",
          showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
          showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
          0.055, 0.23, 0.78, 0.33, 0.15, showGraphs_legendOptions,
          labelText_signal_vs_background, 0.040,
          0.1600, 0.9525, 0.2900, 0.0600,
          10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
          true, 2.1e-4, 9.9e0, "Background Rate", showGraphs_yAxisOffset, 
          "hh_bbwwMEM_dilepton_effectOfFakes_by_mbb_ROC_2graphs_combined.pdf");
        showGraphs(
          showGraphs_canvasSizeX, showGraphs_canvasSizeY,
          graph_ROC_mbbPassed_logScale,             Form("P: %s",     legendEntry_mbbPassed.data()),
          graph_ROC_mbbFailed_logScale,             Form("P: %s",     legendEntry_mbbFailed.data()),
          graph_ROC_missingBJet_mbbFailed_logScale, Form("P_{m}: %s", legendEntry_mbbFailed.data()),
          graph_ROC_combined_logScale,              "Combination",
          showGraphs_colors, showGraphs_markerStyles, showGraphs_markerSizes, 
          showGraphs_lineStyles, showGraphs_lineWidths, showGraphs_drawOptions,
          0.050, 0.23, 0.73, 0.33, 0.22, showGraphs_legendOptions,
          labelText_signal_vs_background, 0.040,
          0.1600, 0.9525, 0.2900, 0.0600,
          10, 0., 1.01, "Signal Efficiency", showGraphs_xAxisOffset,
          true, 2.1e-4, 9.9e0, "Background Rate", showGraphs_yAxisOffset, 
          "hh_bbwwMEM_dilepton_effectOfFakes_by_mbb_ROC_4graphs_combined.pdf");
      }
    }
  }

  TFile* outputFile_delphes = new TFile("histogramsForPaper_delphes.root", "RECREATE");
  outputFile_delphes->cd();
  for ( int idxHistogram = kProbS; idxHistogram <= kMll; ++idxHistogram ) {
    std::string histogramName = getHistogramName(idxHistogram);
    TH1* histogram_signal = histograms[2][-1][kDisabled][kSignal_lo][idxHistogram];
    histogram_signal->SetName(Form("signal_lo_%s_delphes", histogramName.data()));
    histogram_signal->Write();
    TH1* histogram_background = histograms[2][-1][kDisabled][kBackground_lo][idxHistogram];
    histogram_background->SetName(Form("background_lo_%s_delphes", histogramName.data()));
    histogram_background->Write();
    if ( idxHistogram == kLR ) {
      TH1* histogram_signal_memLR_finebin = histograms_memLR_finebin[2][-1][kDisabled][kSignal_lo];
      histogram_signal_memLR_finebin->SetName("signal_lo_memLR_finebin_delphes");
      histogram_signal_memLR_finebin->Write();
      TH1* histogram_background_memLR_finebin = histograms_memLR_finebin[2][-1][kDisabled][kBackground_lo];
      histogram_background_memLR_finebin->SetName("background_lo_memLR_finebin_delphes");
      histogram_background_memLR_finebin->Write();
    }
  }
  delete outputFile_delphes;

  clock.Show("makeMEMPerformancePlotsFromNtuples_bbww_dilepton_delphes");
}
