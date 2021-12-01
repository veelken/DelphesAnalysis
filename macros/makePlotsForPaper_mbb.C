
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <assert.h>

bool makePlots_png  = true;
bool makePlots_pdf  = true;
bool makePlots_root = false;

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

TH1* loadHistogram(TFile* inputFile, const std::string& histogramName)
{  
  TH1* histogram = (TH1*)inputFile->Get(histogramName.data());
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName.data() << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  if ( !histogram->GetSumw2N() ) histogram->Sumw2();
  double integral = histogram->Integral(1, histogram->GetNbinsX()); // CV: exclude underflow and overflow bins
  if ( integral > 0. ) histogram->Scale(1./integral);
  return histogram;
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
  histogram1->SetFillColor(0);
  histogram1->SetFillStyle(0);
  histogram1->SetLineColor(colors[0]);
  histogram1->SetLineStyle(lineStyles[0]);
  histogram1->SetLineWidth(lineWidths[0]);
  histogram1->SetMarkerColor(colors[0]);
  histogram1->SetMarkerStyle(markerStyles[0]);
  histogram1->SetMarkerSize(markerSizes[0]);

  assert(histogram2);
  histogram2->SetFillColor(0);
  histogram2->SetFillStyle(0);
  histogram2->SetLineColor(colors[1]);
  histogram2->SetLineStyle(lineStyles[1]);
  histogram2->SetLineWidth(lineWidths[1]);
  histogram2->SetMarkerColor(colors[1]);
  histogram2->SetMarkerStyle(markerStyles[1]);
  histogram2->SetMarkerSize(markerSizes[1]);

  if ( histogram3 ) {
    histogram3->SetFillColor(0);
    histogram3->SetFillStyle(0);
    histogram3->SetLineColor(colors[2]);
    histogram3->SetLineStyle(lineStyles[2]);
    histogram3->SetLineWidth(lineWidths[2]);
    histogram3->SetMarkerColor(colors[2]);
    histogram3->SetMarkerStyle(markerStyles[2]);
    histogram3->SetMarkerSize(markerSizes[2]);
  }

  if ( histogram4 ) {
    histogram4->SetFillColor(0);
    histogram4->SetFillStyle(0);
    histogram4->SetLineColor(colors[3]);
    histogram4->SetLineStyle(lineStyles[3]);
    histogram4->SetLineWidth(lineWidths[3]);
    histogram4->SetMarkerColor(colors[3]);
    histogram4->SetMarkerStyle(markerStyles[3]);
    histogram4->SetMarkerSize(markerSizes[3]);
  }
  
  TAxis* xAxis = histogram1->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);
  xAxis->SetTitleSize(60);
  xAxis->SetTitleFont(43);
  xAxis->SetRangeUser(50., 200.);
  //xAxis->SetLabelOffset(-0.01);
  xAxis->SetLabelSize(0.050);
  xAxis->SetLabelFont(42);
  xAxis->SetTickLength(0.040);
  xAxis->SetNdivisions(505);

  TAxis* yAxis = histogram1->GetYaxis();
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

  histogram1->SetTitle("");
  histogram1->SetStats(false);

  histogram1->Draw(Form("%ssame", drawOptions[0].data()));
  histogram2->Draw(Form("%ssame", drawOptions[1].data()));
  if ( histogram3 ) histogram3->Draw(Form("%ssame", drawOptions[2].data()));
  if ( histogram4 ) histogram4->Draw(Form("%ssame", drawOptions[3].data()));
  histogram1->Draw("axissame");

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
    legend->AddEntry(histogram1, legendEntry1.data(), legendOptions[0].data());
    legend->AddEntry(histogram2, legendEntry2.data(), legendOptions[1].data());
    if ( histogram3 ) legend->AddEntry(histogram3, legendEntry3.data(), legendOptions[2].data());
    if ( histogram4 ) legend->AddEntry(histogram4, legendEntry4.data(), legendOptions[3].data());
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
  delete canvas;
}

void dumpSignalEfficiency_and_BackgroundRate(TH1* histogram_signal, TH1* histogram_background)
{
  double xMin_fit =  50.;
  double xMax_fit = 200.;

  TAxis* xAxis = histogram_signal->GetXaxis();
  xAxis->SetRangeUser(xMin_fit, xMax_fit);

  double mean = histogram_signal->GetMean();
  double rms = histogram_signal->GetRMS();

  xAxis->SetRange(); // reset x-axis range

  std::string fitFunctionFormula = "[0] + [1]*x + [2]*TMath::Gaus(x, [3], [4])";
  std::string fitFunctionName = "fitFunction";
  TF1* fitFunction = new TF1(fitFunctionName.data(), fitFunctionFormula.data(), xMin_fit, xMax_fit);
  fitFunction->SetParameter(0, 0.);
  fitFunction->SetParameter(1, 0.);
  fitFunction->SetParameter(2, histogram_signal->Integral());
  fitFunction->SetParameter(3, mean);
  fitFunction->SetParameter(4, rms);
  histogram_signal->Fit(fitFunction);

  std::cout << "fitted Gaussian: mean = " << fitFunction->GetParameter(3) << ", rms = " << fitFunction->GetParameter(4) << std::endl;

  for ( int numSigma = 1; numSigma <= 2; ++numSigma ) { 
    double xMin_window = mean - numSigma*rms;
    double xMax_window = mean + numSigma*rms;
    double integral_signal_window = histogram_signal->Integral(histogram_signal->FindBin(xMin_window), histogram_signal->FindBin(xMax_window));
    double integral_signal = histogram_signal->Integral(0, histogram_signal->GetNbinsX() + 1);
    //std::cout << "signal: integral = " << integral_signal << ", within mass-window = " << integral_signal_window << std::endl;
    double integral_background_window = histogram_background->Integral(histogram_background->FindBin(xMin_window), histogram_background->FindBin(xMax_window));
    double integral_background = histogram_background->Integral(0, histogram_background->GetNbinsX() + 1);
    //std::cout << "background: integral = " << integral_background << ", within mass-window = " << integral_background_window << std::endl;
    std::cout << "results for " << numSigma << " sigma (" << xMin_window << " < mbb < " << xMax_window << "):" 
              << " signal efficiency = " << integral_signal_window/integral_signal << ","
              << " background rate = " << integral_background_window/integral_background << std::endl;
  }

  delete fitFunction;
}

void makePlotsForPaper_mbb()
{
  gROOT->SetBatch(true);

  TH1::AddDirectory(false);

  std::string inputFilePath = "/home/veelken/CMSSW_11_1_2/CMSSW_11_1_2/src/hhAnalysis/";
  std::string inputFileName = "bbwwMEMPerformanceStudies/macros/histogramsForPaper.root";
  std::string inputFileName_delphes = "DelphesAnalysis/macros/histogramsForPaper_delphes.root";

  std::string histogramName_mbb_signal_unsmeared     = "signal_lo_mbb_unsmeared";
  std::string histogramName_mbb_background_unsmeared = "background_lo_mbb_unsmeared";
  std::string histogramName_mbb_signal_smeared       = "signal_lo_mbb_smeared";
  std::string histogramName_mbb_background_smeared   = "background_lo_mbb_smeared";
  std::string histogramName_mbb_signal_delphes       = "signal_lo_mbb_delphes";
  std::string histogramName_mbb_background_delphes   = "background_lo_mbb_delphes";

  TFile* inputFile = openFile(inputFilePath, inputFileName);
  TH1* histogram_mbb_signal_unsmeared     = loadHistogram(inputFile, histogramName_mbb_signal_unsmeared);
  TH1* histogram_mbb_background_unsmeared = loadHistogram(inputFile, histogramName_mbb_background_unsmeared);
  TH1* histogram_mbb_signal_smeared       = loadHistogram(inputFile, histogramName_mbb_signal_smeared);  
  TH1* histogram_mbb_background_smeared   = loadHistogram(inputFile, histogramName_mbb_background_smeared);

  TFile* inputFile_delphes = openFile(inputFilePath, inputFileName_delphes);
  TH1* histogram_mbb_signal_delphes       = loadHistogram(inputFile_delphes, histogramName_mbb_signal_delphes);
  TH1* histogram_mbb_background_delphes   = loadHistogram(inputFile_delphes, histogramName_mbb_background_delphes);

  std::string labelText_signal = "HH #rightarrow b#bar{b} WW^{*} #rightarrow b#bar{b} l^{+}#nu l^{-}#bar{#nu}";
  std::string labelText_background = "t#bar{t} #rightarrow bW #bar{b}W #rightarrow b l^{+}#nu #bar{b} l^{-}#bar{#nu}";
  //std::string labelText_signal_vs_background = Form("%s vs %s", labelText_signal.data(), labelText_background.data());
  std::string labelText_signal_vs_background = "";

  int showHistograms_canvasSizeX = 1050;
  int showHistograms_canvasSizeY =  950;
  double showHistograms_xAxisOffset = 0.96;
  double showHistograms_yAxisOffset = 1.21;
  int showHistograms_colors[4]       = { kGreen - 6, kBlack, kBlue - 7, 28 };
  int showHistograms_markerStyles[4] = { 20, 24, 21, 25 };
  int showHistograms_markerSizes[4]  = { 2, 2, 2, 2 };
  int showHistograms_lineStyles[4]   = { 1, 1, 1, 7 };
  int showHistograms_lineWidths[4]   = { 3, 2, 2, 3 };
  std::vector<std::string> showHistograms_drawOptions = { "hist", "ep", "ep", "hist" };
  std::vector<std::string> showHistograms_legendOptions = { "l", "p", "p", "l" };

  showHistograms(
    showHistograms_canvasSizeX, showHistograms_canvasSizeY,
    histogram_mbb_signal_unsmeared, "MC truth",
    //histogram_mbb_signal_smeared, "E_{b} smearing",
    histogram_mbb_signal_delphes, "Delphes",
    nullptr, "",
    nullptr, "",
    showHistograms_colors, showHistograms_markerStyles, showHistograms_markerSizes, 
    showHistograms_lineStyles, showHistograms_lineWidths, showHistograms_drawOptions,
    0.055, 0.23, 0.78, 0.33, 0.15, showHistograms_legendOptions,
    "", 0.055,
    0.1800, 0.9525, 0.2900, 0.0900,
    "m_{bb} [GeV]", showHistograms_xAxisOffset,
    true, 3.1e-5, 1.9e0, "dN/dm_{bb} [1/GeV]", showHistograms_yAxisOffset, 
    Form("makePlotsForPaper_mbb.pdf"));

  std::cout << "Signal efficiency and background rate for mass-window cut @ MC truth:" << std::endl;
  dumpSignalEfficiency_and_BackgroundRate(histogram_mbb_signal_unsmeared, histogram_mbb_background_unsmeared);
  std::cout << "Signal efficiency and background rate for mass-window cut @ E_{b} smearing:" << std::endl;
  dumpSignalEfficiency_and_BackgroundRate(histogram_mbb_signal_smeared, histogram_mbb_background_smeared);
  std::cout << "Signal efficiency and background rate for mass-window cut @ Delphes:" << std::endl;
  dumpSignalEfficiency_and_BackgroundRate(histogram_mbb_signal_delphes, histogram_mbb_background_delphes);

  delete inputFile;
  delete inputFile_delphes;
}
