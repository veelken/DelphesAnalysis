
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TBenchmark.h>

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <assert.h>

// CV: The fit function is taken from the presentation
//       https://indico.in2p3.fr/event/6020/contributions/36753/attachments/29528/36363/Talk_beamer.pdf
//    (slide 12)
//     The term ln(x)/x has been replaced by ln(x)/x^c, in order to increase the flexibility of the fit
std::string fitFunctionFormula = "[0] + [1]*log(x)/pow(x, [2])";

bool makePlots_png  = true;
bool makePlots_pdf  = true;
bool makePlots_root = true;

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

TProfile2D* loadHistogram(TFile* inputFile, const std::string& directory, const std::string& histogramName)
{  
  TString histogramName_full = directory.data();
  if ( !histogramName_full.EndsWith("/") ) histogramName_full.Append("/");
  histogramName_full.Append(histogramName.data());
  TProfile2D* histogram = (TProfile2D*)inputFile->Get(histogramName_full.Data());
  if ( !histogram ) {
    std::cerr << "Failed to load histogram = " << histogramName.data() << " from file = " << inputFile->GetName() << " !!" << std::endl;
    assert(0);
  }
  return histogram;
}

void makeControlPlot(double canvasSizeX, double canvasSizeY,
		     TH1* histogram_signal,     TF1* fitFunction_signal,     const std::string& legendEntry_signal,
                     TH1* histogram_background, TF1* fitFunction_background, const std::string& legendEntry_background,
                     double legendTextSize, double legendPosX, double legendPosY, double legendSizeX, double legendSizeY,
                     const std::string& labelText, double labelTextSize, double labelPosX, double labelPosY, double labelSizeX, double labelSizeY,
                     bool useLogScaleX, const std::string& xAxisTitle, double xAxisOffset,
                     bool useLogScaleY, double yMin, double yMax, const std::string& yAxisTitle, double yAxisOffset,
		     const std::string& outputFileName)
{
  const double margin_top      = 0.025;
  const double margin_left     = 0.155;
  const double margin_bottom   = 0.170;
  const double margin_right    = 0.035;

  TCanvas* canvas = new TCanvas("canvas", "canvas", canvasSizeX, canvasSizeY);
  canvas->SetFillColor(10);
  canvas->SetBorderSize(2);
  canvas->SetTopMargin(margin_top);
  canvas->SetLeftMargin(margin_left);
  canvas->SetBottomMargin(margin_bottom);
  canvas->SetRightMargin(margin_right); 
  canvas->SetLogx(useLogScaleX);
  canvas->SetLogy(useLogScaleY);
  canvas->Draw();
  canvas->cd();
  
  assert(histogram_signal);
  histogram_signal->SetFillColor(0);
  histogram_signal->SetFillStyle(0);
  histogram_signal->SetLineColor(1);
  histogram_signal->SetLineStyle(1);
  histogram_signal->SetLineWidth(1);
  histogram_signal->SetMarkerColor(1);
  histogram_signal->SetMarkerStyle(24);
  histogram_signal->SetMarkerSize(2);

  assert(fitFunction_signal);
  fitFunction_signal->SetLineColor(1);
  fitFunction_signal->SetLineStyle(1);
  fitFunction_signal->SetLineWidth(2);

  assert(histogram_background);
  histogram_background->SetFillColor(0);
  histogram_background->SetFillStyle(0);
  histogram_background->SetLineColor(2);
  histogram_background->SetLineStyle(1);
  histogram_background->SetLineWidth(1);
  histogram_background->SetMarkerColor(2);
  histogram_background->SetMarkerStyle(20);
  histogram_background->SetMarkerSize(2);

  assert(fitFunction_background);
  fitFunction_background->SetLineColor(2);
  fitFunction_background->SetLineStyle(1);
  fitFunction_background->SetLineWidth(2);
  
  TH1* histogram_ref = histogram_background;

  TAxis* xAxis = histogram_ref->GetXaxis();
  xAxis->SetTitle(xAxisTitle.data());
  xAxis->SetTitleOffset(xAxisOffset);
  xAxis->SetTitleSize(60);
  xAxis->SetTitleFont(43);
  xAxis->SetLabelSize(0.050);
  xAxis->SetLabelFont(42);
  xAxis->SetTickLength(0.040);
  xAxis->SetNdivisions(505);

  TAxis* yAxis = histogram_ref->GetYaxis();
  yAxis->SetTitle(yAxisTitle.data());
  yAxis->SetTitleOffset(yAxisOffset);
  yAxis->SetTitleSize(60);
  yAxis->SetTitleFont(43);
  if ( yMax > yMin ) {
    yAxis->SetRangeUser(yMin, yMax);
  }
  yAxis->SetLabelSize(0.055);
  yAxis->SetLabelFont(42);
  yAxis->SetTickLength(0.040);  
  yAxis->SetNdivisions(505);

  histogram_ref->SetTitle("");
  histogram_ref->SetStats(false);

  histogram_background->Draw("e1p");
  histogram_signal->Draw("e1psame");
  fitFunction_background->Draw("same");
  fitFunction_signal->Draw("same");

  histogram_ref->Draw("axissame");

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
  if ( legendEntry_signal != "" && legendEntry_background != "" ) {
    legend = new TLegend(legendPosX, legendPosY, legendPosX + legendSizeX, legendPosY + legendSizeY, NULL, "brNDC");
    legend->SetFillColor(10);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(legendTextSize);
    legend->SetTextColor(1);
    legend->SetMargin(0.20);
    legend->AddEntry(histogram_signal,                      legendEntry_signal.data(),      "p");
    legend->AddEntry(fitFunction_signal,     Form("%s fit", legendEntry_signal.data()),     "l");
    legend->AddEntry(histogram_background,                  legendEntry_background.data(),  "p");
    legend->AddEntry(fitFunction_background, Form("%s fit", legendEntry_background.data()), "l");
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

void addToCalibrationFile(ostream* calibrationFile, double min_absEta, double max_absEta, double min_pt, double max_pt, TF1* fitFunction, bool isLast)
{
  std::string text_bin = Form("(abs(eta) > %1.2f && abs(eta) <= %1.2f) * (pt > %1.2f && pt <= %1.2f)", min_absEta, max_absEta, min_pt, max_pt);
  TString text_fitFunction = fitFunctionFormula;
  int numParameter = fitFunction->GetNpar();
  for ( int idxParameter = 0; idxParameter < numParameter; ++idxParameter ) {
    text_fitFunction = text_fitFunction.ReplaceAll(Form("[%i]", idxParameter), Form("%1.5e", fitFunction->GetParameter(idxParameter)));
  }
  text_fitFunction.ReplaceAll("x", "pt");
  text_fitFunction.ReplaceAll(" + -", " - ");
  assert(calibrationFile);
  (*calibrationFile) << text_bin << " * 1./(" << text_fitFunction.Data() << ")";
  if ( !isLast ) (*calibrationFile) << " + ";
  (*calibrationFile) << std::endl;
}

void process(TProfile2D* histogram_signal, TProfile2D* histogram_background, const std::string& histogramLabel, ostream* calibrationFile)
{
  assert(histogram_signal->GetNbinsY() == histogram_background->GetNbinsY());
  TAxis* yAxis = histogram_signal->GetYaxis();
  int numBinsEta = yAxis->GetNbins();
  for ( int idxBinEta = 1; idxBinEta <= numBinsEta; ++idxBinEta ) {
    double min_absEta = yAxis->GetBinLowEdge(idxBinEta);
    double max_absEta = yAxis->GetBinUpEdge(idxBinEta);
    TString etaLabel;
    std::string labelTextEta;
    if ( min_absEta > 0. && max_absEta < 5. ) {
      etaLabel = Form("absEta%1.2fto%1.2f", min_absEta, max_absEta);
      labelTextEta = Form("%1.2f < #eta < %1.2f", min_absEta, max_absEta);
    } else if ( min_absEta > 0. ) {
      etaLabel = Form("absEtaGt%1.2f", min_absEta);
      labelTextEta = Form("#eta > %1.2f", min_absEta);
    } else if ( max_absEta < 5. ) {
      etaLabel = Form("absEtaLt%1.2f", max_absEta);
      labelTextEta = Form("#eta < %1.2f", max_absEta);
    }
    etaLabel = etaLabel.ReplaceAll(".", "_");

    std::string histogramNamePt_signal = Form("%s_%s", histogram_signal->GetName(), etaLabel.Data());
    TH1* histogramPt_signal = histogram_signal->ProfileX(histogramNamePt_signal.data(), idxBinEta, idxBinEta);

    std::string histogramNamePt_background = Form("%s_%s", histogram_background->GetName(), etaLabel.Data());
    TH1* histogramPt_background = histogram_background->ProfileX(histogramNamePt_background.data(), idxBinEta, idxBinEta);

    assert(histogramPt_signal->GetNbinsX() == histogramPt_background->GetNbinsX());
    TAxis* xAxis = histogramPt_signal->GetXaxis();

    double min_pt = xAxis->GetXmin();
    double max_pt = xAxis->GetXmax();

    std::string fitFunctionName_signal = Form("fitFunction_signal_%s", etaLabel.Data());
    TF1* fitFunction_signal = new TF1(fitFunctionName_signal.data(), fitFunctionFormula.data(), min_pt, max_pt);
    fitFunction_signal->SetParameter(0,  1.0);
    fitFunction_signal->SetParameter(1, -0.1);
    fitFunction_signal->SetParameter(2,  1.0);
    histogramPt_signal->Fit(fitFunction_signal);

    std::string fitFunctionName_background = Form("fitFunction_background_%s", etaLabel.Data());
    TF1* fitFunction_background = new TF1(fitFunctionName_background.data(), fitFunctionFormula.data(), min_pt, max_pt);
    fitFunction_background->SetParameter(0,  1.0);
    fitFunction_background->SetParameter(1, -0.1);
    fitFunction_background->SetParameter(2,  1.0);
    histogramPt_background->Fit(fitFunction_background);

    std::string outputFileName = Form("compJetCalibration_%s_%s.png", histogramLabel.data(), etaLabel.Data());
    makeControlPlot(1050, 950,
                    histogramPt_signal,     fitFunction_signal,     "Signal",
                    histogramPt_background, fitFunction_background, "Background",
                    0.050, 0.62, 0.22, 0.31, 0.24,
                    labelTextEta, 0.050, 0.205, 0.87, 0.31, 0.06,
                    true, "p_{T} [GeV]", 1.20,
                    false, 0.70, 1.15, "Response", 1.20,
                    outputFileName);

    addToCalibrationFile(calibrationFile, min_absEta, max_absEta, min_pt, max_pt, fitFunction_signal, idxBinEta == numBinsEta);

    delete fitFunction_signal;
    delete fitFunction_background;
    delete histogramPt_signal;
    delete histogramPt_background;
  }
}

void compJetCalibration()
{
  gROOT->SetBatch(true);

  TH1::AddDirectory(false);
 
  std::string inputFilePath = "/home/veelken/CMSSW_11_1_2/CMSSW_11_1_2/src/hhAnalysis/DelphesAnalysis/test/";
  std::string inputFileName = "calibrate_jets_delphes_all.root";
  TFile* inputFile = openFile(inputFilePath, inputFileName);

  std::string directory_signal          = "signal";
  std::string directory_background      = "background";

  std::string histogramName_jetResponse = "jetPtResponse";
  std::string histogramName_jetVisFrac  = "jetPtVisFrac";
  
  TProfile2D* jetResponse_signal     = loadHistogram(inputFile, directory_signal,     histogramName_jetResponse);
  TProfile2D* jetVisFrac_signal      = loadHistogram(inputFile, directory_signal,     histogramName_jetVisFrac);
  TProfile2D* jetResponse_background = loadHistogram(inputFile, directory_background, histogramName_jetResponse);
  TProfile2D* jetVisFrac_background  = loadHistogram(inputFile, directory_background, histogramName_jetVisFrac);

  std::ofstream* calibrationFile_jetResponse = new std::ofstream("compJetCalibration_jetResponse.txt");
  process(jetResponse_signal, jetResponse_background, "jetResponse", calibrationFile_jetResponse);
  delete calibrationFile_jetResponse;

  std::ofstream* calibrationFile_jetVisFrac = new std::ofstream("compJetCalibration_jetVisFrac.txt");
  process(jetVisFrac_signal,  jetVisFrac_background,  "jetVisFrac",  calibrationFile_jetVisFrac);
  delete calibrationFile_jetVisFrac;

  delete inputFile;
}

