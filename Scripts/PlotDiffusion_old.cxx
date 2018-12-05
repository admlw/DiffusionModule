#include "PlotDiffusion.h"

void increaseError(TH1D *h);
std::pair<int, int> getBinXError(TH1D *driftHisto);

// Define necessary WaveformFunctions (see ../Algorithms/WaveformFunctions.*)
namespace diffmod {

    double WaveformFunctions::convertXToTicks(double xPosition, int WAVEFORM_DRIFT_START_TICK, int WAVEFORM_DRIFT_SIZE, double X_WIDTH){

      //    double tick = WAVEFORM_START_TICK + 17.9427*xPosition;
      double tick = WAVEFORM_DRIFT_START_TICK + (((double)WAVEFORM_DRIFT_SIZE/X_WIDTH)*xPosition);
      return tick;

    }

    double WaveformFunctions::convertTicksToX(int tick, int WAVEFORM_DRIFT_START_TICK, int WAVEFORM_DRIFT_SIZE, double X_WIDTH){

      //double xPos = (tick - 812)/17.9427;
      double xPos = (tick - WAVEFORM_DRIFT_START_TICK)/(WAVEFORM_DRIFT_SIZE/X_WIDTH);
      return xPos;

    }
}

void PlotDiffusion() {
  TFile *fin = new TFile("diffmod.root", "READ");
  return;
}

void increaseError(TH1D* h){
  for (int i = 0; i < h->GetNbinsX(); i++){
    h->SetBinError(i+1, h->GetBinError(i+1)*1.02);
  }
}

std::pair<int, int> getBinXError(TH1D* driftHisto){

  std::pair<int, int> errLowHigh;

  double lowBinVal = driftHisto->GetBinContent(1);
  double highBinVal = driftHisto->GetBinContent(driftHisto->GetNbinsX());
  double midBinVal = driftHisto->GetBinContent((int)driftHisto->GetNbinsX()/2);

  // deal with U-shaped bins
  if (lowBinVal > midBinVal && highBinVal > midBinVal){

    errLowHigh.first = driftHisto->GetBinLowEdge(1);
    errLowHigh.second = driftHisto->GetBinLowEdge(driftHisto->GetNbinsX());

    return errLowHigh;

  }

  // find max bin nearest mean
  double driftMean = driftHisto->GetMean();
  double driftStdDev = driftHisto->GetStdDev();
  double maximum = 0;
  double maximumbin = 0;
  for (int i = 0; i < driftHisto->GetNbinsX(); i++){

    double binVal = driftHisto->GetBinContent(i);
    if (binVal > maximum){
      maximum = binVal;
      maximumbin = i;
    }
    if (binVal == maximum){
      if (std::abs(driftMean - driftHisto->GetBinCenter(i)) < std::abs(driftMean - driftHisto->GetBinCenter(maximumbin))){
        maximum = binVal;
        maximumbin = i;
      }
    }

  }

  int driftMaxBin = driftHisto->GetMaximumBin();
  double driftMaxBinContent = driftHisto->GetBinContent(driftMaxBin);
  double driftIntegral = driftHisto->Integral();

  int lowBin = driftMaxBin;
  int highBin = driftMaxBin;

  double fractionOfTotal = driftMaxBinContent/driftIntegral; 

  while (fractionOfTotal < 0.68){

    // find whether to step left or right
    lowBinVal = driftHisto->GetBinContent(lowBin-1);
    highBinVal = driftHisto->GetBinContent(highBin+1);

    if (lowBinVal > highBinVal && lowBin > 1){
      lowBin = lowBin-1;
      fractionOfTotal+=(lowBinVal/driftIntegral);
    }
    else if (highBin < driftHisto->GetNbinsX()){
      highBin = highBin+1;
      fractionOfTotal+=(highBinVal/driftIntegral);
    }
    else if(highBinVal == 0  && lowBinVal == 0){
      lowBin = lowBin-1;
      highBin = highBin+1;
    }
    else{
      std::cout << "---ERRORS ARE NOT GOOD" << std::endl;
      std::cout << "lowBinVal: " << lowBinVal << " highBinVal: " << highBinVal << std::endl;
      std::cout << "lowBin: " << lowBin << " highBin: " << highBin << std::endl;
      std::cout << "nBins: " << driftHisto->GetNbinsX() << std::endl;
    }

    std::cout << "lowbin: " << lowBin << " highbin: " << highBin << std::endl; 

  }

  errLowHigh.first = driftHisto->GetBinLowEdge(lowBin+1);
  errLowHigh.second = driftHisto->GetBinLowEdge(highBin-1);


  return errLowHigh;

}

void prodPlot(TString* inputFileName){

  //
  // setup
  //

  int WAVEFORM_DRIFT_START_TICK=800;
  int WAVEFORM_DRIFT_SIZE=5400; // End tick
  int NUMBER_DRIFT_BINS=25;
  double X_WIDTH=256.0;
  double DRIFT_VELOCITY=0.111436;

  TFile* fInput = new TFile(*inputFileName, "READ");
  if (!fInput) {
    std::cout << "Bad file" << std::endl;
    return;
  }
  TFile *fOutput = new TFile("fits_"+*inputFileName, "RECREATE");
  if (!fOutput) {
    std::cout << "Bad output file" << std::endl;
  }

  TH1D* histo;
  //TH1D* histoNWvfms = (TH1D*)fInput->Get("h_nWvfmsInBin");
  TH1D* histoCorrectedTicks = (TH1D*)fInput->Get("h_correctedTicks");
  TH1D* histoChisq = new TH1D("histoChisq", ";drift bin; chisq", NUMBER_DRIFT_BINS, 0, NUMBER_DRIFT_BINS);
  //h_wire_in_window
  //h_wire_baseline_corrected

  double driftDistances[NUMBER_DRIFT_BINS];
  double driftDistancesErrsLow[NUMBER_DRIFT_BINS];
  double driftDistancesErrsHigh[NUMBER_DRIFT_BINS];
  double rmsVals[NUMBER_DRIFT_BINS]; 
  double rmsValsErrsLow[NUMBER_DRIFT_BINS];
  double rmsValsErrsHigh[NUMBER_DRIFT_BINS];
  //diffusionFunctions waveFuncs;
  diffmod::WaveformFunctions waveFuncs;
  double lowConv, highConv; // used for conversion from ticks to microseconds
  int NUMBER_TICKS_PER_BIN = WAVEFORM_DRIFT_SIZE / NUMBER_DRIFT_BINS;
  double textSize = 0.045;
  double labelSize = 0.045;

  //
  // logic
  //

  for (int i = 0; i < NUMBER_DRIFT_BINS; i++){

    //if (histoNWvfms->GetBinContent(i+1) < 500) continue;

    fOutput->cd();
    histo = (TH1D*)fInput->Get(Form("diffusionmodule/histo_bin_%i", i));

    // ensure histogram is filled
    if (histo->Integral() == 0){

      std::cout << "bin " << i << " is empty!" << std::endl;

      driftDistances[i] = -1; 
      rmsVals[i] = 0;
      continue;

    }

    std::cout << "Subdrifthisto" << std::endl;
    TH1D* subDriftHisto = new TH1D("subDriftHisto", "", histoCorrectedTicks->GetNbinsX()/NUMBER_DRIFT_BINS, histo->GetBinLowEdge(1)-WAVEFORM_DRIFT_START_TICK, histo->GetBinLowEdge(histo->GetNbinsX()+1)-WAVEFORM_DRIFT_START_TICK);


    for (int j = 0; j < subDriftHisto->GetNbinsX(); j++){
      subDriftHisto->SetBinContent(j, histoCorrectedTicks->GetBinContent((i*46)+j));
    }


    // get histogram for error calculation
    std::cout << "Error histo" << std::endl;
    histoCorrectedTicks->GetXaxis()->SetRangeUser(histo->GetBinLowEdge(0) - WAVEFORM_DRIFT_START_TICK, histo->GetBinLowEdge(histo->GetNbinsX()-1)-WAVEFORM_DRIFT_START_TICK);

    std::pair<int, int> errLowHigh = getBinXError(subDriftHisto);
    subDriftHisto->Delete();

    // Calulate average tick value for the central points of the waveforms
    double binMean = histoCorrectedTicks->GetMean();
    double binStdDev = histoCorrectedTicks->GetStdDev();
    double binNumberWaveforms = histoCorrectedTicks->Integral();

    // perform baseline correction
    std::cout << "Baseline correction" << std::endl;
    TH1D* histo_baselineCorrected = histo;//new TH1D("histo_baselineCorrected", "", histo->GetNbinsX(), 0, histo->GetNbinsX());
    histo_baselineCorrected->SetName(histo->GetName());

    // convert to microseconds
    std::cout << "Convert to us" << std::endl;
    lowConv = histo_baselineCorrected->GetBinLowEdge(1) * 0.5;
    highConv = histo_baselineCorrected->GetBinLowEdge(histo_baselineCorrected->GetNbinsX()+1) * 0.5;
    histo_baselineCorrected->GetXaxis()->SetLimits(lowConv, highConv);

    // find first, last bins above 0.1
    double lowFit =0;
    double highFit = 0;
    double fitLimit = histo_baselineCorrected->GetMaximum() * 0.1;
    for (int i = histo_baselineCorrected->GetMaximumBin(); i >0; i = i-1){

      double evalInBin = histo_baselineCorrected->GetBinContent(i);
      if (evalInBin < fitLimit){

        lowFit = histo_baselineCorrected->GetBinCenter(i);
        break;
      }

    }
    for (int i = histo_baselineCorrected->GetMaximumBin(); i < histo_baselineCorrected->GetNbinsX(); i++){

      double evalInBin = histo_baselineCorrected->GetBinContent(i);
      if (evalInBin < fitLimit){

        highFit = histo_baselineCorrected->GetBinCenter(i);
        break;
      }

    }

    double chisqNdf = 10;
    TF1* fitFunction;
    while (chisqNdf > 1){
      histo_baselineCorrected->Fit("gaus", "q");//, "", lowFit, highFit);
      fitFunction = histo_baselineCorrected->GetFunction("gaus");
      chisqNdf = fitFunction->GetChisquare()/fitFunction->GetNDF();
      increaseError(histo_baselineCorrected);
    }

    std::cout << "chisqNDF:" << chisqNdf << std::endl;
    //if (fitFunction->GetNDF() > 0 && histoNWvfms->GetBinContent(i+1) > 0)
    if (fitFunction->GetNDF() > 0)
      histoChisq->Fill(i, chisqNdf);

    //if (/*chisqNdf > 1.1 ||*/ histoNWvfms->GetBinContent(i+1) < 500) continue;

    histo_baselineCorrected->Write();

    //
    // deal with errors
    //

    // extra 2 here is to convert back to ticks in order to use conversion function
    int centralTick = (i * NUMBER_TICKS_PER_BIN) + WAVEFORM_DRIFT_START_TICK + (NUMBER_TICKS_PER_BIN / 2);
    driftDistances[i] = waveFuncs.convertTicksToX(binMean, 0, WAVEFORM_DRIFT_SIZE, X_WIDTH);
    driftDistancesErrsLow[i] = waveFuncs.convertTicksToX(binStdDev, 0, WAVEFORM_DRIFT_SIZE, X_WIDTH); // driftDistances[i] - waveFuncs.convertTicksToX(errLowHigh.first, 0, WAVEFORM_DRIFT_SIZE, X_WIDTH); 
    driftDistancesErrsHigh[i] = waveFuncs.convertTicksToX(binStdDev,0,WAVEFORM_DRIFT_SIZE, X_WIDTH); // waveFuncs.convertTicksToX(errLowHigh.second, 0, WAVEFORM_DRIFT_SIZE, X_WIDTH) - driftDistances[i];

    rmsVals[i] = std::pow(histo_baselineCorrected->GetFunction("gaus")->GetParameter(2),2);
    //rmsValsErrs[i] = std::pow(histo_baselineCorrected->GetFunction("gaus")->GetParError(2),2);
    rmsValsErrsLow[i] = std::pow(histo_baselineCorrected->GetFunction("gaus")->GetParameter(2),2)*2*(histo_baselineCorrected->GetFunction("gaus")->GetParError(2))/(histo_baselineCorrected->GetFunction("gaus")->GetParameter(2));
    rmsValsErrsHigh[i] = std::pow(histo_baselineCorrected->GetFunction("gaus")->GetParameter(2),2)*2*(histo_baselineCorrected->GetFunction("gaus")->GetParError(2))/(histo_baselineCorrected->GetFunction("gaus")->GetParameter(2));
/*
    if (i == 1) {
      rmsVals[i] = 0;
      rmsValsErrsLow[i] = 0;
      rmsValsErrsHigh[i] = 0;
      driftDistances[i] = 0;

    }
*/
    std::cout << "-- ERRORS FOR POINT " << i << std::endl;
    std::cout << "Drift Distance: " << driftDistances[i] << " - " << driftDistancesErrsLow[i] << " + " << driftDistancesErrsHigh[i] << std::endl;
    std::cout << "sigma: " << rmsVals[i] << " - " << rmsValsErrsLow[i] << " + " << rmsValsErrsHigh[i] << std::endl;

    if (histo_baselineCorrected->GetMaximum() < 1  || driftDistances[i] > X_WIDTH){
      driftDistances[i] = -1;
      driftDistancesErrsLow[i] = 0; 
      driftDistancesErrsHigh[i] = 0; 
      rmsVals[i] = 0;
      rmsValsErrsLow[i] = 0;
      rmsValsErrsHigh[i] = 0;

    }

  }

  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
  //TPad *topPad = new TPad("topPad", "", 0.005, 0.6, 0.995, 0.995);
  //TPad *midPad = new TPad("midPad", "", 0.005, 0.3, 0.995, 0.6);
  //TPad *bottomPad = new TPad("bottomPad", "", 0.005, 0.005, 0.995, 0.3);
  TPad *topPad = new TPad("topPad", "", 0.005, 0.4, 0.995, 0.995);
  TPad *midPad = new TPad("midPad", "", 0.005, 0.005, 0.995, 0.3);
  //TPad *bottomPad = new TPad("bottomPad", "", 0.005, 0.005, 0.995, 0.3);
  topPad->SetBottomMargin(0.01);
  topPad->SetLeftMargin(0.15);
  midPad->SetLeftMargin(0.15);
  midPad->SetTopMargin(0.05);
  midPad->SetBottomMargin(0.05);
  midPad->SetGridy(1);
  //bottomPad->SetTopMargin(0.035);
  //bottomPad->SetLeftMargin(0.15);
  //bottomPad->SetBottomMargin(0.35);
  //bottomPad->SetLineStyle(2);
  topPad->Draw();
  midPad->Draw();
  //bottomPad->Draw();

  topPad->cd();

  TGraphAsymmErrors *gr = new TGraphAsymmErrors(NUMBER_DRIFT_BINS, driftDistances, rmsVals, driftDistancesErrsLow, driftDistancesErrsHigh, rmsValsErrsLow, rmsValsErrsHigh);
  gr->SetTitle("");
  gr->SetName("diffGraph");
  gr->Draw("ap");
  gr->SetMarkerStyle(8);
  gr->SetMarkerSize(0.8);

  gr->Write();
  gr->Fit("pol1", "", "", 0.1, X_WIDTH);
  gr->GetYaxis()->SetRangeUser(0.0001, 8.0);
  gr->GetXaxis()->SetRangeUser(0, X_WIDTH);

  std::cout << "Chisq: " << gr->GetFunction("pol1")->GetChisquare() << std::endl;
  std::cout << "NDF  : " << gr->GetFunction("pol1")->GetNDF() << std::endl;
  std::cout << "Chisq/NDF: " << gr->GetFunction("pol1")->GetChisquare()/gr->GetFunction("pol1")->GetNDF() << std::endl;

  // modify axes
  TAxis *xaxis = gr->GetXaxis();
  xaxis->SetTitle("");
  xaxis->SetLabelSize(0);

  TAxis *yaxis = gr->GetYaxis();
  yaxis->SetTitle("#sigma^{2} (#mus^{2})");
  yaxis->SetTitleSize(textSize * 1.0/0.4);
  yaxis->SetLabelSize(labelSize * 1.0/0.4);
  yaxis->SetTitleOffset(0.7);

  //
  // acutal diffusion calculation
  //

  double diffusionValue = gr->GetFunction("pol1")->GetParameter(1) * DRIFT_VELOCITY * DRIFT_VELOCITY * DRIFT_VELOCITY * 1000000/2;
  double diffusionValueErr = gr->GetFunction("pol1")->GetParError(1) * DRIFT_VELOCITY * DRIFT_VELOCITY * DRIFT_VELOCITY * 1000000/2;
  std::cout << "Diffusion value: " << diffusionValue << " +/- " << diffusionValueErr << std::endl;

  // include TPaveText on plot
  TPaveText *pt = new TPaveText(0.17, 0.73, 0.5, 0.88, "NDC");
  pt->SetTextAlign(12);
  pt->SetFillStyle(0);
  pt->SetTextSize(0.03*1.0/0.4);
  pt->SetBorderSize(0);
  TString diffusions = Form("Measured D_{L}: %0.2f +/- %0.2f cm^{2}/s", diffusionValue, diffusionValueErr);
  pt->AddText(diffusions);
  pt->Draw("same");

  midPad->cd();
  double rmsValsBottom[NUMBER_DRIFT_BINS]; 
  double rmsValsErrsLowBottom[NUMBER_DRIFT_BINS];
  double rmsValsErrsHighBottom[NUMBER_DRIFT_BINS];
  TF1* f2 = (TF1*)gr->GetFunction("pol1");

  for (int i = 0; i < NUMBER_DRIFT_BINS; i++){

    rmsValsBottom[i] = (rmsVals[i] - (f2->Eval(driftDistances[i])))/f2->Eval(driftDistances[i]);
    rmsValsErrsLowBottom[i] = (double)rmsValsErrsLow[i]/f2->Eval(driftDistances[i]);
    rmsValsErrsHighBottom[i] = (double)rmsValsErrsHigh[i]/f2->Eval(driftDistances[i]);

  }

  TGraphAsymmErrors* gr2 = new TGraphAsymmErrors(NUMBER_DRIFT_BINS, driftDistances, rmsValsBottom, driftDistancesErrsLow, driftDistancesErrsHigh, rmsValsErrsLowBottom, rmsValsErrsHighBottom);
  gr2->GetYaxis()->SetRangeUser(-0.1, 0.1);
  gr2->SetMarkerStyle(8);
  gr2->SetMarkerSize(0.8);
  gr2->SetTitle("");
  gr2->GetYaxis()->SetTitle("(#sigma^{2} - Fit)/Fit");
  gr2->GetXaxis()->SetLimits(0, 256);
  gr2->GetXaxis()->SetTitleSize(0);
  gr2->GetXaxis()->SetLabelSize(0);
  gr2->GetYaxis()->SetNdivisions(505);
  gr2->GetYaxis()->SetTitleOffset(0.5);
  gr2->GetYaxis()->SetTitleSize(textSize * 1.0/0.3);
  gr2->GetYaxis()->SetLabelSize(labelSize * 1.0/0.3);

  gr2->Draw("ap");

  TF1* f3 = new TF1("f3", "0.", 0, 256.);
  f3->Draw("same");

  /*
  bottomPad->cd();
  gStyle->SetOptStat(0);
  histoNWvfms->GetXaxis()->SetTitle("Drift Distance (cm)");
  histoNWvfms->SetFillColor(kBlack);
  histoNWvfms->SetLineColor(kBlack);
  histoNWvfms->GetXaxis()->SetTitleSize(textSize * 1.0/0.3);
  histoNWvfms->GetXaxis()->SetLabelSize(labelSize * 1.0/0.3);
  histoNWvfms->GetYaxis()->SetTitle("#waveforms");
  histoNWvfms->GetYaxis()->SetNdivisions(304);
  histoNWvfms->GetYaxis()->SetTitleOffset(0.55);
  histoNWvfms->SetMinimum(0);
  histoNWvfms->GetYaxis()->SetTitleSize(textSize * 1.0/0.3);
  histoNWvfms->GetYaxis()->SetLabelSize(labelSize * 1.0/0.3);
  histoNWvfms->GetXaxis()->SetLimits(0, 256);

  histoNWvfms->Draw();
  */

  c1->SaveAs("DiffusionPlot.pdf");

  TCanvas *c2 = new TCanvas("c2", "c2", 500, 500);
  c2->SetLogy();
  histoChisq->Draw();
  c2->SaveAs("chisq.png");

  fOutput->Close();
  fInput->Close();

}

int main(int argv, char** argc){

  TString* inputFileName = new TString(argc[1]);

  prodPlot(inputFileName);
  return 0;

}
