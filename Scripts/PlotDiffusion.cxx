#include "PlotDiffusion.h"
using std::cout;
using std::endl;

void increaseError(TH1D *h);

// Define necessary WaveformFunctions (see ../Algorithms/WaveformFunctions.*)
namespace diffmod {

    double WaveformFunctions::convertXToTicks(double xPosition){

      const int WAVEFORM_DRIFT_START_TICK=800;
      const int WAVEFORM_DRIFT_SIZE=4600; // End tick (5400 - 800)
      const double X_WIDTH=256.0;

      double tick = WAVEFORM_DRIFT_START_TICK + (((double)WAVEFORM_DRIFT_SIZE/X_WIDTH)*xPosition);
      return tick;

    }

    double WaveformFunctions::convertTicksToX(int tickVal) {

      const int WAVEFORM_DRIFT_START_TICK=800;
      const int WAVEFORM_DRIFT_SIZE=4600; // End tick (5400 - 800)
      const double X_WIDTH=256.0;

      double xPos = (tickVal - WAVEFORM_DRIFT_START_TICK) / (WAVEFORM_DRIFT_SIZE/X_WIDTH); 
      return xPos;

    }
}

void increaseError(TH1D* h){
  for (int i = 0; i < h->GetNbinsX(); i++){
    h->SetBinError(i+1, h->GetBinError(i+1)*1.02);
  }
}


void makePlot(TString* inputFileName){

  // Initialization 
  const int WAVEFORM_DRIFT_START_TICK=800;
  const int WAVEFORM_DRIFT_SIZE=4600; // End tick (5400 - 800)
  const double X_WIDTH=256.0;
  const int NUMBER_DRIFT_BINS=25;
  const double DRIFT_VELOCITY=0.1098;

  TFile* fInput = new TFile(*inputFileName, "READ");
  if (!fInput) {
    std::cout << "Bad file" << std::endl;
    return;
  }
  TFile *fOutput = new TFile("fits_"+*inputFileName, "RECREATE");
  if (!fOutput) {
    std::cout << "Bad output file" << std::endl;
  }

  //TH1D* histoNWvfms = (TH1D*)fInput->Get("h_nWvfmsInBin");
  //TH1D* h_peakTime = (TH1D*)fInput->Get("hit_peak_time");
  //TH1D* histoChisq = new TH1D("histoChisq", ";drift bin; chisq", NUMBER_DRIFT_BINS, 0, NUMBER_DRIFT_BINS);
  //h_wire_in_window
  //h_wire_baseline_corrected

  double driftDistances[NUMBER_DRIFT_BINS];
  //double driftDistancesErrsLow[NUMBER_DRIFT_BINS];
  //double driftDistancesErrsHigh[NUMBER_DRIFT_BINS];
  double driftDistancesErrs[NUMBER_DRIFT_BINS]; // bin width/2 for now
  double sigmaVals[NUMBER_DRIFT_BINS]; 
  //double sigmaValsErrsLow[NUMBER_DRIFT_BINS];
  //double sigmaValsErrsHigh[NUMBER_DRIFT_BINS];
  double sigmaValsErrs[NUMBER_DRIFT_BINS]; // Error on Gaussian fit? Maybe?

  diffmod::WaveformFunctions waveFuncs;
  double lowConv, highConv; // Used for conversion from ticks to microseconds
  int NUMBER_TICKS_PER_BIN = WAVEFORM_DRIFT_SIZE / NUMBER_DRIFT_BINS;
  double titleSize = 0.045;
  double labelSize = 0.045;

  TH1D* waveformHist;

  // Loop over bins, do the things
  for (int i = 0; i < NUMBER_DRIFT_BINS; i++){

    std::cout << "************* BIN NUMBER " << i << " **********************" << std::endl;

    //if (histoNWvfms->GetBinContent(i+1) < 500) continue;

    fOutput->cd();
    waveformHist = (TH1D*)fInput->Get(Form("diffusionmodule/histo_bin_%i", i));
    // Ensure histogram is filled
    if (waveformHist->Integral() == 0){
      std::cout << "bin " << i << " is empty!" << std::endl;
      driftDistances[i] = -1; 
      sigmaVals[i] = 0;
      continue;
    }

    // Convert tick histogram to microseconds, find fit range 
    double lowConv = waveformHist->GetBinLowEdge(1)*0.5;
    double highConv = waveformHist->GetBinLowEdge(waveformHist->GetNbinsX()+1)*0.5;
    waveformHist->GetXaxis()->SetLimits(lowConv, highConv);
    double lowFit = 0, highFit = 0;
    
    // Stop fit at 10% of maximum value (arbitrary?)
    double fitLimit = waveformHist->GetMaximum()*0.1;
    for (int i = waveformHist->GetMaximumBin(); i > 0; i--) {
      if (waveformHist->GetBinContent(i) < fitLimit) {
        lowFit = waveformHist->GetBinCenter(i);
        break;
      }
    }
    for (int i = waveformHist->GetMaximumBin(); i < waveformHist->GetNbinsX(); i++) {
      if (waveformHist->GetBinContent(i) < fitLimit) {
        highFit = waveformHist->GetBinCenter(i);
        break;
      }
    }

    // Error inflation
    double chisqNdf = 10;
    TF1 *gausfit = new TF1("gausfit", "gaus");
    // Chi2 inflation
    while (chisqNdf > 1) {
      waveformHist->Fit(gausfit, "q", "", lowFit, highFit);
      chisqNdf = waveformHist->GetFunction("gausfit")->GetChisquare()/gausfit->GetNDF();
      increaseError(waveformHist);
    }
    
    TCanvas *c_test = new TCanvas("c_test", "c_test", 1000, 700);
    gStyle->SetOptFit(1);
    waveformHist->Draw("hist");
    TF1* fitted_function = waveformHist->GetFunction("gausfit");
    if (!fitted_function) {
        std::cout << "Bad fit function in test histogram" << std::endl;
        break;
    }
    fitted_function->SetLineColor(kRed);
    fitted_function->Draw("same");
    TString bin = Form("%i", i);
    c_test->SaveAs("testHist_"+bin+".png", "PNG");
    delete c_test;

    // Get drift distance; factor of 2 to get microsecond value into ticks for conversion
    double driftDistance = waveFuncs.convertTicksToX(gausfit->GetParameter(1)*2 );
    driftDistances[i] = driftDistance;
    //std::cout << "Drift distance: " << driftDistance << std::endl;
    driftDistancesErrs[i] = 0.5*X_WIDTH/NUMBER_DRIFT_BINS; // x-error is 1/2 bin width for now

    // Get sigma^2 (pulse width squared)
    sigmaVals[i] = std::pow(gausfit->GetParameter(2), 2); 
    std::cout << "Sigma^2 val: " << sigmaVals[i] << std::endl;
    // Note: factor of 2 in error calc from squaring sigma 
    // (see http://ipl.physics.harvard.edu/wp-uploads/2013/03/PS3_Error_Propagation_sp13.pdf)
    sigmaValsErrs[i] = 2*gausfit->GetParError(2); 

  }

  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 700);
  TPad *topPad = new TPad("topPad", "", 0.005, 0.005, 0.995, 0.995);
  //TPad *midPad = new TPad("midPad", "", 0.005, 0.3, 0.995, 0.6);
  //TPad *botPad = new TPad("botPad", "", 0.005, 0.005, 0.995, 0.3);
  //topPad->cd();
  // Diffusion plot
  TGraphErrors *g1 = new TGraphErrors(NUMBER_DRIFT_BINS,
      driftDistances, sigmaVals, 
      driftDistancesErrs, sigmaValsErrs
  );
  std::cout << "Making TGraph" << std::endl;
  g1->SetTitle("");
  g1->GetXaxis()->SetTitle("Drift Distance (cm)");
  g1->GetXaxis()->SetTitleSize(titleSize);
  g1->GetXaxis()->SetLabelSize(labelSize);
  g1->GetYaxis()->SetTitle("#sigma_{t}^{2} (#mus)");
  g1->GetYaxis()->SetTitleSize(titleSize);
  g1->GetYaxis()->SetTitleOffset(0.9);
  g1->GetYaxis()->SetLabelSize(labelSize);
  g1->Draw("ap");
  g1->SetMarkerStyle(8);
  g1->SetMarkerSize(0.8);
  g1->GetYaxis()->SetRangeUser(0.001, 10.);
  g1->GetXaxis()->SetRangeUser(0, X_WIDTH);

  // Make sigma plot for verification
  TCanvas *c_sig = new TCanvas("c_sig", "c_sig", 1000, 700);
  TH2D *h_sig = (TH2D*)fInput->Get("diffusionmodule/h_driftVsigma");
  //h_sig->SetMarkerColor(kAzure+2);
  h_sig->Draw("colz");
  g1->Draw("p same");
  c_sig->SaveAs("sigmaPlot.png", "PNG");
  delete c_sig;


  // Linear fit to diffusion plot
  TF1* polFit = new TF1("polfit", "pol1");
  g1->Fit("polfit", "", "", 0.1, X_WIDTH);
  g1->GetFunction("polfit")->SetLineColor(kRed);
  g1->GetFunction("polfit")->Draw("same");

  // Get diffusion value from slope of linear fit
  double diffusionValue = polFit->GetParameter(1)*std::pow(DRIFT_VELOCITY, 3) * 1000000/2;
  std::cout << "Slope = " << polFit->GetParameter(1) << std::endl;
  std::cout << "Slope * v^3 = " << polFit->GetParameter(1)*std::pow(DRIFT_VELOCITY, 3) << std::endl; 
  double diffusionValueErr = polFit->GetParError(1)*std::pow(DRIFT_VELOCITY, 3) * 1000000/2;
  cout << "Diffusion value: " << diffusionValue << " +/- " << diffusionValueErr << endl; 

  TPaveText *pt = new TPaveText(0.12, 0.75, 0.5, 0.9, "NDC");
  pt->SetTextAlign(12);
  pt->SetFillStyle(0);
  pt->SetTextSize(0.045);
  pt->SetBorderSize(0);
  TString diffTextMeas = Form("Measured D_{L}: %0.2f +/- %0.2f cm^{2}/s", diffusionValue, diffusionValueErr);
  TString diffTextInp = Form("Input D_{L}: 6.20 cm^{2}/s");
  pt->AddText(diffTextMeas);
  pt->AddText(diffTextInp);
  pt->Draw("same");

  c1->SaveAs("DiffusionPlot.png", "PNG");
  delete c1;

  fOutput->Close();
  fInput->Close();

}

int main(int argv, char** argc){

  TString* inputFileName = new TString(argc[1]);

  makePlot(inputFileName);
  return 0;

}
