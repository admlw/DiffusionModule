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
  const double DRIFT_VELOCITY=0.111436;

  TFile* fInput = new TFile(*inputFileName, "READ");
  if (!fInput) {
    std::cout << "Bad file" << std::endl;
    return;
  }
  TFile *fOutput = new TFile("fits_"+*inputFileName, "RECREATE");
  if (!fOutput) {
    std::cout << "Bad output file" << std::endl;
  }

  TH1D* waveformHist;
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
  double textSize = 0.045;
  double labelSize = 0.045;

  // Loop over bins, do the things
  for (int i = 0; i < NUMBER_DRIFT_BINS; i++){

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

    // Find fit range for histogram 
    // Stop fit at 10% of maximum value (arbitrary?)
    TH1D *waveformHist_us = waveformHist;
    /*
    TH1D *waveformHist_us = new TH1D("waveformHist_us", "", 
        waveformHist->GetNbinsX(), 
        waveformHist->GetBinLowEdge(1)*0.5,
        waveformHist->GetBinLowEdge(waveformHist->GetNbinsX()+1)*0.5
    );
    cout << "Making new hist with bounds " << waveformHist->GetBinLowEdge(1)*0.5 
      << " to " << waveformHist->GetBinLowEdge(waveformHist->GetNbinsX()+1)*0.5 << endl;
    for (int i = 0; i < waveformHist_us->GetNbinsX(); i++) {
      waveformHist_us->Fill(waveformHist->GetXaxis()->GetBinCenter(i+1)*0.5);
      cout << "Filling with " << waveformHist->GetXaxis()->GetBinLowEdge(i+1)*0.5 << endl;
    }
    */
    waveformHist_us->SetName(waveformHist->GetName() );
    double lowConv = waveformHist_us->GetBinLowEdge(1)*0.5;
    double highConv = waveformHist_us->GetBinLowEdge(waveformHist_us->GetNbinsX()+1)*0.5;
    waveformHist_us->GetXaxis()->SetLimits(lowConv, highConv);
    double lowFit = 0, highFit = 0;
    double fitLimit = waveformHist_us->GetMaximum()*0.1;

    for (int i = waveformHist_us->GetMaximumBin(); i > 0; i--) {
      if (waveformHist_us->GetBinContent(i) < fitLimit) {
        lowFit = waveformHist_us->GetBinCenter(i);
        break;
      }
    }
    for (int i = waveformHist_us->GetMaximumBin(); i < waveformHist_us->GetNbinsX(); i++) {
      if (waveformHist_us->GetBinContent(i) < fitLimit) {
        highFit = waveformHist_us->GetBinCenter(i);
        break;
      }
    }

    // Define fit range around the peak, else ROOT might derp out
    double chisqNdf = 10;
    //TF1 *fit = new TF1("fit", "gaus", lowFit, highFit);
    TF1 *fit;
    cout << "lowFit = " << lowFit << endl;
    cout << "highFit = " << highFit << endl;
    cout << "Max = " << waveformHist_us->GetXaxis()->GetBinCenter(waveformHist_us->GetMaximumBin() ) << endl;
    // Chi2 inflation
    while (chisqNdf > 1) {
      waveformHist_us->Fit("gaus", "qr");
      fit = waveformHist_us->GetFunction("gaus");
      chisqNdf = fit->GetChisquare()/fit->GetNDF();
      increaseError(waveformHist_us);
    }
    
    TCanvas *c = new TCanvas;
    waveformHist_us->Draw("hist");
    fit->SetLineColor(kRed);
    fit->Draw("same");
    TString bin = Form("%i", i);
    c->SaveAs("testHist_"+bin+".png", "PNG");

    cout << "Fit param 0: " << fit->GetParameter(0) << endl;
    cout << "Fit param 1: " << fit->GetParameter(1) << endl;
    cout << "Fit param 2: " << fit->GetParameter(2) << endl;

    // Get drift distance
    //double driftDistance = waveFuncs.convertTicksToX(peakTickVal);
    double driftDistance = waveFuncs.convertTicksToX(fit->GetParameter(1) );
    driftDistances[i] = driftDistance;
    driftDistancesErrs[i] = 0.5*X_WIDTH/NUMBER_DRIFT_BINS; // x-error is 1/2 bin width for now
    cout << "Drift distance = " << driftDistances[i] << " +/- " << driftDistancesErrs[i] << endl;

    // Get sigma^2 (pulse width squared)
    sigmaVals[i] = std::pow(fit->GetParameter(2)/2, 2); // Divide by two to convert to microseconds
    sigmaValsErrs[i] = 2*fit->GetParError(2); 
    // Note: factor of 2 from squaring sigma 
    // (see http://ipl.physics.harvard.edu/wp-uploads/2013/03/PS3_Error_Propagation_sp13.pdf)
    cout << "Sigma^2 = " << sigmaVals[i] << " +/- " << sigmaValsErrs[i] << endl;

  }

  fOutput->Close();
  fInput->Close();

}

int main(int argv, char** argc){

  TString* inputFileName = new TString(argc[1]);

  makePlot(inputFileName);
  return 0;

}
