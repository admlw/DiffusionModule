#include "PlotDiffusion.h"
using std::cout;
using std::endl;

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

        //std::cout << "lowbin: " << lowBin << " highbin: " << highBin << std::endl; 
 
    }
  
    errLowHigh.first = driftHisto->GetBinLowEdge(lowBin+1);
    errLowHigh.second = driftHisto->GetBinLowEdge(highBin-1);
  
    return errLowHigh;
}


void makePlot(TString* inputFileName){

  // Initialization 
  const int waveformDriftStartTick=800;
  const int waveformDriftEndTick=5400;
  const int WAVEFORM_DRIFT_SIZE=waveformDriftEndTick-waveformDriftStartTick; // End tick (5400 - 800)
  const int NUMBER_DRIFT_BINS=25;
  int NUMBER_TICKS_PER_BIN = WAVEFORM_DRIFT_SIZE / NUMBER_DRIFT_BINS;
  const int minTime = waveformDriftStartTick/2; // 400 microseconds
  const int maxTime = waveformDriftEndTick/2; // 2700 microseconds
  const double DRIFT_VELOCITY=0.1098;

  TString dirname = "/uboone/data/users/amogan/v08_00_00_19/output_diffmod_files/";
  TString *dir = &dirname;
  TFile *fInput = new TFile(*dir+*inputFileName, "READ");
  if (!fInput) {
    std::cout << "Bad file" << std::endl;
    return;
  }
    
  TFile *fOutput = new TFile(*dir+"fits_"+*inputFileName, "RECREATE");
  if (!fOutput) {
    std::cout << "Bad output file" << std::endl;
    return;
  }

  TCanvas *c = new TCanvas();
  TH1D *h_correctedTicks = new TH1D("h_correctedTicks", "", 1150, 800, 5400); // Range is 800 to 5400
  TTree *t = (TTree*)fInput->Get("DiffusionModule/difftree");
  if (!t) {
      cout << "Bad tree" << endl;
      return;
  }
  c->cd();
  t->Draw("hit_peak_time>>h_correctedTicks");
  c->SaveAs("hit_peak_time.png", "PNG");

  TH1D *histoNWvfms = (TH1D*)fInput->Get("DiffusionModule/h_nWvfmsInBin");
  if (!histoNWvfms) {
    std::cout << "Bad waveform hist" << std::endl;
    return;
  }

  double driftTimes[NUMBER_DRIFT_BINS];
  double driftTimesErrs[NUMBER_DRIFT_BINS]; // bin width/2 for now
  double sigmaVals[NUMBER_DRIFT_BINS]; 
  double sigmaValsErrs[NUMBER_DRIFT_BINS]; // Error on Gaussian fit? Maybe?

  diffmod::WaveformFunctions waveFuncs;
  double lowConv, highConv; // Used for conversion from ticks to microseconds
  double titleSize = 0.1;
  double textSize = 0.07;
  double labelSize = 0.07;

  TH1D* waveformHist;

  // Loop over bins, do the things
  for (int i = 0; i < NUMBER_DRIFT_BINS; i++){

    std::cout << "-------------------------------" << std::endl;
    std::cout << "BIN " << i << std::endl;

    if (histoNWvfms->GetBinContent(i+1) < 500) {
        std::cout << "Skipping bin " << i << " due to low stats" << std::endl;
        continue;
    }

    fOutput->cd();
    waveformHist = (TH1D*)fInput->Get(Form("DiffusionModule/summed_waveform_bin_%i", i));
    // Ensure histogram is filled
    if (waveformHist->Integral() == 0){
      std::cout << "bin " << i << " is empty!" << std::endl;
      driftTimes[i] = -1; 
      sigmaVals[i] = 0;
      continue;
    }

    // Convert tick histogram to microseconds, find fit range 
    double lowConv = waveformHist->GetBinLowEdge(1)*0.5;
    double highConv = waveformHist->GetBinLowEdge(waveformHist->GetNbinsX()+1)*0.5;
    waveformHist->GetXaxis()->SetLimits(lowConv, highConv);
    double lowFit = 0, highFit = 0;
    
    // Stop fit at 10% of maximum value 
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

    TF1 *gausfit = new TF1("gausfit", "gaus");
    waveformHist->Fit(gausfit, "q", "", lowFit, highFit);

    // Error inflation
    double chisqNdf = 10;
    // Chi2 inflation
    while (chisqNdf > 1) {
      waveformHist->Fit(gausfit, "q", "", lowFit, highFit);
      chisqNdf = waveformHist->GetFunction("gausfit")->GetChisquare()/gausfit->GetNDF();
      increaseError(waveformHist);
    }

    double chisq = waveformHist->GetFunction("gausfit")->GetChisquare();
    double Ndf = gausfit->GetNDF();
    std::cout << "chi^2 = " <<  chisq << std::endl;
    std::cout << "NDF = " <<  Ndf << std::endl;
    std::cout << "chi^2/NDF = " <<  (double)chisq/Ndf << std::endl;

    /*
    std::cout << "Setting waveformHist bounds to " << waveformHist->GetBinLowEdge(0) - waveformDriftStartTick << 
                 " and " << waveformHist->GetBinLowEdge(waveformHist->GetNbinsX()-1)-waveformDriftStartTick << std::endl;
    h_correctedTicks->GetXaxis()->SetRangeUser(waveformHist->GetBinLowEdge(0) - waveformDriftStartTick, 
                                               waveformHist->GetBinLowEdge(waveformHist->GetNbinsX()-1)-waveformDriftStartTick);
    */
    // Factors of two to convert
    h_correctedTicks->GetXaxis()->SetLimits(lowConv, highConv);
    std::cout << "Setting waveformHist bounds to " << waveformHist->GetBinLowEdge(1) << 
                 " and " << waveformHist->GetBinLowEdge(waveformHist->GetNbinsX()-1) << std::endl;
    h_correctedTicks->GetXaxis()->SetRangeUser(waveformHist->GetBinLowEdge(1), 
                                               waveformHist->GetBinLowEdge(waveformHist->GetNbinsX()-1));
    double binMean = h_correctedTicks->GetMean();
    double binStdDev = h_correctedTicks->GetStdDev();
    int N = h_correctedTicks->GetEntries();
    
    // Get drift time from truncated mean of sigma distribution in each bin
    // Error is (1/sqrt(N)) * 0.5/2 (half a tick width, if not using hit information)
    driftTimes[i] = binMean;
    driftTimesErrs[i] = (1/sqrt(N) ) * (0.5/2.);

    // Get sigma^2 (pulse width squared)
    sigmaVals[i] = std::pow(gausfit->GetParameter(2), 2); 
    sigmaValsErrs[i] = sqrt(2)*sigmaVals[i]*(waveformHist->GetFunction("gausfit")->GetParError(2))/(waveformHist->GetFunction("gausfit")->GetParameter(2));

  }

  // For checking fit range
  /*
  for (int k = 12; k < NUMBER_DRIFT_BINS; k++) {
      if (sigmaVals[k]!=0) {
        sigmaVals[k] = 0;
        sigmaValsErrs[k] = 0;
      }
      if (driftTimes[k]!=-1) {
        driftTimes[k] = -1;
        driftTimesErrs[k] = -1;
      }
  }
  */

  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
  gStyle->SetTextFont(22);
  TPad *topPad = new TPad("topPad", "", 0.005, 0.58, 0.995, 0.995);
  TPad *midPad = new TPad("midPad", "", 0.005, 0.33, 0.995, 0.57);
  TPad *botPad = new TPad("botPad", "", 0.005, 0.005, 0.995, 0.33);
  topPad->SetBottomMargin(0.01);
  topPad->SetLeftMargin(0.15);
  midPad->SetLeftMargin(0.15);
  midPad->SetTopMargin(0.07);
  midPad->SetBottomMargin(0.05);
  midPad->SetGridy(1);
  botPad->SetTopMargin(0.07);
  botPad->SetLeftMargin(0.15);
  botPad->SetBottomMargin(0.35);
  botPad->SetLineStyle(2);
  topPad->Draw();
  midPad->Draw();
  botPad->Draw();

  topPad->cd();

  // Diffusion plot
  TGraphErrors *gr1 = new TGraphErrors(NUMBER_DRIFT_BINS,
      driftTimes, sigmaVals, 
      driftTimesErrs, sigmaValsErrs
  );
  gr1->SetTitle("");
  //gr1->GetXaxis()->SetTitle("Drift Time (#mus)");
  gr1->GetXaxis()->SetTitleSize(titleSize);
  gr1->GetXaxis()->SetLabelSize(labelSize);
  gr1->GetYaxis()->SetTitle("#sigma_{t}^{2} (#mus)");
  gr1->GetYaxis()->SetTitleSize(titleSize);
  gr1->GetYaxis()->SetTitleOffset(0.5);
  gr1->GetYaxis()->SetLabelSize(labelSize);
  gr1->Draw("ap");
  gr1->SetMarkerStyle(8);
  gr1->SetMarkerSize(0.8);
  //gr1->GetYaxis()->SetRangeUser(0.001, 10.);
  gr1->GetYaxis()->SetRangeUser(0.001, 10.);
  gr1->GetXaxis()->SetRangeUser(minTime, maxTime);

  // Linear fit to diffusion plot
  TF1* polFit = new TF1("polfit", "pol1");
  //gr1->Fit("polfit", "", "", minTime, 1450);
  gr1->Fit("polfit", "", "", minTime, maxTime);
  gr1->GetFunction("polfit")->SetLineColor(kRed);
  gr1->GetFunction("polfit")->Draw("same");

  // Get diffusion value from slope of linear fit
  double diffusionValue = polFit->GetParameter(1)*DRIFT_VELOCITY*DRIFT_VELOCITY* 1000000/2;
  std::cout << "Slope = " << polFit->GetParameter(1) << " +/- " << polFit->GetParError(1) << " microseconds" << std::endl;
  //std::cout << "Slope * v^2 = " << polFit->GetParameter(1)*std::pow(DRIFT_VELOCITY, 2) << std::endl; 
  double diffusionValueErr = polFit->GetParError(1)*std::pow(DRIFT_VELOCITY, 2) * 1000000/2;
  cout << "Diffusion value: " << diffusionValue << " +/- " << diffusionValueErr << endl; 

  TPaveText *pt = new TPaveText(0.16, 0.55, 0.6, 0.87, "NDC");
  pt->SetTextAlign(12);
  pt->SetFillStyle(0);
  pt->SetTextSize(textSize);
  pt->SetBorderSize(0);
  TString diffTextMeas = Form("Measured D_{L}: %0.2f +/- %0.2f cm^{2}/s", diffusionValue, diffusionValueErr);
  TString diffTextInp = Form("Input D_{L}: 6.40 cm^{2}/s");
  TString chi2 = Form("#chi^{2}/NDF: %0.2f/%i = %0.2f", polFit->GetChisquare(), polFit->GetNDF(), polFit->GetChisquare()/polFit->GetNDF() );
  TString sigma0 = Form("Measured #sigma_{0}^{2}: %0.2f +/- %0.2f #mus^{2}", polFit->Eval(400), polFit->GetParError(0) );
  pt->AddText(diffTextMeas);
  pt->AddText(diffTextInp);
  pt->AddText(chi2);
  pt->AddText(sigma0);
  pt->Draw("same");

  midPad->cd();
  double sigmaValsBottom[NUMBER_DRIFT_BINS]; 
  double sigmaValsErrsBottom[NUMBER_DRIFT_BINS];
  TF1* f2 = (TF1*)gr1->GetFunction("polfit");

  for (int i = 0; i < NUMBER_DRIFT_BINS; i++){

    sigmaValsBottom[i] = (sigmaVals[i] - (f2->Eval(driftTimes[i])))/f2->Eval(driftTimes[i]);
    sigmaValsErrsBottom[i] = (double)sigmaValsErrs[i]/f2->Eval(driftTimes[i]);
    //cout << "Sigma errs " << i << ": " << sigmaValsErrs[i] << std::endl;
    //cout << "Sigma errs bottom " << i << ": " << sigmaValsErrsBottom[i] << std::endl;

  }

  TGraphErrors* gr2 = new TGraphErrors(NUMBER_DRIFT_BINS, 
                                       driftTimes, sigmaValsBottom, 
                                       driftTimesErrs, sigmaValsErrsBottom);
  gr2->GetYaxis()->SetRangeUser(-0.05, 0.05);
  gr2->SetMarkerStyle(8);
  gr2->SetMarkerSize(0.8);
  gr2->SetTitle("");
  gr2->GetYaxis()->SetTitle("(#sigma_{t}^{2} - Fit)/Fit");
  gr2->GetXaxis()->SetLimits(minTime, maxTime);
  gr2->GetXaxis()->SetTitleSize(0);
  gr2->GetXaxis()->SetLabelSize(0);
  gr2->GetYaxis()->SetNdivisions(505);
  gr2->GetYaxis()->SetTitleOffset(0.35);
  gr2->GetYaxis()->SetTitleSize(titleSize*1.3);
  gr2->GetYaxis()->SetLabelSize(labelSize*1.3);

  gr2->Draw("ap");

  TF1* f3 = new TF1("f3", "0.", minTime, maxTime);
  f3->Draw("same");

  botPad->cd();
  gStyle->SetOptStat(0);
  histoNWvfms->GetXaxis()->SetTitle("Drift Time (#mus)");
  histoNWvfms->SetFillColor(kBlack);
  histoNWvfms->SetLineColor(kBlack);
  histoNWvfms->GetXaxis()->SetTitleSize(titleSize);
  histoNWvfms->GetXaxis()->SetLabelSize(labelSize);
  histoNWvfms->GetXaxis()->SetNdivisions(505);
  histoNWvfms->GetYaxis()->SetTitle("No. Waveforms");
  histoNWvfms->GetYaxis()->SetMaxDigits(3); 
  histoNWvfms->GetYaxis()->SetNdivisions(304); 
  histoNWvfms->GetYaxis()->SetTitleOffset(0.5);
  histoNWvfms->SetMinimum(0);
  histoNWvfms->GetYaxis()->SetTitleSize(titleSize);
  histoNWvfms->GetYaxis()->SetLabelSize(labelSize);
  histoNWvfms->GetXaxis()->SetLimits(minTime, maxTime);

  histoNWvfms->Draw();

  c1->cd();
  TString output_plot_dir = "DiffusionPlots/";
  c1->SaveAs(output_plot_dir+"DiffusionPlot.png", "PNG");
  c1->SaveAs(output_plot_dir+"DiffusionPlot.pdf", "PDF");
  delete c1;

  fOutput->Close();
  fInput->Close();

}

int main(int argv, char** argc){

  TString* inputFileName = new TString(argc[1]);
  if (!inputFileName) {
    cout << "No input file specified" << endl;
    return 1;
  }

  makePlot(inputFileName);
  return 0;

}
