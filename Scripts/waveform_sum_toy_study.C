#include "StylePlots.h"

double findXCorrection(TH1D*, TH1D*, int);
double getRms2(TH1D*);
std::vector<double> getSigma(TH1D*);

void waveform_sum_toy_study() {
  SetGenericStyle();

  TFile *fout = new TFile("waveform_sum_toy_study_plots.root", "WRITE");

  // 25 drift time bins over 4600 tick readout
  gStyle->SetLegendTextSize(0.04);
  int  nbins     = 184;
  bool useAnode  = false; // Use cathode-like mean, sigma if false

  // Gaus and Sigma estimated using summed_waveform->GetMean() and summed_waveform->GetStdDev()
  Double_t gausMean, gausSigma; 
  if (useAnode) {
    gausMean  = 891.45; 
    gausSigma = 1.41858;
  }
  else {
    gausMean  = 5123.36; 
    gausSigma = 3.80828;
  }
  Double_t rangeLow  = gausMean - nbins/2, rangeHigh = gausMean + nbins/2;

  // Make Gaussian resembling waveforms
  TRandom3 r;
  TH1D *h_wvfm = new TH1D("h_wvfm", "", nbins, rangeLow, rangeHigh);
  for (int i = 0; i < 1e6; i++) {
    h_wvfm->Fill(r.Gaus(gausMean, gausSigma));
  }
  std::cout << "Start mean: "   << h_wvfm->GetMean() << std::endl;
  std::cout << "Start stddev: " << h_wvfm->GetStdDev() << std::endl;

  // Initialize clone (to be added) and sum (running summed waveform)
  TH1D *h_shifted     = (TH1D*)h_wvfm->Clone("h_shifted");
  TH1D *h_sum         = (TH1D*)h_wvfm->Clone("h_sum");
  TH1D *h_sum_noShift = (TH1D*)h_wvfm->Clone("h_sum_noShift");

  // Some handy stuff
  int    numAdditions = 1000;
  double xCorr        = 0.;
  double shift        = 0.;
  double maxShift     = 0.5;
  TRandom3 r_shift;

  //TCanvas *c_test = new TCanvas("c_test", "c_test");
  for (int i = 0; i < numAdditions; i++) {

    // Generate Gaussian with shifted mean, relative to h_wvfm
    h_shifted->Reset();
    shift = r_shift.Uniform(-maxShift, maxShift);
    std::cout << "Shift val: " << shift << std::endl;
    for (int i = 0; i < 1e6; i++) {
      h_shifted->Fill(r.Gaus(gausMean+shift, gausSigma));
    }

    std::cout << "------Add waveform no. " << i << "---------" << std::endl;
    std::cout << "Running summed mean: "   << h_sum->GetMean()   << std::endl;
    std::cout << "Running summed stddev: " << h_sum->GetStdDev() << std::endl;

    std::cout << "Shifted mean PRE-CORR: "   << h_shifted->GetMean()   << std::endl;
    std::cout << "Shifted stddev PRE-CORR: " << h_shifted->GetStdDev() << std::endl;

    xCorr = findXCorrection(h_sum, h_shifted, nbins);
    std::cout << "xCorr: " << xCorr << std::endl;

    // Apply x correction
    TH1D *h_corr = new TH1D("h_corr" , "", nbins, rangeLow, rangeHigh);
    for (int ntick = 1; ntick <= h_shifted->GetNbinsX(); ntick++) {
      h_corr->SetBinContent(ntick, h_shifted->GetBinContent(ntick - xCorr));
    }

    std::cout << "x-corrected mean: "   << h_corr->GetMean()   << std::endl;
    std::cout << "x-corrected stddev: " << h_corr->GetStdDev() << std::endl;

    // Clone summed waveform pre-addition for validation plot
    TH1D *h_temp = (TH1D*)h_sum->Clone();
    h_sum->Add(h_corr);

    // Add h_wvfm to itself for comparison
    h_sum_noShift->Add(h_wvfm);

    /*
    c_test->cd();
    h_sum->GetXaxis()->SetRangeUser(gausMean-20, gausMean+20);
    h_sum->GetXaxis()->SetTitle("Time (ticks)"); 
    h_sum->SetLineColor(kGreen+2);
    h_sum->Draw();

    h_shifted->SetLineColor(kBlack);
    h_shifted->Draw("same");

    h_corr->SetLineColor(kRed+1);
    h_corr->Draw("same");

    h_temp->SetLineColor(kAzure+1);
    h_temp->Draw("same");

    TLegend *l_test = new TLegend(0.55, 0.7, 0.75, 0.9);
    l_test->AddEntry(h_temp    , "Initial Wvfm"    , "l");
    l_test->AddEntry(h_shifted , "Uncorr. Added Wvfm", "l");
    l_test->AddEntry(h_corr    , "Corr. Added Wvfm", "l");
    l_test->AddEntry(h_sum     , "Summed Wvfm"  , "l");
    l_test->Draw("same");

    TString outHistName = Form("test_bin_%i.pdf", i);
    c_test->SaveAs(outHistName, "PDF");
    */

    delete h_corr;

  }

  std::cout << "Summed mean: "   << h_sum->GetMean()   << std::endl;
  std::cout << "Summed stddev: " << h_sum->GetStdDev() << std::endl;

  TFile *f_comp = new TFile("/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionAna/diffmod_run3_crt_Aug2020_newFV_bugFix.root", "READ");
  TH1D *h_sum_data = (TH1D*)f_comp->Get("DiffusionModule/plane2/summed_waveform_bin_0_plane2");
  TCanvas *c = new TCanvas("c", "c", 500, 500);
  c->cd();
  //c->SetLogy();

  h_sum_noShift->GetXaxis()->SetRangeUser(gausMean-20, gausMean+20);
  h_sum_noShift->GetXaxis()->SetTitle("Time (ticks)");
  h_sum_noShift->GetXaxis()->CenterTitle();
  h_sum_noShift->SetLineColor(kPTVibrantCyan);
  h_sum_noShift->GetYaxis()->SetRangeUser(0, h_sum_noShift->GetMaximum()*1.25);
  h_sum_noShift->Draw();
  //h_wvfm->DrawNormalized();

  h_sum->SetLineColor(kPTVibrantMagenta);
  h_sum->Draw("same");
  //h_sum->DrawNormalized("same");

  TLegend *l_comp = new TLegend(0.55, 0.75, 0.75, 0.85);
  l_comp->SetLineWidth(0);
  l_comp->SetFillStyle(0);
  l_comp->AddEntry(h_sum_noShift, "Un-Shifted Sum", "l");
  l_comp->AddEntry(h_sum        , "Shifted Sum"   , "l");
  l_comp->Draw("same");

  std::cout << "Start Mean: "   << h_sum_noShift->GetMean()   << std::endl;
  std::cout << "Summed Mean: "     << h_sum     ->GetMean()   << std::endl;
  std::cout << "Start StdDev: " << h_sum_noShift->GetStdDev() << std::endl;
  std::cout << "Summed StdDev: "   << h_sum     ->GetStdDev() << std::endl;

  c->SaveAs("toystudy.pdf");

  /*
  h_wvfm->SetLineColor(kAzure+1);
  h_wvfm->DrawNormalized("same");

  h_shifted->SetLineColor(kRed+1);
  h_shifted->Draw("same");
  */

}

// Function copies from WaveformAlgorithms.cxx (with minor edits)
double getRms2(TH1D* h){

  double mean = 0;
  double mean2 = 0;
  double totalCharge = 0;
  double rms2 = 0;
  double threshold = 10.0 * h->GetMaximum()/100.0;

  for (int j = 1; j < h->GetNbinsX(); j++) {

    h->SetBinContent(j, h->GetBinContent(j) - threshold);

    if (h->GetBinContent(j) >= 0){

      totalCharge = totalCharge + h->GetBinContent(j);
      // 0.5 because ticks -> us(?)
      mean  = mean  + (double)h->GetBinCenter(j) * 0.5 * (double)h->GetBinContent(j);
      mean2 = mean2 + (double)h->GetBinCenter(j) * 0.5 * (double)h->GetBinCenter(j)  * 0.5 * (double)h->GetBinContent(j);
    }

    h->SetBinContent(j, h->GetBinContent(j) + threshold);

  }

  if (totalCharge != 0){

    mean  = mean /totalCharge;
    mean2 = mean2/totalCharge;
    rms2 = sqrt(mean2 - mean*mean) * sqrt(mean2-mean*mean);

  }

  return rms2;

}

std::vector<double> getSigma(TH1D* h_rawDCorrected){

  std::vector<double> returnVector;
  if (h_rawDCorrected->Integral() == 0 || h_rawDCorrected->GetMaximum() <= 0 ) {
    returnVector = {-1.0, -1.0, -1.0};
    return returnVector;
  }

  // loop over histogram, find first point below 0.1 ADC, fit only between these regions
  int maxBin = h_rawDCorrected->GetMaximumBin();
  double maxVal = h_rawDCorrected->GetMaximum();
  double fitLowerLimit = 0;
  double fitHigherLimit = 0;
  double percentageOfHeight = 10.;
  double cutOff = maxVal* percentageOfHeight/100.;

  for (int i = 0; i < 200; i++){

    double binValueLow = h_rawDCorrected->GetBinContent(maxBin-i);
    if (binValueLow < cutOff){

      fitLowerLimit = h_rawDCorrected->GetBinCenter(maxBin-i); 
      break;
    }

  }

  for (int i = 0; i < 200; i++){

    double binValueHigh = h_rawDCorrected->GetBinContent(maxBin+i);
    if (binValueHigh < cutOff){

      fitHigherLimit = h_rawDCorrected->GetBinCenter(maxBin+i); 
      break;
    }

  }


  h_rawDCorrected->Fit("gaus", "0q", "", fitLowerLimit, fitHigherLimit);

  double sigma;
  double chisq;
  double mean;
  //double mean_err;

  if (h_rawDCorrected->GetFunction("gaus")){

    mean =  h_rawDCorrected->GetFunction("gaus")->GetParameter(1);
    sigma = h_rawDCorrected->GetFunction("gaus")->GetParameter(2);
    chisq = h_rawDCorrected->GetFunction("gaus")->GetChisquare();

  }
  else {
    mean = -1;
    sigma = -1;
    chisq = -1;
  }

  returnVector = {mean, sigma, chisq};

  return returnVector;
}

double findXCorrection(TH1D* summedWaveform, TH1D* h, int nTicksPerBin){

  int centerBin = -1;
  if (summedWaveform->GetMaximum() == 0){
    centerBin = nTicksPerBin/2;
  }
  else {
    // First element of vector returned by getSigma gives the Gaussian fit mean
    std::vector<double> tmp_v = getSigma(summedWaveform);
    centerBin = summedWaveform->FindBin(tmp_v.at(0));
  }

  int distanceToBinCenter = 0;
  int testXCorr = 0;
  double rms2 = 1000;
  double fit_mean = getSigma(h).at(0);

  for (int i = -5; i <= 5; i++){
    TH1D* h_summedClone = (TH1D*)summedWaveform->Clone("h_summedClone");
    TH1D* h_clone       = (TH1D*)             h->Clone("h_clone");

    testXCorr = centerBin - h->FindBin(fit_mean) + i;

    for (int ntick = 1; ntick <= h->GetNbinsX(); ntick++){
      h_clone->SetBinContent(ntick, h->GetBinContent(ntick - testXCorr));
      //h_clone->SetBinContent(ntick, h->GetBinContent(ntick + testXCorr));
    }

    // Treat errors correctly...? Maybe?
    h_summedClone->Sumw2();
    h_summedClone->Add(h_clone);

    double rms2Test = getRms2(h_summedClone);
    if (rms2Test < rms2){

      rms2 = rms2Test;
      distanceToBinCenter = testXCorr;

    }

    h_clone->Delete();
    h_summedClone->Delete();

  }

  return distanceToBinCenter;
}







