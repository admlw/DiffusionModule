#include "StylePlots.h"

double findXCorrection(TH1D*, TH1D*, int);
double getRms2(TH1D*);
std::vector<double> getSigma(TH1D*);
std::pair<float, float> getFitRange(TH1D*);

void waveform_sum_toy_study() {

  SetGenericStyle();

  TFile *fout = new TFile("waveform_sum_toy_study_plots.root", "WRITE");

  // 25 drift time bins over 4600 tick readout
  gStyle->SetOptStat(0);
  gStyle->SetLegendTextSize(0.05);
  int  nbins     = 184;
  bool useAnode  = false; // If false, use cathode-like mean and sigma 

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

  TF1 *gaus_start = new TF1("gaus_start", "gaus");
  std::pair<float, float> fitRange_start = getFitRange(h_wvfm);
  h_wvfm->Fit(gaus_start, "0q", "", fitRange_start.first, fitRange_start.second);
  
  std::cout << "Start Gaus mean: "     << gaus_start->GetParameter(1) << std::endl;
  std::cout << "Start Gaus stddev: "   << gaus_start->GetParameter(2) << std::endl;
  std::cout << "start Gaus Chi2: "     << gaus_start->GetChisquare() << std::endl;
  std::cout << "start Gaus NDF: "      << gaus_start->GetNDF() << std::endl;
  std::cout << "start Gaus Chi2/NDF: " << gaus_start->GetChisquare()/gaus_start->GetNDF() << std::endl;

  // Initialize clone (to be added) and sum (running summed waveform)
  TH1D *h_shifted     = (TH1D*)h_wvfm->Clone("h_shifted");
  TH1D *h_sum         = (TH1D*)h_wvfm->Clone("h_sum");
  TH1D *h_sum_noShift = (TH1D*)h_wvfm->Clone("h_sum_noShift");

  // Some handy stuff
  int const numAdditions = 1000;
  double    xCorr        = 0.;
  double    shift        = 0.;
  double    maxShift     = 0.5;
  TRandom3 r_shift;
  double chisqNdfVals[numAdditions] = {0.};
  double xVals       [numAdditions] = {0.};

  //TCanvas *c_test = new TCanvas("c_test", "c_test");
  for (int i = 0; i < numAdditions; i++) {

    double chisqNdf = 10;
    double chisq    = 10;
    double ndf      = 1;
    double sigma    = 100;
    double sigmaErr = 100;

    // Generate Gaussian with shifted mean, relative to h_wvfm
    h_shifted->Reset();
    shift = r_shift.Uniform(-maxShift, maxShift);
    //std::cout << "Shift val: " << shift << std::endl;
    for (int i = 0; i < 1e6; i++) {
      h_shifted->Fill(r.Gaus(gausMean+shift, gausSigma));
    }

    /*
    std::cout << "------Add waveform no. " << i << "---------" << std::endl;
    std::cout << "Running summed mean: "   << h_sum->GetMean()   << std::endl;
    std::cout << "Running summed stddev: " << h_sum->GetStdDev() << std::endl;

    std::cout << "Shifted mean PRE-CORR: "   << h_shifted->GetMean()   << std::endl;
    std::cout << "Shifted stddev PRE-CORR: " << h_shifted->GetStdDev() << std::endl;
    */

    xCorr = findXCorrection(h_sum, h_shifted, nbins);
    //std::cout << "xCorr: " << xCorr << std::endl;

    // Apply x correction
    TH1D *h_corr = new TH1D("h_corr" , "", nbins, rangeLow, rangeHigh);
    for (int ntick = 1; ntick <= h_shifted->GetNbinsX(); ntick++) {
      h_corr->SetBinContent(ntick, h_shifted->GetBinContent(ntick - xCorr));
    }

    //std::cout << "x-corrected mean: "   << h_corr->GetMean()   << std::endl;
    //std::cout << "x-corrected stddev: " << h_corr->GetStdDev() << std::endl;

    // Clone summed waveform pre-addition for validation plot
    TH1D *h_temp = (TH1D*)h_sum->Clone();
    h_sum->Add(h_corr);

    // Add h_wvfm to itself for comparison
    h_sum_noShift->Add(h_wvfm);

    // Gaussian fit to running sum
    // See if chi^2 increases with each addition
    TF1 *gaus_tmp = new TF1("gaus_tmp", "gaus");
    //h_sum->Fit(gaus_tmp, "q", "", 
    //           gausMean-20., gausMean+20.);

    chisq    = gaus_tmp->GetChisquare();
    ndf      = gaus_tmp->GetNDF();
    sigma    = gaus_tmp->GetParameter(2);
    sigmaErr = gaus_tmp->GetParError(2);
    chisqNdf = chisq/ndf;

    chisqNdfVals[i] = chisqNdf;
    xVals[i]        = i;

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
    delete gaus_tmp;

  }

  std::cout << "Summed mean: "   << h_sum->GetMean()   << std::endl;
  std::cout << "Summed stddev: " << h_sum->GetStdDev() << std::endl;

  TFile *f_comp = new TFile("/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionAna/diffmod_run3_crt_Aug2020_newFV_bugFix.root", "READ");
  TH1D *h_sum_data = (TH1D*)f_comp->Get("DiffusionModule/plane2/summed_waveform_bin_0_plane2");

  TCanvas *c = new TCanvas;
  c->cd();
  //c->SetLogy();

  h_sum_noShift->GetXaxis()->SetRangeUser(gausMean-20, gausMean+20);
  //h_sum_noShift->GetYaxis()->SetRangeUser(0., 0.5);
  h_sum_noShift->GetXaxis()->SetTitle("Time (ticks)");
  h_sum_noShift->GetXaxis()->SetTitleSize(0.05);
  h_sum_noShift->SetLineColor(kPTVibrantMagenta);
  //h_sum_noShift->SetLineStyle(7);
  h_sum_noShift->SetLineWidth(2);
  h_sum_noShift->Draw();
  //h_sum_noShift->DrawNormalized("l");

  h_sum->SetLineColor(kPTDarkBlue);
  //h_sum->SetLineStyle(9);
  h_sum->SetLineWidth(2);
  h_sum->Draw("same");
  //h_sum->DrawNormalized("same l");

  TLegend *l_comp = new TLegend(0.58, 0.65, 0.85, 0.85);
  l_comp->AddEntry(h_sum_noShift, "Un-Shifted Sum", "l");
  l_comp->AddEntry(h_sum        , "Shifted Sum"   , "l");
  l_comp->Draw("same");

  std::cout << "Start Mean: "   << h_sum_noShift->GetMean()   << std::endl;
  std::cout << "Summed Mean: "     << h_sum     ->GetMean()   << std::endl;
  std::cout << "Start StdDev: " << h_sum_noShift->GetStdDev() << std::endl;
  std::cout << "Summed StdDev: "   << h_sum     ->GetStdDev() << std::endl;

  std::pair<float, float> fitRange_shift = getFitRange(h_sum);
  TF1 *gaus_final = new TF1("gaus_final", "gaus");
  //h_sum->Fit(gaus_final, "q", "", gausMean-20., gausMean+20.);
  h_sum->Fit(gaus_final, "0q", "", fitRange_shift.first, fitRange_shift.second);
  gaus_final->SetLineColor(kPTRed);
  gaus_final->SetLineStyle(9);
  //gaus_final->Draw("same");
  
  std::cout << "final Gaus mean: "     << gaus_final->GetParameter(1) << std::endl;
  std::cout << "final Gaus stddev: "   << gaus_final->GetParameter(2) << std::endl;
  std::cout << "final Gaus Chi2: "     << gaus_final->GetChisquare() << std::endl;
  std::cout << "final Gaus NDF: "      << gaus_final->GetNDF() << std::endl;
  std::cout << "final Gaus Chi2/NDF: " << gaus_final->GetChisquare()/gaus_final->GetNDF() << std::endl;

  std::pair<float, float> fitRange_noShift = getFitRange(h_sum_noShift);
  TF1 *gaus_final_noShift = new TF1("gaus_final_noShift", "gaus");
  //h_sum_noShift->Fit(gaus_final_noShift, "q", "", gausMean-20., gausMean+20.);
  h_sum_noShift->Fit(gaus_final_noShift, "0q", "", fitRange_noShift.first, fitRange_noShift.second);
  gaus_final_noShift->SetLineColor(kPTLightBlue);
  gaus_final_noShift->SetLineStyle(5);
  //gaus_final_noShift->Draw("same");

  std::cout << "final_noShift Gaus mean: "     << gaus_final_noShift->GetParameter(1) << std::endl;
  std::cout << "final_noShift Gaus stddev: "   << gaus_final_noShift->GetParameter(2) << std::endl;
  std::cout << "final_noShift Gaus Chi2: "     << gaus_final_noShift->GetChisquare() << std::endl;
  std::cout << "final_noShift Gaus NDF: "      << gaus_final_noShift->GetNDF() << std::endl;
  std::cout << "final_noShift Gaus Chi2/NDF: " << gaus_final_noShift->GetChisquare()/gaus_final_noShift->GetNDF() << std::endl;


  TString output_plot_name;
  if (useAnode) output_plot_name = "waveform_sum_anode.pdf";
  else          output_plot_name = "waveform_sum_cathode.pdf";
  c->SaveAs(output_plot_name, "PDF");
  
  /*
  TCanvas *c_chi2 = new TCanvas("c_chi2", "", 800, 600);
  TGraph *g_chi2 = new TGraph(numAdditions, xVals, chisqNdfVals);
  g_chi2->SetTitle("");
  g_chi2->GetXaxis()->SetTitle("Num Additions");
  g_chi2->GetYaxis()->SetTitle("Gaus. #chi^{2}/NDF");
  g_chi2->SetMarkerStyle(8);
  g_chi2->SetMarkerSize(0.5);
  g_chi2->Draw();

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

    //h_clone->Sumw2();
    //h_summedClone->Sumw2();
    h_summedClone->Add(h_clone);
    //h_summedClone->Sumw2();

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


std::pair<float, float> getFitRange(TH1D* wvfm){

  std::pair<float, float> fitRanges(0.,0.);

  double lowConv = wvfm->GetBinLowEdge(1);
  double highConv = wvfm->GetBinLowEdge(wvfm->GetNbinsX()+1);
  wvfm->GetXaxis()->SetLimits(lowConv, highConv);

  // Stop fit at 10% of maximum value 
  double fitLimit = wvfm->GetMaximum()*0.1;
  for (int i = wvfm->GetMaximumBin(); i > 0; i--) {
    if (wvfm->GetBinContent(i) < fitLimit) {
      fitRanges.first = wvfm->GetBinCenter(i);
      break;
    }
  }
  for (int i = wvfm->GetMaximumBin(); i < wvfm->GetNbinsX(); i++) {
    if (wvfm->GetBinContent(i) < fitLimit) {
      fitRanges.second = wvfm->GetBinCenter(i);
      break;
    }
  }

  return fitRanges;
}





