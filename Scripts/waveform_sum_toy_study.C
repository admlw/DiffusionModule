double findXCorrection(TH1D*, TH1D*, int, double);
double getRms2(TH1D*);
std::vector<double> getSigma(TH1D*);

void waveform_sum_toy_study() {

  // 25 drift time bins over 4600 tick readout
  Int_t nbins = 184;
  Int_t meanDrift = 900; // 900 ticks = near anode
  Int_t rangeLow = meanDrift - nbins/2, rangeHigh = meanDrift + nbins/2;
  // Gaus and Sigma estimated using summed_waveform->GetMean() and summed_waveform->GetStdDev()
  Double_t gausMean = 900.; // Near anode
  Double_t gausSigma = 1.41858;

  // Make Gaussian resembling waveforms
  TRandom3 r;
  TH1D *h_wvfm = new TH1D("h_wvfm", "", nbins, rangeLow, rangeHigh);
  for (int i = 0; i < 1e6; i++) {
    h_wvfm->Fill(r.Gaus(gausMean, gausSigma));
  }

}

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

double findXCorrection(TH1D* summedWaveform, TH1D* h, int nTicksPerBin, double fit_mean){

  int centerBin = -1;
  if (summedWaveform->GetMaximum() == 0){
    centerBin = nTicksPerBin/2;
  }
  else {
    // First element of vector returned by getSigma gives the fit_mean
    std::vector<double> tmp_v = getSigma(summedWaveform);
    centerBin = summedWaveform->FindBin(tmp_v.at(0));
  }

  int distanceToBinCenter = 0;
  int testXCorr;
  double rms2 = 1000;
  for (int i = -5; i <= 5; i++){
    TH1D* h_summedClone = (TH1D*)summedWaveform->Clone("h_summedClone");
    TH1D* h_clone = (TH1D*)h->Clone("h_clone");

    testXCorr = centerBin - h->FindBin(fit_mean) + i;

    for (int ntick = 1; ntick <= h->GetNbinsX(); ntick++){

      h_clone->SetBinContent(ntick, h->GetBinContent(ntick + testXCorr));

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
