#include "WaveformFunctions.h"

namespace diffmod {

  bool WaveformFunctions::passesHitSelection(art::Ptr< recob::Hit > hit, 
                                             double hitGOFCut, 
                                             int hitMultiplicityCut){

    if (hit->Multiplicity() <= hitMultiplicityCut && 
        hit->GoodnessOfFit() < hitGOFCut) {
      return true;
    }
    else {
      return false;
    }
  }

  double WaveformFunctions::convertXToTicks(double xPosition, int wvfmDriftStartTick, int wvfmDriftSize, double xWidth){

    //    double tick = WAVEFORM_START_TICK + 17.9427*xPosition;
    double tick = wvfmDriftStartTick + (((double)wvfmDriftSize/xWidth)*xPosition);
    return tick;

  }

  double WaveformFunctions::convertTicksToX(int tick, int wvfmDriftStartTick, int wvfmDriftSize, double xWidth){

    //double xPos = (tick - 812)/17.9427;
    double xPos = (tick - wvfmDriftStartTick)/(wvfmDriftSize/xWidth);
    return xPos;

  }

  TH1D* WaveformFunctions::applyGlobalBaselineCorrection(TH1D* h_rawD, TH1D* h_rawDCorrected){

    // first loop over the bins and find the beginning and end of the ROI
    std::pair<int, int> firstLastBinROI = std::make_pair(-1, -1);
    std::pair<int, int> clippedSignalBins = std::make_pair(-1,-1);
    double firstBinContent = h_rawD->GetBinContent(1);
    double lastBinContent = h_rawD->GetBinContent(h_rawD->GetNbinsX()); 
    double binContent;

    for (int i = 1; i <= h_rawD->GetNbinsX(); i++) {

      binContent = h_rawD->GetBinContent(i);
      if (binContent != firstBinContent){

        firstLastBinROI.first = i+5;
        break;
      }

    }

    for (int i = h_rawD->GetNbinsX(); i >=1; i--){ 

      binContent = h_rawD->GetBinContent(i);
      if (binContent != lastBinContent) {

        firstLastBinROI.second = i-5;
        break;
      }
    }

    // clip out signal region;
    clippedSignalBins.first = h_rawD->GetMaximumBin()-20;
    clippedSignalBins.second = h_rawD->GetMaximumBin()+20;

    // loop over bins and find baseline within ROI (outside of signal)
    double cumulativeSum = 0;
    int tickCounter = 0;
    for (int i = firstLastBinROI.first; i <= firstLastBinROI.second; i++){

      if (i >= clippedSignalBins.first && i <= clippedSignalBins.second) continue;

      cumulativeSum = cumulativeSum + h_rawD->GetBinContent(i);
      tickCounter++;

    }

    double baseline  = 0;
    if (tickCounter !=0) baseline = cumulativeSum/tickCounter;
    //std::cout << "[WAVEFORMFUNCS] Baseline = " << baseline << std::endl;

    // correct baseline and return histogram
    for (int i = 0; i <= h_rawDCorrected->GetNbinsX(); i++){ 
      h_rawDCorrected->SetBinContent(i, h_rawD->GetBinContent(i) - baseline);
    }

    return h_rawDCorrected;

  }

    /*
       TH1D* WaveformFunctions::rebin(TH1D* h_rawD, TH1D* h_rebinned, int NUMBER_BINS_PER_BIN, int wvfmDriftSize, int NUMBER_DRIFT_BINS){
  
  //
  // rebin function
  //
  
  // no standard way to split single bin in to multiple bins in root
  // used to change every 1-tick bin in to a 1/NUMBER_BINS_PER_BIN-tick bin
  // note that the ADC value is not scaled by 1/NUMBER_BINS_PER_BIN
  
  h_rebinned->SetBins(wvfmDriftSize/NUMBER_DRIFT_BINS * NUMBER_BINS_PER_BIN, h_rawD->GetXaxis()->GetXmin(), h_rawD->GetXaxis()->GetXmax());
  
  for (int i=1; i <= h_rawD->GetNbinsX(); i++){
  
  for (int j = 1; j <= NUMBER_BINS_PER_BIN; j++){
  
  h_rebinned->SetBinContent(((i-1)*NUMBER_BINS_PER_BIN)+j, h_rawD->GetBinContent(i));
  
  }
  
  }
  
  return h_rebinned;
  }
  */
  
  double WaveformFunctions::findXCorrection(TH1D* summedWaveform, TH1D* h, int nTicksPerBin, double fit_mean){
  
    int centerBin = -1;
    if (summedWaveform->GetMaximum() == 0){
      centerBin = nTicksPerBin/2;
    }
    else {
      // First element of vector returned by getSigma gives the fit_mean
      centerBin = summedWaveform->FindBin(this->getSigma(summedWaveform).at(0));
    }
  
    int distanceToBinCenter = 0;
    int testXCorr;
    double rms2 = 1000;
    for (int i = -5; i <= 5; i++){
      TH1D* h_summedClone = (TH1D*)summedWaveform->Clone("h_summedClone");
      TH1D* h_clone       = (TH1D*)             h->Clone("h_clone");
  
      testXCorr = centerBin - h->FindBin(fit_mean) + i;
  
      for (int ntick = 1; ntick <= h->GetNbinsX(); ntick++){
  
        h_clone->SetBinContent(ntick, h->GetBinContent(ntick + testXCorr));
  
      }

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
  
  std::vector<double> WaveformFunctions::getSigma(TH1D* h_rawDCorrected){
  
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
  
  double WaveformFunctions::getMedian(TH1D* h){
  
    double quantile = 0.5;
    double median; 
  
    if(h->Integral() == 0) {
      median = 0.0;
      return median;
    }
  
    h->GetQuantiles(1,&median,&quantile);                                
  
    return median;
  
  }
  
  double WaveformFunctions::getRms2(TH1D* h){
  
    double mean = 0;
    double mean2 = 0;
    double totalCharge = 0;
    double rms2 = 0;
    double threshold = 10.0 * h->GetMaximum()/100.0;
  
    for (int j = 1; j < h->GetNbinsX(); j++) {
  
      h->SetBinContent(j, h->GetBinContent(j) - threshold);
  
      if (h->GetBinContent(j) >= 0){
  
        totalCharge = totalCharge + h->GetBinContent(j);
        mean = mean + (double)h->GetBinCenter(j) * 0.5 * (double)h->GetBinContent(j);
        mean2 = mean2 + (double)h->GetBinCenter(j) * 0.5 * (double)h->GetBinCenter(j) *0.5 * (double)h->GetBinContent(j);
      }
  
      h->SetBinContent(j, h->GetBinContent(j) + threshold);
  
    }
  
    if (totalCharge !=0){
  
      mean = mean/totalCharge;
      mean2 = mean2/totalCharge;
      rms2 = sqrt(mean2 - mean*mean) * sqrt(mean2-mean*mean);
  
    }
  
    return rms2;
  
  }

}
