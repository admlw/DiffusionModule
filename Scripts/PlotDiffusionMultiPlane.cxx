#include "PlotDiffusionMultiPlane.h"

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
  std::ostringstream out;
  out.precision(n);
  out << std::fixed << a_value;
  return out.str();
}

void drawPaveText(std::string plane, double diffV, double diffE, double chisq, double ndf, double sig0, double sig0err){

  std::string diffString =
    "D_{L}: " 
    + to_string_with_precision(diffV,2)
    + " +/- "
    + to_string_with_precision(diffE,2)
    + "cm^{2}/s";

  std::string chisqString =
    "#chi^{2}/ndf: "
    + to_string_with_precision(chisq,2) 
    + "/" 
    + to_string_with_precision(ndf  ,2);

  std::string sigma0String =
    "#sigma_{0}^{2}: "
    + to_string_with_precision(sig0     , 2) 
    + "+/-" 
    + to_string_with_precision(sig0err  , 2);
    + "#mus^{2}";

  //TString sigma0 = Form("Measured #sigma_{0}^{2}: %0.2f +/- %0.2f #mus^{2}", polFit->Eval(0), polFit->GetParError(0) );

  TPaveText* tpv = new TPaveText(0.12, 0.60, 0.7, 0.90, "NDC");
  /*
=======
  TPaveText* tpv = new TPaveText(0.16, 0.69, 0.7, 0.89, "NDC");
>>>>>>> b6207e378f586eae25ad337694756eb941b67434
*/
  tpv->SetTextAlign(11);
  tpv->SetFillStyle(0);
  tpv->SetLineWidth(0);
  tpv->SetTextFont(43);
  tpv->SetTextSize(18);
  tpv->AddText(plane.c_str());
  tpv->AddText(diffString.c_str());
  tpv->AddText(chisqString.c_str());
  tpv->AddText(sigma0String.c_str());
  tpv->DrawClone("same");
}

void styleGraph(TGraph* h, float minx, float maxx, float miny, float maxy){
  int fontSize = 20;

  h->GetXaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleFont(43);
  h->GetXaxis()->SetTitleSize(fontSize);
  h->GetYaxis()->SetTitleSize(fontSize);
  h->GetXaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(fontSize);
  h->GetYaxis()->SetLabelSize(fontSize);
  h->GetXaxis()->SetTitleOffset(3.5);
  h->GetYaxis()->SetTitleOffset(1.8);
  h->GetXaxis()->SetRangeUser(minx, maxx);
  h->GetYaxis()->SetRangeUser(miny, maxy);
  h->GetYaxis()->SetNdivisions(505);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.4);
}

std::pair<float, float> getFitRange(TH1D* wvfm){

  std::pair<float, float> fitRanges(0.,0.);

  // Convert tick histogram to microseconds, find fit range 
  double lowConv = wvfm->GetBinLowEdge(1)*0.5;
  double highConv = wvfm->GetBinLowEdge(wvfm->GetNbinsX()+1)*0.5;
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

void increaseError(TH1D* h){
  for (int i = 0; i < h->GetNbinsX(); i++){
    h->SetBinError(i+1, h->GetBinError(i+1)*1.02);
  }
}

std::pair<int, int> getBinXError(TH1D* driftHisto){

    std::pair<int, int> errLowHigh;
    double lowBinVal  = driftHisto->GetBinContent(1);
    double highBinVal = driftHisto->GetBinContent(driftHisto->GetNbinsX());
    double midBinVal  = driftHisto->GetBinContent((int)driftHisto->GetNbinsX()/2);

    // deal with U-shaped bins
    if (lowBinVal > midBinVal && highBinVal > midBinVal){
        errLowHigh.first  = driftHisto->GetBinLowEdge(1);
        errLowHigh.second = driftHisto->GetBinLowEdge(driftHisto->GetNbinsX());
     
        return errLowHigh;
    }
   
    // find max bin nearest mean
    double driftMean   = driftHisto->GetMean();
    double driftStdDev = driftHisto->GetStdDev();
    double maximum     = 0;
    double maximumbin  = 0;

    for (int i = 0; i < driftHisto->GetNbinsX(); i++){
       double binVal = driftHisto->GetBinContent(i);
       if (binVal > maximum){
           maximum    = binVal;
           maximumbin = i;
       }
       if (binVal == maximum){
           if (std::abs(driftMean - driftHisto->GetBinCenter(i)) < std::abs(driftMean - driftHisto->GetBinCenter(maximumbin))){
               maximum    = binVal;
               maximumbin = i;
           }
       }

    }
   
    int    driftMaxBin        = driftHisto->GetMaximumBin();
    double driftMaxBinContent = driftHisto->GetBinContent(driftMaxBin);
    double driftIntegral      = driftHisto->Integral();
    
    int lowBin  = driftMaxBin;
    int highBin = driftMaxBin;
     
    double fractionOfTotal = driftMaxBinContent/driftIntegral; 
      
    while (fractionOfTotal < 0.68){
       
        // find whether to step left or right
        lowBinVal  = driftHisto->GetBinContent(lowBin-1);
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
            lowBin  = lowBin-1;
            highBin = highBin+1;
        }
        else{
            std::cout << "---ERRORS ARE NOT GOOD" << std::endl;
            std::cout << "lowBinVal: " << lowBinVal 
                      << " highBinVal: " << highBinVal << std::endl;
            std::cout << "lowBin: " << lowBin 
                      << " highBin: " << highBin << std::endl;
            std::cout << "nBins: " << driftHisto->GetNbinsX() << std::endl;
        }
    }
  
    errLowHigh.first  = driftHisto->GetBinLowEdge(lowBin+1);
    errLowHigh.second = driftHisto->GetBinLowEdge(highBin-1);
  
    return errLowHigh;
}


void makePlot(std::string inputFileName){

  // Initialization 
  // TODO: check these numbers are ok!
  const int  waveformDriftStartTick = 800;
  const int  waveformDriftEndTick   = 5400;
  const int  waveformDriftSize      = waveformDriftEndTick
                                      -waveformDriftStartTick; // End tick (5400 - 800)
  const int  numWaveformCut         = 0;                     // Cut bins with fewer than this number of waveforms
  const int  numDriftBins           = 25;
  int        nuTicksPerBin          = waveformDriftSize 
                                     /numDriftBins;
  const int  minTime                = waveformDriftStartTick/2; // 400 microseconds
  const int  maxTime                = waveformDriftEndTick/2;   // 2700 microseconds
  const bool isData                 = true;
  const bool isMakeWaveformPlots    = false;
  double     driftVelocity;

  // For data measurement, use drift velocity at anode. For MC, use
  // drift velocity at nominal E-field of 273 V/cm. Why? Basically 
  // because that's what the simulation does, even though we think
  // the anode drift velocity is the better thing to use
  if (isData) driftVelocity = 0.10762;  // Anode drift velocity
  else        driftVelocity = 0.1098;   // Average drift velocity

  TFile *fInput = new TFile(inputFileName.c_str(), "READ");
  if (!fInput) {
    std::cout << "Bad file" << std::endl;
    return;
  }

  // from input file, determine which planes to use
  std::vector< bool > isUsePlane = {false, false, false};
  TDirectoryFile* topDir = (TDirectoryFile*)fInput->Get("DiffusionModule");
  if (topDir->Get("plane0"))
    isUsePlane.at(0) = true;
  if (topDir->Get("plane1"))
    isUsePlane.at(1) = true;
  if (topDir->Get("plane2"))
    isUsePlane.at(2) = true;

  // one each for the 3 planes
  double driftTimes      [3][numDriftBins];
  double driftTimesErrs  [3][numDriftBins]; 
  double sigmaSqrVals    [3][numDriftBins]; 
  double sigmaSqrValsErrs[3][numDriftBins];


  // define output file
  TFile *fOutput = new TFile(("fits_"+(inputFileName).substr((inputFileName).find_last_of("/")+1,
                                                             (inputFileName).length())).c_str(), 
                             "RECREATE");
  if (!fOutput) {
    std::cout << "Bad output file" << std::endl;
    return;
  }

  // pick up tree
  TTree *t = (TTree*)fInput->Get("DiffusionModule/difftree");
  if (!t) {
    std::cout << "Bad tree" << std::endl;
      return;
  }

  std::vector<TH1D*> nWvfmsVec;
	float waveformHistXLow;
	float waveformHistXHigh;
  for (int ip = 0; ip < isUsePlane.size(); ip++){
  
    if (isUsePlane.at(ip) == false) 
      continue;
    else std::cout << "Using Plane " << ip << std::endl;
    
    std::string planeName = "plane"+std::to_string(ip);
   
    // get some basic histograms
    TCanvas *c = new TCanvas(("c"+planeName).c_str());
    TH1D* hPeakTime = new TH1D(("hPeakTime"+planeName).c_str(),
                               ";hit peak time;",
                               1150, 
                               waveformDriftStartTick, 
                               waveformDriftEndTick);

    c->cd();

    std::string drawString = "hit_peak_time >> hPeakTime" + planeName;
    std::string cutString  = "hit_view == "+std::to_string(ip);
    t->Draw(drawString.c_str(), cutString.c_str());

    c->SaveAs(("hitPeakTime"+planeName+".png").c_str());

    std::string hNWvfmLoc = "DiffusionModule/"
                            + planeName
                            + "/h_nWvfmsInBin"
                            + planeName;

    std::cout << "printing name: " << hNWvfmLoc << std::endl;

    TH1D* hNWvfms = (TH1D*)fInput->Get(hNWvfmLoc.c_str());
    nWvfmsVec.push_back(hNWvfms);
    
    if (hNWvfms == nullptr)
      throw std::logic_error ("bad waveform histogram");

    // loop over the drift bins for this plane
    TH1D* waveformHist;
    for (int idb = 0; idb < numDriftBins; ++idb){
      std::cout 
        << "Looping plane " 
        << ip 
        << ", drift bin " 
        << idb
        << std::endl;

      if (hNWvfms->GetBinContent(idb+1) < numWaveformCut) {
        std::cout 
          << "Skipping, fewer than " 
          << numWaveformCut 
          << " waveforms" 
          << std::endl;
        
        sigmaSqrVals    [ip][idb] = 0;
        sigmaSqrValsErrs[ip][idb] = 0;
        driftTimes      [ip][idb] = -1;
        driftTimesErrs  [ip][idb] = -1;
        continue;
      }
      
      // Change fit range/binning
      /*
      for (int k = 5; k < 21; k++) {
          if (sigmaSqrVals[ip][k] != 0) {

              sigmaSqrVals[ip][k]     = 0;
              sigmaSqrValsErrs[ip][k] = 0;

          }
          if (driftTimes[ip][k] != -1) {

              driftTimes    [ip][k] = -1;
              driftTimesErrs[ip][k] = -1;

          }
      }
      */

      fOutput->cd();
      std::string histLoc = "DiffusionModule/"
                            + planeName
                            + "/summed_waveform_bin_"
                            + std::to_string(idb)
                            + "_" + planeName;

      waveformHist = (TH1D*)fInput->Get(histLoc.c_str());

      // Ensure histogram is filled
      if (waveformHist->Integral() == 0){
        std::cout << "bin " << idb << " is empty!" << std::endl;
        driftTimes  [ip][idb]  = -1; 
        sigmaSqrVals[ip][idb]  = 0;
        continue;
      }

      std::pair<float, float> fitRanges = getFitRange(waveformHist);

      TF1 *gausfit = new TF1("gausfit", "gaus");
      waveformHist->Fit(gausfit, "q", "", 
                        fitRanges.first, fitRanges.second);

      std::cout << "Gaus fit range " << fitRanges.first  << 
                   " to "            << fitRanges.second << std::endl;
      // Error inflation
      double chisqNdf = 10;
      double chisq    = 10;
      double ndf      = 1;
      double sigma    = 100;
      double sigmaErr = 100;
      // Chi2 inflation
      int iter = 0;
      while (chisqNdf > 1) {
        waveformHist->Fit(gausfit, "q", "", 
                          fitRanges.first, fitRanges.second);

        chisq    = gausfit->GetChisquare();
        ndf      = gausfit->GetNDF();
        sigma    = gausfit->GetParameter(2);
        sigmaErr = gausfit->GetParError(2);
        chisqNdf = chisq/ndf;
        increaseError(waveformHist);
        iter++;
      }
      iter=0;
      
      std::cout << "Num iter  = " << iter              << std::endl;
      std::cout << "chi^2     = " << chisq             << std::endl;
      std::cout << "NDF       = " << ndf               << std::endl;
      std::cout << "chi^2/NDF = " << (double)chisq/ndf << std::endl;
      std::cout << "sigma     = " << sigma             << std::endl;
      std::cout << "sigma err = " << sigmaErr          << std::endl;

      // this is to make plots for single waveforms
      if (isMakeWaveformPlots){
        TCanvas *c_test = new TCanvas("c_test", "", 750, 550);
        c_test->cd();
        gStyle->SetOptFit(1);
        waveformHist->Draw("e1");
        TString waveformHistName = Form("waveformHist_%i", idb);
        TString testDir          = "waveformHistPlots/";
        c_test->SaveAs(testDir+waveformHistName+planeName+".pdf");
        c_test->SaveAs(testDir+waveformHistName+planeName+".png");
      }

      hPeakTime->GetXaxis()->SetLimits(fitRanges.first, fitRanges.second);
  
      // get useful information from histogram
      waveformHistXLow  = waveformHist->GetBinLowEdge(1);
      waveformHistXHigh = 
        waveformHist->GetBinLowEdge(waveformHist->GetNbinsX()-1);

      std::cout << "Setting waveformHist bounds to " 
                << waveformHistXLow 
                << " and " 
                << waveformHistXHigh 
                << std::endl;
      hPeakTime->GetXaxis()->SetRangeUser(waveformHistXLow, 
                                          waveformHistXHigh);

      double binMean   = hPeakTime->GetMean();
      double binStdDev = hPeakTime->GetStdDev();
      int    N         = hPeakTime->GetEntries();

      // get drift time from truncated mean of sigma distribution
      // in each bin
      // Error is (1/(sqrt(N)) * 0.5/2 (half a tick width, if not using 
      // hit information)
      std::cout << "binMean: " << binMean << " minTime: " << minTime << std::endl;
      driftTimes    [ip][idb] = binMean - minTime;
      driftTimesErrs[ip][idb] = (1/sqrt(N) * (0.5/2));

      // get pulse width squared
      sigmaSqrVals    [ip][idb] = std::pow(sigma,2);
      sigmaSqrValsErrs[ip][idb] = sqrt(2) * sigmaSqrVals[ip][idb] * (sigmaErr/sigma);

    }
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 1000);
  c1->cd();
  gStyle->SetTextFont(43);
  TPad* uCan = new TPad("uCan", "", 0.005, 0.3 , 0.995, 0.995);
  TPad* uRat = new TPad("uRat", "", 0.005, 0.005 , 0.995, 0.3);
  uCan->SetTopMargin(0.07);
  uCan->SetBottomMargin(0.01);
	uCan->SetLeftMargin(0.15);
  uRat->SetTopMargin(0.02);
  uRat->SetBottomMargin(0.3);
	uRat->SetLeftMargin(0.15);
  uCan->Draw();
  uRat->Draw();

	TCanvas *c2 = new TCanvas("c2", "c2", 500, 500);
  TPad* vCan = new TPad("vCan", "", 0.005, 0.3 , 0.995, 0.995);
  TPad* vRat = new TPad("vRat", "", 0.005, 0.005 , 0.995, 0.3);
  vCan->SetTopMargin(0.07);
  vCan->SetBottomMargin(0.01);
	vCan->SetLeftMargin(0.15);
  vRat->SetTopMargin(0.02);
  vRat->SetBottomMargin(0.3);
	vRat->SetLeftMargin(0.15);
  vCan->Draw();
  vRat->Draw();

	TCanvas *c3 = new TCanvas("c3", "c3", 500, 500);
  TPad* yCan = new TPad("yCan", "", 0.005, 0.3 , 0.995, 0.995);
  TPad* yRat = new TPad("yRat", "", 0.005, 0.005 , 0.995, 0.3);
  yCan->SetTopMargin(0.07);
  yCan->SetBottomMargin(0.01);
	yCan->SetLeftMargin(0.15);
  yRat->SetTopMargin(0.02);
  yRat->SetBottomMargin(0.3);
	yRat->SetLeftMargin(0.15);
  yCan->Draw();
  yRat->Draw();

  TGraph* uGr = new TGraphErrors(numDriftBins,
                                 driftTimes[0], sigmaSqrVals[0],
                                 driftTimesErrs[0], sigmaSqrValsErrs[0]);

  TGraph* vGr = new TGraphErrors(numDriftBins,
                                 driftTimes[1], sigmaSqrVals[1],
                                 driftTimesErrs[1], sigmaSqrValsErrs[1]);

  TGraph* yGr = new TGraphErrors(numDriftBins,
                                 driftTimes[2], sigmaSqrVals[2],
                                 driftTimesErrs[2], sigmaSqrValsErrs[2]);

  // perform fits
  TF1* linFit = new TF1("linFit", "pol1");

  uGr->Fit("linFit", "", "", 0, maxTime-minTime);
  vGr->Fit("linFit", "", "", 0, maxTime-minTime);
  yGr->Fit("linFit", "", "", 0, maxTime-minTime);

  std::vector<TF1*> fitVec;
  fitVec.push_back(uGr->GetFunction("linFit"));
  fitVec.push_back(vGr->GetFunction("linFit"));
  fitVec.push_back(yGr->GetFunction("linFit"));

  // get diffusion values
  double uGrad    = fitVec.at(0)->GetParameter(1);
  double uGradErr = fitVec.at(0)->GetParError (1);
  double uChiSqr  = fitVec.at(0)->GetChisquare();
  double uNdf     = fitVec.at(0)->GetNDF();
  double uSig0    = fitVec.at(0)->Eval(0);
  double uSig0Err = fitVec.at(0)->GetParError(0);
  double vGrad    = fitVec.at(1)->GetParameter(1);
  double vGradErr = fitVec.at(1)->GetParError (1);
  double vChiSqr  = fitVec.at(1)->GetChisquare();
  double vNdf     = fitVec.at(1)->GetNDF();
  double vSig0    = fitVec.at(1)->Eval(0);
  double vSig0Err = fitVec.at(1)->GetParError(0);
  double yGrad    = fitVec.at(2)->GetParameter(1);
  double yGradErr = fitVec.at(2)->GetParError (1);
  double yChiSqr  = fitVec.at(2)->GetChisquare();
  double yNdf     = fitVec.at(2)->GetNDF();
  double ySig0    = fitVec.at(2)->Eval(0);
  double ySig0Err = fitVec.at(2)->GetParError(0);

  double uDiffV = uGrad    * std::pow(driftVelocity,2) * 1e6/2.;
  double uDiffE = uGradErr * std::pow(driftVelocity,2) * 1e6/2.;
  double vDiffV = vGrad    * std::pow(driftVelocity,2) * 1e6/2.;
  double vDiffE = vGradErr * std::pow(driftVelocity,2) * 1e6/2.;
  double yDiffV = yGrad    * std::pow(driftVelocity,2) * 1e6/2.;
  double yDiffE = yGradErr * std::pow(driftVelocity,2) * 1e6/2.;

  // loop all of the information we have and get values for ratios
  double sigmaSqrValsRatio     [3][numDriftBins];
  double sigmaSqrValsErrsRatio [3][numDriftBins];

  for (int ipl = 0; ipl < 3; ++ipl){
    for (int idb = 0; idb < numDriftBins; ++idb){

      double eval = fitVec.at(ipl)->Eval(driftTimes[ipl][idb]);
      sigmaSqrValsRatio     [ipl][idb] = (sigmaSqrVals[ipl][idb] - eval)/eval;
      sigmaSqrValsErrsRatio [ipl][idb] = (sigmaSqrValsErrs [ipl][idb]/eval);
    }
  }

  TGraph* uGrR = new TGraphErrors(numDriftBins,
                                  driftTimes[0], sigmaSqrValsRatio[0],
                                  driftTimesErrs[0], sigmaSqrValsErrsRatio[0]);

  TGraph* vGrR = new TGraphErrors(numDriftBins,
                                  driftTimes[1], sigmaSqrValsRatio[1],
                                  driftTimesErrs[1], sigmaSqrValsErrsRatio[1]);

  TGraph* yGrR = new TGraphErrors(numDriftBins,
                                  driftTimes[2], sigmaSqrValsRatio[2],
                                  driftTimesErrs[2], sigmaSqrValsErrsRatio[2]);


  // set histogram styles
  styleGraph(uGr , 0.0, 2300, 0.0001   , 5       );
  styleGraph(uGrR, 0.0, 2300, -0.09999 , 0.09999 );
  styleGraph(vGr , 0.0, 2300, 0.0001   , 5       );
  styleGraph(vGrR, 0.0, 2300, -0.09999 , 0.09999 );
  styleGraph(yGr , 0.0, 2300, 0.0001   , 5       );
  styleGraph(yGrR, 0.0, 2300, -0.09999 , 0.09999 );

  // set fit styles
  fitVec.at(0)->SetLineColor(kAzure+1);
  fitVec.at(0)->SetNpx(1000);
  fitVec.at(1)->SetLineColor(kGreen+1);
  fitVec.at(1)->SetNpx(1000);
  fitVec.at(2)->SetLineColor(kRed);
  fitVec.at(2)->SetNpx(1000);

  uCan->cd();
  uGr->SetTitle(";;#sigma^{2}");
  uGr->Draw("ap");
  fitVec.at(0)->Draw("same");
  if (nWvfmsVec.at(0)->Integral() > 0){
    nWvfmsVec.at(0)->Scale(2./nWvfmsVec.at(0)->GetMaximum());
    nWvfmsVec.at(0)->SetFillColor(kGray);
    nWvfmsVec.at(0)->SetLineWidth(0);
    nWvfmsVec.at(0)->GetXaxis()->SetLimits(0,maxTime-minTime);
    nWvfmsVec.at(0)->Draw("hist same");
  }
  drawPaveText("U Plane", uDiffV, uDiffE, uChiSqr, uNdf, uSig0, uSig0Err);

  uRat->cd();
  uRat->SetGridy();
  uGrR->SetTitle(";Drift time (#mus);(#sigma^{2}-Fit)/#sigma^{2}");
  uGrR->Draw("ap");

  TF1* lin = new TF1("lin", "0.", 0, maxTime-minTime); 
  lin->SetLineColor(kAzure+1);
  lin->DrawClone("same");

  vCan->cd();
  vGr->SetTitle(";;#sigma^{2}");
  vGr->Draw("ap");
  fitVec.at(1)->Draw("same");
  if (nWvfmsVec.at(1)->Integral() > 0){
    nWvfmsVec.at(1)->Scale(2./nWvfmsVec.at(1)->GetMaximum());
    nWvfmsVec.at(1)->SetFillColor(kGray);
    nWvfmsVec.at(1)->SetLineWidth(0);
    nWvfmsVec.at(1)->GetXaxis()->SetLimits(0,maxTime-minTime);
    nWvfmsVec.at(1)->Draw("hist same");
  }
  drawPaveText("V Plane", vDiffV, vDiffE, vChiSqr, vNdf, vSig0, vSig0Err);

  vRat->cd();
  vRat->SetGridy();
  vGrR->SetTitle(";Drift time (#mus);(#sigma^{2}-Fit)/#sigma^{2}");
  vGrR->Draw("ap");
  lin->SetLineColor(kGreen+1);
  lin->DrawClone("same");

  yCan->cd();
  yGr->SetTitle(";;#sigma^{2}");
  yGr->Draw("ap");
  fitVec.at(2)->Draw("same");
  if (nWvfmsVec.at(2)->Integral() > 0){
    nWvfmsVec.at(2)->Scale(2./nWvfmsVec.at(2)->GetMaximum());
    nWvfmsVec.at(2)->SetFillColor(kGray);
    nWvfmsVec.at(2)->SetLineWidth(0);
    nWvfmsVec.at(2)->GetXaxis()->SetLimits(0,maxTime-minTime);
    nWvfmsVec.at(2)->Draw("hist same");
    yCan->RedrawAxis();
  }
  drawPaveText("Y Plane", yDiffV, yDiffE, yChiSqr, yNdf, ySig0, ySig0Err);

  yRat->cd();
  yRat->SetGridy();
  yGrR->SetTitle(";Drift time (#mus);(#sigma^{2}-Fit)/#sigma^{2}");
  yGrR->Draw("ap");
  lin->SetLineColor(kRed);
  lin->DrawClone("same");

	c1->cd();
  c1->SaveAs("output_uplane.png");
	c2->cd();
  c2->SaveAs("output_vplane.png");
	c3->cd();
  c3->SaveAs("output_yplane.png");

  fOutput->Close();
  fInput->Close();

}

int main(int argv, char** argc){

  std::string inputFileName(argc[1]);
  makePlot(inputFileName);

  return 0;

}
