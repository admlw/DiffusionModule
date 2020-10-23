//.............................................................................
void centerAndArbitrizeHistogram(TH1D* h){

  // this centers the summed waveform in the histogram bounds
  // may want to replace this with a gaussian fit?
  int maxBin    = h->GetMaximumBin();
  int halfNBins = h->GetNbinsX()/2.;

  if (maxBin >= halfNBins){
    for (int j = 0; j <= h->GetNbinsX(); ++j){
      double newVal = h->GetBinContent(j);
      h->SetBinContent(j + halfNBins - maxBin, newVal);
    }
  }

  if (maxBin < halfNBins){
    for (int j = h->GetNbinsX(); j >= 0; j--){
      double newVal = h->GetBinContent(j);
      h->SetBinContent(j + halfNBins - maxBin, newVal);
    }
  }

  h->GetXaxis()->SetLimits(-92, 92);
  h->GetXaxis()->SetRangeUser(-20, 20);

}

//.............................................................................
std::vector<TH1D*> getHistograms(TFile* file, int planeNo, int color){

  std::vector<TH1D*> returnHists;

  for (int i = 0; i < 25; ++i){
    std::string histName = "DiffusionModule/plane" + std::to_string(planeNo) +
                           "/summed_waveform_bin_" + std::to_string(i) +
                           "_plane" + std::to_string(planeNo);

    if ((TH1D*)file->Get(histName.c_str()))
      returnHists.push_back((TH1D*)file->Get(histName.c_str()));
    else{
      std::cout << "-- No summed histogram in bin" << std::endl;
      continue;
    }

    returnHists[i]->SetLineColor(color);
    returnHists[i]->SetLineWidth(2);
    returnHists[i]->SetTitle(";Arbitary Time (ticks); ADCs (norm.)");

    returnHists[i]->GetXaxis()->CenterTitle();
    returnHists[i]->GetXaxis()->SetTitleFont(43);
    returnHists[i]->GetXaxis()->SetLabelFont(43);
    returnHists[i]->GetXaxis()->SetTitleSize(30);
    returnHists[i]->GetXaxis()->SetLabelSize(30);
    returnHists[i]->GetXaxis()->SetTitleOffset(2.6);

    returnHists[i]->GetYaxis()->CenterTitle();
    returnHists[i]->GetYaxis()->SetTitleFont(43);
    returnHists[i]->GetYaxis()->SetLabelFont(43);
    returnHists[i]->GetYaxis()->SetTitleSize(30);
    returnHists[i]->GetYaxis()->SetLabelSize(30);
    returnHists[i]->GetYaxis()->SetTitleOffset(2.6);

    centerAndArbitrizeHistogram(returnHists[i]);
  }

  return returnHists;
}

//.............................................................................
std::vector<TH1D*> getOffsetHistograms(std::vector<TH1D*>& hists){

  std::vector<TH1D*> returnerHists;

  float baseline = 2300/50;
  for (int i = 0; i < hists.size(); ++i){
    returnerHists.push_back((TH1D*)hists[i]->Clone());
    returnerHists[i]->Scale(1000./(returnerHists[i]->Integral()));
    for (int j = 0; j <= returnerHists[i]->GetNbinsX(); ++j){

      returnerHists[i]->SetBinContent(j, 
                                      returnerHists[i]->GetBinContent(j) + baseline); 
    }
    baseline += (2300./25);
  }

  return returnerHists;
}

//.............................................................................
std::vector<TF1*> getLines(std::vector<TH1D*>& hists){

  std::vector<TF1*> returnFuncs;

  float baseline = 2300/50;
  for (int i = 0; i < 25; ++i){
    std::string funcName = "line" + std::to_string(i);
    returnFuncs.push_back(new TF1(funcName.c_str(), "[0]", -92, 92));
    returnFuncs[i]->SetParameter(0, baseline);
    returnFuncs[i]->SetLineColor(kGray+2);
    returnFuncs[i]->SetLineWidth(1);
    returnFuncs[i]->SetLineStyle(2);

    baseline += (2300./25);
  }

  return returnFuncs;
}

//.............................................................................
void make_classic_plot(){

  int planeNo = 2;

  // style opts
  gStyle->SetOptStat(0);

  // input files
  TFile* fData   = new TFile("diffmod_run3_crt_Aug2020_newFV_bugFix.root", "READ");
  TFile* fMC6p4  = new TFile("diffmod_single_mu_6p4.root" , "READ");
  TFile* fMC3p23 = new TFile("diffmod_single_mu_3p23.root", "READ");

  // pull out histograms
  std::vector<TH1D*> hData   = getHistograms(fData   , planeNo, kBlack);
  std::vector<TH1D*> hMC6p4  = getHistograms(fMC6p4  , planeNo, kAzure+1);
  std::vector<TH1D*> hMC3p23 = getHistograms(fMC3p23 , planeNo, kGreen+1);

  std::vector<TH1D*> hOffsetData = getOffsetHistograms(hData);
  std::vector<TF1*>  lines       = getLines(hOffsetData);

  // setup canvas
  TCanvas *c1 = new TCanvas("c1","",600,1000);                  
  c1->SetLeftMargin(0.2);
  c1->SetTopMargin(0.05);
  c1->SetRightMargin(0.05);

  TCanvas* c2   = new TCanvas("c2", "", 600, 1000);
  TPad *rtPad   = new TPad("rtPad"  , "", 0.00, 0.68, 1.00, 1.00);
  TPad *rmPad   = new TPad("rmPad"  , "", 0.00, 0.37, 1.00, 0.68);
  TPad *rbPad   = new TPad("rbPad"  , "", 0.00, 0.03 , 1.00, 0.37);
  rtPad->SetBottomMargin(0);
  rmPad->SetTopMargin(0);
  rmPad->SetBottomMargin(0);
  rbPad->SetTopMargin(0);
  rbPad->SetBottomMargin(0.17);
  rtPad->SetLeftMargin(0.2);
  rmPad->SetLeftMargin(0.2);
  rbPad->SetLeftMargin(0.2);

  rtPad->Draw(); rmPad->Draw(); rbPad->Draw();

  c1->cd();
  TH2D* bg = new TH2D("bg", ";Arbitary Time (ticks);Drift Time (#mus)",
                      40, -20, 20, 100, 0, 2800);

  bg->GetXaxis()->CenterTitle();
  bg->GetYaxis()->CenterTitle();
  bg->GetXaxis()->SetTitleSize(0.06);
  bg->GetYaxis()->SetTitleSize(0.06);
  bg->GetXaxis()->SetLabelSize(0.06);
  bg->GetYaxis()->SetLabelSize(0.06);
  bg->GetXaxis()->SetTitleOffset(0.8);
  bg->GetYaxis()->SetTitleOffset(1.6);

  bg->Draw();
    
  for (int i = 0; i < hData.size(); ++i){
    hOffsetData[i]->Draw("hist same");
    lines[i]->Draw("same");
  }

  TPaveText* pt1 = new TPaveText(0.225, 0.85, 0.70, 0.95, "NDC");
  pt1->AddText("MicroBooNE Run 3 Cosmic Data");
  pt1->AddText("E = 273 V/cm");
  pt1->AddText("Y Plane");
  pt1->SetFillStyle(0);
  pt1->SetBorderSize(0);
  pt1->SetTextSize(0.04);
  pt1->SetTextAlign(12);
  pt1->Draw("same");    

  rtPad->cd();

  // the 3.23 one for bin24 need shifted right so deal with that here
  //for (int i = 0; i < hMC3p23[24]->GetNbinsX(); ++i){
  //    hMC3p23[24]->SetBinContent(i,
  //                               hMC3p23[24]->GetBinContent(i));
  //}


  hMC3p23[24]->DrawNormalized("");
  hMC6p4 [24]->DrawNormalized("same hist");
  hData  [24]->DrawNormalized("same hist");
  
  TPaveText* ptt = new TPaveText(0.7, 0.80, 0.89, 0.85, "NDC");
  ptt->AddText("Drift Time = 2254 #mus");
  ptt->SetFillStyle(0);
  ptt->SetBorderSize(0);
  ptt->SetTextSize(0.06);
  ptt->SetTextAlign(32);
  ptt->Draw("same");    

  rmPad->cd();
  hMC3p23[12]->DrawNormalized("");
  hMC6p4 [12]->DrawNormalized("same hist");
  hData  [12]->DrawNormalized("same hist");
  
  TPaveText* ptm = new TPaveText(0.7, 0.89, 0.89, 0.94, "NDC");
  ptm->AddText("Drift Time = 1150 #mus");
  ptm->SetFillStyle(0);
  ptm->SetBorderSize(0);
  ptm->SetTextSize(0.06);
  ptm->SetTextAlign(32);
  ptm->Draw("same");    

  rbPad->cd();
  hMC3p23[0]->DrawNormalized("");
  hMC6p4 [0]->DrawNormalized("same hist");
  hData  [0]->DrawNormalized("same hist");

  TPaveText* ptb = new TPaveText(0.7, 0.89, 0.89, 0.94, "NDC");
  ptb->AddText("Drift Time = 46 #mus");
  ptb->SetFillStyle(0);
  ptb->SetBorderSize(0);
  ptb->SetTextSize(0.06);
  ptb->SetTextAlign(32);
  ptb->Draw("same");    

  TLegend* leg = new TLegend(0.58, 0.67, 0.87, 0.87);
  leg->AddEntry(hMC3p23[0], "MC, D_{L} = 3.23 cm^{2}/s", "l");
  leg->AddEntry(hMC6p4[0] , "MC, D_{L} = 6.40 cm^{2}/s", "l");
  leg->AddEntry(hData[0], "MicroBooNE Data", "l");
  leg->SetLineWidth(0);
  leg->SetFillStyle(0);
  leg->Draw("same");

  c1->SaveAs("ClassicPlot.pdf");
  c1->SaveAs("ClassicPlot.png");
  c2->SaveAs("ThreeWaveformComparison.pdf");
  c2->SaveAs("ThreeWaveformComparison.png");

}
