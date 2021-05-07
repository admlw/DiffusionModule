#include "StylePlots.h"

void NormalizeAndDraw(TH1D* h, std::string opt){
  h->Scale(1./h->Integral());
  h->GetYaxis()->SetRangeUser(-0.008, 0.175);
  h->GetXaxis()->SetRangeUser(-19, 19);
  h->GetYaxis()->SetTitleOffset(1000);
  h->GetXaxis()->SetTitleOffset(10000);
  h->GetXaxis()->SetNdivisions(305);
  h->GetYaxis()->SetNdivisions(505);
  h->Draw(opt.c_str());
}

// partition canvas
void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
   if (!C) return;

   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;

   Float_t hSpacing = 0.0;
   Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;

   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;

   for (Int_t i=0;i<Nx;i++) {

      if (i==0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr-hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx-1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr-hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }

      for (Int_t j=0;j<Ny;j++) {

         if (j==0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu-vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         } else if (j == Ny-1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = tMargin / (vposu-vposd);
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }

         C->cd(0);

         char name[16];
         sprintf(name,"pad_%i_%i",i,j);
         TPad *pad = (TPad*) gROOT->FindObject(name);
         if (pad) delete pad;
         pad = new TPad(name,"",hposl,vposd,hposr,vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);

         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);

         pad->Draw();
      }
   }
}

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
std::vector<TH1D*> getHistograms(TFile* file, int planeNo, int color, int style=1, int width=2){

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
    returnHists[i]->SetLineStyle(style);
    returnHists[i]->SetLineWidth(width);
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

  std::vector<std::vector<TH1D*>> hDataTP;
  std::vector<std::vector<TH1D*>> hMC6p4TP;
  std::vector<std::vector<TH1D*>> hMC3p23TP;

  SetGenericStyle();

  for (int planeNo = 0; planeNo < 3; ++planeNo){

    // style opts
    gStyle->SetOptStat(0);

    // input files
    std::string dir = "/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionAna/";
    TFile* fData   = new TFile(std::string(dir+"diffmod_run3_crt_Aug2020_newFV_bugFix.root").c_str(), "READ");
    TFile* fMC6p4  = new TFile(std::string(dir+"diffmod_single_muons_angular_June2020.root").c_str(), "READ");
    //TFile* fMC3p23 = new TFile(std::string("/uboone/app/users/amogan/diffusion_mcc9/workdir/diffmod_single_muons_DL_down_50percent.root").c_str(), "READ");
    TFile* fMC3p23 = new TFile(std::string("/pnfs/uboone/persistent/users/alister1/diffusion_paper/diffmod_sim_DL_3.74.root").c_str(), "READ");


    // pull out histograms
    std::vector<TH1D*> hData   = getHistograms(fData   , planeNo, kBlack);
    std::vector<TH1D*> hMC6p4  = getHistograms(fMC6p4  , planeNo, kPTVibrantCyan, kDashed, 3);
    std::vector<TH1D*> hMC3p23 = getHistograms(fMC3p23 , planeNo, kPTVibrantMagenta, kDashed, 3);

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
    TH2D* bg = new TH2D("bg", ";Time (ticks, Zero-centered);Drift Time (#mus)",
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

    TPaveText* pt1 = new TPaveText(0.225, 0.83, 0.70, 0.93, "NDC");
    pt1->AddText("MicroBooNE Cosmic Data");
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

    hDataTP  .push_back(hData);
    hMC3p23TP.push_back(hMC3p23);
    hMC6p4TP .push_back(hMC6p4);

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

    c1->SaveAs(std::string("ClassicPlot"+std::to_string(planeNo)+".pdf").c_str());
    c1->SaveAs(std::string("ClassicPlot"+std::to_string(planeNo)+".png").c_str());
    c2->SaveAs(std::string("ThreeWaveformComparison"+std::to_string(planeNo)+".pdf").c_str());
    c2->SaveAs(std::string("ThreeWaveformComparison"+std::to_string(planeNo)+".png").c_str());

  }

  TCanvas* tp = new TCanvas("tp", "tp", 1000, 1000);
  CanvasPartition(tp, 3, 3, 0.1, 0.1, 0.15, 0.1);
  //pad_0_0
  TPad* pad_0_0 = ((TPad*)gROOT->FindObject("pad_0_0"));
  pad_0_0->SetRightMargin(0.04);
  pad_0_0->cd();

  std::string drawopt = "hist same l";
  NormalizeAndDraw(hDataTP[0][0], drawopt);
  NormalizeAndDraw(hMC6p4TP[0][0],drawopt);
  NormalizeAndDraw(hMC3p23TP[0][0], drawopt);
  //pad_0_1
  TPad* pad_0_1 = ((TPad*)gROOT->FindObject("pad_0_1"));
  pad_0_1->SetRightMargin(0.04);
  pad_0_1->cd();
  NormalizeAndDraw(hDataTP[0][12],drawopt);
  NormalizeAndDraw(hMC6p4TP[0][12],drawopt);
  NormalizeAndDraw(hMC3p23TP[0][12],drawopt);

  //pad_0_2
  TPad* pad_0_2 = ((TPad*)gROOT->FindObject("pad_0_2"));
  pad_0_2->SetRightMargin(0.04);
  pad_0_2->cd();
  NormalizeAndDraw(hDataTP[0][24],drawopt);
  NormalizeAndDraw(hMC6p4TP[0][24],drawopt);
  NormalizeAndDraw(hMC3p23TP[0][24],drawopt);

  //pad_1_0
  TPad* pad_1_0 = ((TPad*)gROOT->FindObject("pad_1_0"));
  pad_1_0->SetRightMargin(0.04);
  pad_1_0->cd();
  NormalizeAndDraw(hDataTP[1][0], drawopt);
  NormalizeAndDraw(hMC6p4TP[1][0],drawopt);
  NormalizeAndDraw(hMC3p23TP[1][0], drawopt);

  //pad_1_1
  TPad* pad_1_1 = ((TPad*)gROOT->FindObject("pad_1_1"));
  pad_1_1->SetRightMargin(0.04);
  pad_1_1->cd();
  NormalizeAndDraw(hDataTP[1][12], drawopt);
  NormalizeAndDraw(hMC6p4TP[1][12],drawopt);
  NormalizeAndDraw(hMC3p23TP[1][12], drawopt);

  //pad_1_2
  TPad* pad_1_2 = ((TPad*)gROOT->FindObject("pad_1_2"));
  pad_1_2->SetRightMargin(0.04);
  pad_1_2->cd();
  NormalizeAndDraw(hDataTP[1][24], drawopt);
  NormalizeAndDraw(hMC6p4TP[1][24],drawopt);
  NormalizeAndDraw(hMC3p23TP[1][24], drawopt);

  //pad_2_0
  TPad* pad_2_0 = ((TPad*)gROOT->FindObject("pad_2_0"));
  //pad_2_0->SetRightMargin(0.04);
  pad_2_0->cd();
  NormalizeAndDraw(hDataTP[2][0], drawopt);
  NormalizeAndDraw(hMC6p4TP[2][0],drawopt);
  NormalizeAndDraw(hMC3p23TP[2][0], drawopt);

  //pad_2_1
  TPad* pad_2_1 = ((TPad*)gROOT->FindObject("pad_2_1"));
  //pad_2_1->SetRightMargin(0.04);
  pad_2_1->cd();
  NormalizeAndDraw(hDataTP[2][12], drawopt);
  NormalizeAndDraw(hMC3p23TP[2][12], drawopt);
  NormalizeAndDraw(hMC6p4TP[2][12],drawopt);

  //pad_2_2
  TPad* pad_2_2 = ((TPad*)gROOT->FindObject("pad_2_2"));
  //pad_2_2->SetRightMargin(0.04);
  pad_2_2->cd();

  NormalizeAndDraw(hDataTP[2][24], drawopt);
  NormalizeAndDraw(hMC6p4TP[2][24],drawopt);
  NormalizeAndDraw(hMC3p23TP[2][24], drawopt);

  tp->cd();
  TLatex* xaxis = new TLatex(0.49, 0.1, "Time (ticks, Zero-centered)");
  xaxis->SetNDC();
  xaxis->SetTextSize(0.03);
  xaxis->SetTextAlign(22);
  xaxis->Draw();

  TLatex* yaxis = new TLatex(0.02, 0.55, "Waveform Amplitude (area normalized)");
  yaxis->SetNDC();
  yaxis->SetTextAngle(90);
  yaxis->SetTextSize(0.03);
  yaxis->SetTextAlign(22);
  yaxis->Draw();

  TLatex *upl = new TLatex(0.18, 0.91, "U Plane");
  upl->SetNDC();
  upl->SetTextSize(0.03);
  upl->Draw();
  TLatex *vpl = new TLatex(0.45, 0.91, "V Plane");
  vpl->SetNDC();
  vpl->SetTextSize(0.03);
  vpl->Draw();
  TLatex *ypl = new TLatex(0.72, 0.91, "Y Plane");
  ypl->SetNDC();
  ypl->SetTextSize(0.03);
  ypl->Draw();

  TLatex *fd = new TLatex(0.91, 0.9, "Drift Time = 2254 #mus");
  fd->SetNDC();
  fd->SetTextSize(0.03);
  fd->SetTextAngle(270);
  fd->Draw();
  TLatex *md = new TLatex(0.91, 0.58, "1150 #mus");
  md->SetNDC();
  md->SetTextSize(0.03);
  md->SetTextAngle(270);
  md->Draw();
  TLatex *nd = new TLatex(0.91, 0.32, "46 #mus");
  nd->SetNDC();
  nd->SetTextSize(0.03);
  nd->SetTextAngle(270);
  nd->Draw();

  TLegend* leg2 = new TLegend(0.11, 0.81, 0.31, 0.9);
  leg2->AddEntry(hDataTP[2][24], "MicroBooNE Data", "l");
  leg2->AddEntry(hMC3p23TP[2][24], "Sim. D_{L} = 3.74 cm^{2}/s", "l");
  leg2->AddEntry(hMC6p4TP[2][24], "Sim. D_{L} = 6.40 cm^{2}/s", "l");
  leg2->SetLineWidth(0);
  leg2->SetFillStyle(0);
  leg2->Draw("same");



  tp->SaveAs("tp_comparison.pdf");

}
