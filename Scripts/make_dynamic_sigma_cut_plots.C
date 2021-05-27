#include "StylePlots.h"

// partition canvas

void areaNormaliseX(TH2* h){

  for (int i = 0; i < h->GetNbinsX()+1; ++i){

    float content = 0.0;
    for (int j = 0; j < h->GetNbinsY()+1; ++j){
      content += h->GetBinContent(i,j);
    }

    for (int j = 0; j < h->GetNbinsY()+1; ++j){
      h->SetBinContent(i,j,h->GetBinContent(i,j)/content);
    }

  }

}

void SaveCanvas(TH2D* h, TH1D* h2, Int_t textColor){
  std::cout << "SaveCanvas" << std::endl;

  std::string name = (std::string)h->GetName();
  TCanvas* c1 = new TCanvas(name.c_str(), name.c_str(), 500, 600);
  TPad *topPad = new TPad("topPad", "", 0.005, 0.36, 0.995, 0.995);
  TPad *bottomPad = new TPad("bottomPad", "", 0.005, 0.005, 0.995, 0.36);
  topPad->SetTopMargin(0.1);
  topPad->SetBottomMargin(0.005);
  topPad->SetRightMargin(0.15);
  topPad->SetLeftMargin(0.15);
  bottomPad->SetTopMargin(0.005);
  bottomPad->SetBottomMargin(0.25);
  bottomPad->SetRightMargin(0.15);
  bottomPad->SetLeftMargin(0.15);
  topPad->Draw();
  bottomPad->Draw();

  topPad->cd();
  areaNormaliseX(h);
  h->SetContour(1000);
  h->GetYaxis()->SetLimits(0,5);
  h->GetYaxis()->SetRangeUser(0.55,2.99);
  h->GetZaxis()->SetRangeUser(0.000005,0.22);
  h->GetXaxis()->SetLimits(0,2300);
  h->GetXaxis()->SetTitle("Drift Time (#mus)");
  h->GetYaxis()->SetTitle("Hit Width (#mus)");
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(22);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(22);

  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitleOffset(10000);
  h->GetXaxis()->SetLabelOffset(10000);
  //h->UseCurrentStyle();
  h->Draw("colz");
  ApplyLabel(DataType::kData, 0.83, -1, textColor);

  bottomPad->cd();
  h2->GetXaxis()->SetLimits(0,2300);
  h2->GetXaxis()->SetTitle("Drift Time (#mus)");
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->CenterTitle();
  h2->GetYaxis()->SetTitleFont(43);
  h2->GetYaxis()->SetTitleSize(22);
  h2->GetYaxis()->SetTitleOffset(1.6);
  h2->GetXaxis()->SetTitleOffset(2.5);
  h2->GetXaxis()->SetNdivisions(505);
  h2->GetYaxis()->SetNdivisions(505);
  h2->GetXaxis()->SetTitleFont(43);
  h2->GetXaxis()->SetTitleSize(22);
  h2->GetYaxis()->SetLabelFont(43);
  h2->GetYaxis()->SetLabelSize(22);
  h2->GetXaxis()->SetLabelFont(43);
  h2->GetXaxis()->SetLabelSize(22);
    h2->Scale(1./(std::pow(10,4)));
    h2->GetYaxis()->SetRangeUser(0, 85.);
    h2->GetYaxis()->SetTitle("x10^{4} Hits");
  h2->SetLineWidth(2);
  h2->SetLineColor(kBlack);
  h2->SetFillColor(kBlack);
  h2->SetFillStyle(3345);
  //h2->UseCurrentStyle();
  h2->Draw("hist");

  c1->SaveAs((name+".pdf").c_str());

}

void make_dynamic_sigma_cut_plots(){

  std::cout << "start" << std::endl;

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBlueGreenYellow);

  std::cout << "set styles" << std::endl;

  TFile* f2 = new TFile("/pnfs/uboone/persistent/diffusion_analysis_final_files/diffmod_run3data_paper.root", "read");

  TTree* t2 = (TTree*)f2->Get("DiffusionModule/difftree");
  std::cout << "getting hists" << std::endl;
  //TH1D* nWvfmsPreCut = new TH1D("h_nWvfmsInBinplane2PostCut", ";;", 25, 0, 25);
  //t2->Draw("wvfm_bin_no >> h_nWvfmsInBinplane2PostCut", "track_avg_trans_dist < 6");

  TH2D* precut  = (TH2D*)f2->Get("DiffusionModule/plane2/h_sigma_v_bin_precutplane2");
  TH2D* postcut = (TH2D*)f2->Get("DiffusionModule/plane2/h_sigma_v_bin_postcutplane2");
  TH1D* nWvfmsPreCut = (TH1D*)precut->ProjectionX();
  TH1D* nWvfmsPostCut = (TH1D*)f2->Get("DiffusionModule/plane2/h_nWvfmsInBinplane2");

  std::cout << "got hists" << std::endl;

  SetGenericStyle();
  SaveCanvas(precut , nWvfmsPreCut, kWhite);
  SaveCanvas(postcut, nWvfmsPostCut, kBlack);

}

