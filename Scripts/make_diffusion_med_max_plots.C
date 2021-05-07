#include "StylePlots.h"

float GetMedian(TH1D* h){

  Double_t x,q;
  q = 0.5;
  h->GetQuantiles(1,&x,&q);
  return x*x;

}

std::vector<float> getMedians(TH2D* h){
  std::vector<float> medians;
  for (int i = 1; i < h->GetNbinsX()+1; ++i){
    std::string name = "h_" + std::string(h->GetName()) + std::to_string(i);
    TH1D* h2 = h->ProjectionY(name.c_str(), i, i);
    medians.push_back(GetMedian(h2));
  }
  return medians;
}

void make_diffusion_med_max_plots(){

  TFile* f2 = new TFile("/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionAna/diffmod_run3_crt_Aug2020_newFV_bugFix.root", "read");
  TH2D* precut  = (TH2D*)f2->Get("DiffusionModule/plane2/h_sigma_v_bin_precutplane2");
  TH2D* postcut = (TH2D*)f2->Get("DiffusionModule/plane2/h_sigma_v_bin_postcutplane2");

  std::vector precutMedians  = getMedians(precut);
  std::vector postcutMedians = getMedians(postcut);

  TH1D* diffs = new TH1D("diffs", "postcut-precut median sigma(mus)", 25, 0, 25);
  for (int i = 0; i < 25; ++i){
    std::cout << "median pre " << precutMedians[i] << " post: " << postcutMedians[i] << std::endl;
    diffs->SetBinContent(i+1, postcutMedians[i]-precutMedians[i]);
  }

  //diffs->GetYaxis()->SetRangeUser(-0.01, 0.025);
  diffs->Draw();

  //SetGenericStyle();
  //gROOT->ForceStyle();



}
