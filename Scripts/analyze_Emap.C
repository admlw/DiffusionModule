#include "StylePlots.h"

void analyze_Emap() {

  SetGenericStyle();
  TString dir = "/uboone/app/users/amogan/diffusion_mcc9/workdir/";
  TFile *fin = new TFile(dir+"Emap-NTT-500-N3-S500_laserdata_v1098_cathode2548_vectorstart_notshiftcurve.root", "read");
  //TFile *fin = new TFile(dir+"Emap-NTT-1-MergedMapsSmoothCosmicAndLaserNoDriftVVolumeSmoothed_newCurve.root", "read");

  TH3 *hin = (TH3*)fin->Get("Distorted_v_X");
  // Slice in X near anode (26 total bins on z axis)
  hin->GetXaxis()->SetRange(1,1);
  TH2D *hproj = (TH2D*)hin->Project3D("yz");
  TH1D *h_vdist = new TH1D("h_vdist", "", 50, 1.0, 1.2);

  // Initialize some stuff
  //double vd = 1.098; // mm/us
  double vd = 1.0762; // mm/us
  double vi = 0.;
  double tmp_v = 0.;
  bool useFV = true;
  TString vdstring = Form("v_{d} = %0.3f mm/#mus", vd);

  // Loop over bins, calculate (v_i - v_0)/v_0
  for (int i = 1; i < hproj->GetNbinsX()+1; i++) {
    for (int j = 1; j < hproj->GetNbinsY()+1; j++) {
      vi = hproj->GetBinContent(i, j); 
      //std::cout << "vi = " << vi << std::endl;
      //std::cout << "i, j  = " << i << ", " << j << std::endl;
      //std::cout << "vd - vi = " << vd - vi << std::endl;
      tmp_v = (vi - vd)/vd * 100.; 
      if (useFV) {
        if (hproj->GetXaxis()->GetBinLowEdge(i) < 400. ||
           (hproj->GetXaxis()->GetBinLowEdge(i) > 675. && hproj->GetXaxis()->GetBinLowEdge(i) < 775.) ||
            hproj->GetXaxis()->GetBinLowEdge(i) > 950.)
            //tmp_v = 0.;
            tmp_v = -999;
      }
      hproj->SetBinContent(i, j, tmp_v);
      if (vi >= 1.0979 && vi <= 1.0981) {
        //std::cout << "Skipping bin " << i << ", " << j << std::endl;
        continue;
      }
      h_vdist->Fill(vi);
    }
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 600);
  gStyle->SetPalette(kRainBow);
  //gStyle->SetTitleAlign(13);
  
  hproj->GetXaxis()->SetTitle("Z (cm)");
  hproj->GetYaxis()->SetTitle("Y (cm)");
  hproj->GetZaxis()->SetRangeUser(-3, 3);
  hproj->SetTitle("(v_{x} - v_{d})/v_{d} [%] at X = 10 cm, "+vdstring);

  hproj->UseCurrentStyle();
  hproj->Draw("colz");
  c1->SaveAs("vmap.pdf", "PDF");

  /*
  TCanvas *c2 = new TCanvas;
  gStyle->SetOptStat(1101);
  h_vdist->SetLineColor(kAzure+2);
  h_vdist->GetXaxis()->SetTitle("v_{d} [mm/#mus]");
  h_vdist->Draw();

  Int_t maxval = h_vdist->GetMaximum();
  TLine *l = new TLine(1.098, 0., 1.098, maxval);
  l->SetLineColor(kRed+1);
  l->SetLineStyle(5);
  l->SetLineWidth(2);
  l->Draw("same");

  TPaveText* tpv = new TPaveText(0.16, 0.69, 0.7, 0.89, "NDC");
  tpv->SetTextAlign(11);
  tpv->SetFillStyle(0);
  tpv->SetLineWidth(0);
  tpv->SetTextFont(43);
  tpv->SetTextSize(24);
  tpv->AddText("X = 10 cm");
  tpv->Draw("same");

  std::cout << "1D peak val: " << h_vdist->GetXaxis()->GetBinUpEdge(h_vdist->GetMaximumBin() ) << std::endl;
  std::cout << "1D mean: " << h_vdist->GetMean() << std::endl;

  c2->SaveAs("vdist.pdf", "PDF");
  */




}
