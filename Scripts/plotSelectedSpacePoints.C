#include "StylePlots.h"

void plotSelectedSpacePoints(){

  SetGenericStyle();

  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBlueGreenYellow);

  TFile* f = new TFile("/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionAna/diffmod_run3_crt_Aug2020_newFV_bugFix.root", "read");

  TTree* t = (TTree*)f->Get("DiffusionModule/difftree");

  TH2D* sps = new TH2D("sps", ";Z (cm);Y (cm);Number of Space Points", 500, 0, 1036, 500, -116.5, 116.5);
  sps->SetContour(1000);
  sps->GetXaxis()->CenterTitle();
  sps->GetYaxis()->CenterTitle();
  sps->GetZaxis()->CenterTitle();
  t->Draw("sp_y:sp_z >> sps", "hit_view == 2");


  TCanvas *c1 = new TCanvas("c1", "c1", 1100, 600);
  c1->SetRightMargin(0.15);
  c1->cd();
  sps->GetZaxis()->SetRangeUser(0,360);
  sps->Draw("colz");
 
 // TLine* ltop = new TLine(0, 116.5, 1036, 116.5);
 // ltop->SetLineWidth(2);
 // ltop->Draw("same");

 // TLine* lbot = new TLine(0, -116.5, 1036, -116.5);
 // lbot->SetLineWidth(2);
 // lbot->Draw("same");

 // TLine* ll = new TLine(0, -116.5, 0, +116.5);
 // ll->SetLineWidth(2);
 // ll->Draw("same");

 // TLine* lr = new TLine(1036, -116.5, 1036, +116.5);
 // lr->SetLineWidth(2);
 // lr->Draw("same");
 ApplyLabel(DataType::kData, 0.39);



  c1->SaveAs("SpacePointPositionsYZ.png");
  c1->SaveAs("SpacePointPositionsYZ.pdf");

}
