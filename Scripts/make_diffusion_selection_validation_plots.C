#include "StylePlots.h"

void make_diffusion_selection_validation_plots(){
  SetGenericStyle();
  TFile* f  = new TFile("/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionFiltered/diffusionSelectionInfo_run3_crt_Aug2020_newFV.root", "read");
  TFile* f2 = new TFile("/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionAna/diffmod_run3_crt_Aug2020_newFV_bugFix.root", "read");

  std::cout << "Starting script" << std::endl;
  TTree* t  = (TTree*)f->Get("diffsel/difffiltertree");
  TTree* t2 = (TTree*)f2->Get("DiffusionModule/difftree");
  //gStyle->SetOptStat(0);

  // define histograms
  TH2D* allTracks_thxz_thyz = new TH2D("allTracks_thxz_thyz", 
                                       ";#theta_{xz} (deg);#theta_{yz} (deg)", 
                                       100, -180, 180, 
                                       100, -180, 180);

  TH2D* lengthTracks_thxz_thyz = new TH2D("lengthTracks_thxz_thyz", 
                                       ";#theta_{xz} (deg);#theta_{yz} (deg)", 
                                       100, -180, 180, 
                                       100, -180, 180);
  
  TH2D* t0Tracks_thxz_thyz  = new TH2D("t0Tracks_thxz_thyz", 
                                       ";#theta_{xz} (deg);#theta_{yz} (deg)", 
                                       100, -180, 180, 
                                       100, -180, 180);

  TH2D* fidVolTracks_thxz_thyz  = new TH2D("fidVolTracks_thxz_thyz", 
                                       ";#theta_{xz} (deg);#theta_{yz} (deg)", 
                                       100, -180, 180, 
                                       100, -180, 180);

  TH2D* selTracks_thxz_thyz = new TH2D("selTracks_thxz_thyz", 
                                       ";#theta_{xz} (deg);#theta_{yz} (deg)", 
                                       100, -180, 180, 
                                       100, -180, 180);

  TH2D* selTracks_xy = new TH2D("selTracks_xy", 
                                ";Track Start x (cm);Track Start y (cm)", 
                                100, -10, 260,
                                100, -120, 120);

  TH2D* selTracks_xz = new TH2D("selTracks_xz", 
                                ";Track Start z (cm);Track Start x (cm)", 
                                100, -10, 1040,
                                100, -10, 260);

  TH2D* selTracks_yz = new TH2D("selTracks_yz", 
                                ";Track Start z (cm);Track Start y (cm)", 
                                100, -10, 1040,
                                100, -120, 120); 

  TH2D* selTracks_SCExy = new TH2D("selTracks_SCExy", 
                                   ";SCE-Corrected Track Start x (cm);SCE-Corrected Track Start y (cm)", 
                                   100, -10, 260,
                                   100, -120, 120);

  TH2D* selTracks_SCExz = new TH2D("selTracks_SCExz", 
                                   ";SCE-Corrected Track Start z (cm);SCE-Corrected Track Start x (cm)", 
                                   100, -10, 1040, 
                                   100, -10, 260);

  TH2D* selTracks_SCEyz = new TH2D("selTracks_SCEyz", 
                                   ";SCE-Corrected Track Start z (cm);SCE-Corrected Track Start y (cm)", 
                                   100, -10, 1040, 
                                   100, -120, 120);

  TH1D* allTracks_length = new TH1D("allTracks_length",
                                    ";Track length (cm); Num. tracks",
                                    20, 0, 1000);

  TH1D* lengthTracks_length = new TH1D("lengthTracks_length",
                                    ";Track length (cm); Num. tracks",
                                    20, 0, 1000);

  TH1D* t0Tracks_length = new TH1D("t0Tracks_length",
                                    ";Track length (cm); Num. tracks",
                                    20, 0, 1000);

  TH1D* fidVolTracks_length = new TH1D("fidVolTracks_length",
                                    ";Track length (cm); Num. tracks",
                                    20, 0, 1000);

  TH1D* selTracks_length = new TH1D("selTracks_length",
                                    ";Track length (cm); Num.tracks",
                                    20, 0, 1000);

  TH1D* defTracks_length = new TH1D("defTracks_length",
                                    ";Track length (cm); Num. tracks",
                                    20, 0, 1000);

  // save plots
  TCanvas* c1 = new TCanvas();
  c1->SetRightMargin(0.18);
  //c1->SetLogz();

  // draw to histograms using TTree::Draw()
  t->Draw("trackThetaYZ:trackThetaXZ >> allTracks_thxz_thyz");
  t->Draw("trackThetaYZ:trackThetaXZ >> lengthTracks_thxz_thyz", "trackLength>50");
  t->Draw("trackThetaYZ:trackThetaXZ >> t0Tracks_thxz_thyz", "trackLength>50 && trackIsHasT0 == 1");
  t->Draw("trackThetaYZ:trackThetaXZ >> fidVolTracks_thxz_thyz", "trackLength>50 && trackIsPassVolumeCut == 1 && trackIsHasT0 == 1");
  t->Draw("trackThetaYZ:trackThetaXZ >> selTracks_thxz_thyz", "trackIsSelected == 1");
  t->Draw("trackStartY:trackStartX_t0Corr >> selTracks_xy", "trackIsSelected == 1");
  t->Draw("trackStartX_t0Corr:trackStartZ >> selTracks_xz", "trackIsSelected == 1");
  t->Draw("trackStartY:trackStartZ        >> selTracks_yz", "trackIsSelected == 1");
  t->Draw("trackStartY_SCEcorr:trackStartX_SCEcorr >> selTracks_SCExy", "trackIsSelected == 1");
  t->Draw("trackStartX_SCEcorr:trackStartZ_SCEcorr >> selTracks_SCExz", "trackIsSelected == 1");
  t->Draw("trackStartY_SCEcorr:trackStartZ_SCEcorr >> selTracks_SCEyz", "trackIsSelected == 1");
  /*
  t->Draw("trackEndY:trackEndX_t0Corr >> selTracks_xy", "trackIsSelected == 1");
  t->Draw("trackEndX_t0Corr:trackEndZ >> selTracks_xz", "trackIsSelected == 1");
  t->Draw("trackEndY:trackEndZ        >> selTracks_yz", "trackIsSelected == 1");
  t->Draw("trackEndY_SCEcorr:trackEndX_SCEcorr >> selTracks_SCExy", "trackIsSelected == 1");
  t->Draw("trackEndX_SCEcorr:trackEndZ_SCEcorr >> selTracks_SCExz", "trackIsSelected == 1");
  t->Draw("trackEndY_SCEcorr:trackEndZ_SCEcorr >> selTracks_SCEyz", "trackIsSelected == 1");
  */
  t->Draw("trackLength >> allTracks_length");
  t->Draw("trackLength >> lengthTracks_length", "trackLength>50");
  t->Draw("trackLength >> t0Tracks_length", "trackLength>50 && trackIsHasT0 == 1");
  t->Draw("trackLength >> fidVolTracks_length", "trackLength>50 && trackIsPassVolumeCut == 1 && trackIsHasT0 == 1");
  t->Draw("trackLength >> selTracks_length", "trackIsSelected ==1");
  t2->Draw("track_length >> defTracks_length", "track_avg_trans_dist < 6");

  allTracks_thxz_thyz->Draw("colz");
  allTracks_thxz_thyz->SetContour(1000);
  c1->SaveAs("allTracks_thxz_thyz.png");
  c1->SaveAs("allTracks_thxz_thyz.pdf");

  lengthTracks_thxz_thyz->Draw("colz");
  lengthTracks_thxz_thyz->SetContour(1000);
  c1->SaveAs("lengthTracks_thxz_thyz.png");
  c1->SaveAs("lengthTracks_thxz_thyz.pdf");

  t0Tracks_thxz_thyz->Draw("colz");
  t0Tracks_thxz_thyz->SetContour(1000);
  c1->SaveAs("t0Tracks_thxz_thyz.png");
  c1->SaveAs("t0Tracks_thxz_thyz.pdf");

  fidVolTracks_thxz_thyz->Draw("colz");
  fidVolTracks_thxz_thyz->SetContour(1000);
  c1->SaveAs("fidVolTracks_thxz_thyz.png");
  c1->SaveAs("fidVolTracks_thxz_thyz.pdf");

  selTracks_thxz_thyz->Draw("colz");
  selTracks_thxz_thyz->SetContour(1000);
  c1->SaveAs("selTracks_thxz_thyz.png");
  c1->SaveAs("selTracks_thxz_thyz.pdf");

  selTracks_xy->Draw("colz");
  selTracks_xy->SetContour(1000);
  selTracks_xy->SetMinimum(-1e-9);
  c1->SaveAs("selTracks_xy.png");
  c1->SaveAs("selTracks_xy.pdf");

  selTracks_xz->Draw("colz");
  selTracks_xz->SetContour(1000);
  selTracks_xz->SetMinimum(-1e-9);
  c1->SaveAs("selTracks_xz.png");
  c1->SaveAs("selTracks_xz.pdf");

  selTracks_yz->Draw("colz");
  selTracks_yz->SetContour(1000);
  selTracks_yz->SetMinimum(-1e-9);
  c1->SaveAs("selTracks_yz.png");
  c1->SaveAs("selTracks_yz.pdf");

  selTracks_SCExy->Draw("colz");
  selTracks_SCExy->SetContour(1000);
  selTracks_SCExy->SetMinimum(-1e-9);
  c1->SaveAs("selTracks_SCExy.png");
  c1->SaveAs("selTracks_SCExy.pdf");

  selTracks_SCExz->Draw("colz");
  selTracks_SCExz->SetContour(1000);
  selTracks_SCExz->SetMinimum(-1e-9);
  c1->SaveAs("selTracks_SCExz.png");
  c1->SaveAs("selTracks_SCExz.pdf");

  selTracks_SCEyz->Draw("colz");
  selTracks_SCEyz->SetContour(1000);
  selTracks_SCEyz->SetMinimum(-1e-9);
  c1->SaveAs("selTracks_SCEyz.png");
  c1->SaveAs("selTracks_SCEyz.pdf");

  c1->SetLogy(1);

  gStyle->SetOptStat(0);
  lengthTracks_length->SetLineColor(kBlack);
  lengthTracks_length->SetLineWidth(2);
  lengthTracks_length->Draw();
  t0Tracks_length->SetLineColor(kPTRed);
  t0Tracks_length->SetLineWidth(2);
  t0Tracks_length->Draw("same");
  fidVolTracks_length->SetLineColor(kPTOrange);
  fidVolTracks_length->SetLineWidth(2);
  fidVolTracks_length->Draw("same");
  selTracks_length->SetLineColor(kPTLightBlue);
  selTracks_length->SetLineWidth(2);
  selTracks_length->Draw("same");
  defTracks_length->SetLineColor(kPTDarkBlue);
  defTracks_length->SetLineWidth(2);
  defTracks_length->Draw("same");
  TLegend *leg = new TLegend(0.50, 0.60, 0.75, 0.9);
  gStyle->SetLegendBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(lengthTracks_length, "Tracks > 50 cm", "l");
  leg->AddEntry(t0Tracks_length, "t_{0}-tagged Tracks", "l");
  leg->AddEntry(fidVolTracks_length, "Throughgoing Tracks", "l");
  leg->AddEntry(selTracks_length, "Angle Cuts", "l");
  leg->AddEntry(defTracks_length, "Track Deflection Cut", "l");
  leg->Draw("same");
  c1->SaveAs("trackLengths.png");
  c1->SaveAs("trackLengths.pdf");

  c1->SetLogz(0);
  c1->SetLogy(0);
  
  TH2D* t0Tracks_thxz_thyz_ratio 
    = (TH2D*)t0Tracks_thxz_thyz->Clone("t0Tracks_thxz_thyz_ratio");
  t0Tracks_thxz_thyz_ratio->Divide(allTracks_thxz_thyz);
  t0Tracks_thxz_thyz_ratio->Draw("colz");
  t0Tracks_thxz_thyz_ratio->SetContour(1000);
  c1->SaveAs("t0Tracks_thxz_thyz_ratio.png");
  c1->SaveAs("t0Tracks_thxz_thyz_ratio.pdf");
  
  TH2D* selTracks_thxz_thyz_ratio 
    = (TH2D*)selTracks_thxz_thyz->Clone("selTracks_thxz_thyz_ratio");
  selTracks_thxz_thyz_ratio->Divide(allTracks_thxz_thyz);
  selTracks_thxz_thyz_ratio->Draw("colz");
  selTracks_thxz_thyz_ratio->SetContour(1000);
  c1->SaveAs("selTracks_thxz_thyz_ratio.png");
  c1->SaveAs("selTracks_thxz_thyz_ratio.pdf");

  std::cout << " n total tracks:              " <<  allTracks_thxz_thyz->GetEntries() << std::endl;
  std::cout << " Tracks > 50 cm:              " <<  lengthTracks_thxz_thyz->GetEntries() << std::endl;
  std::cout << " n t0-tagged tracks > 50 cm:  " <<  t0Tracks_thxz_thyz->GetEntries()  << std::endl;
  std::cout << " % t0-tagged tracks > 50 cm:  " <<  100 * t0Tracks_thxz_thyz->GetEntries()/
                                                          lengthTracks_thxz_thyz->GetEntries() << std::endl;
  std::cout << " fidVol tracks:               " <<  fidVolTracks_thxz_thyz->GetEntries()  << std::endl;
  std::cout << " \% fidVol/Total:             " <<  100 * fidVolTracks_thxz_thyz->GetEntries()/  
                                                          allTracks_thxz_thyz->GetEntries() << std::endl;
  std::cout << " \% fidVol/t0:                " <<  100 * fidVolTracks_thxz_thyz->GetEntries()/  
                                                          t0Tracks_thxz_thyz->GetEntries() << std::endl;
  std::cout << " n selected tracks:           " <<  selTracks_thxz_thyz->GetEntries() << std::endl;
  std::cout << " \% selected tracks:          " <<  100 * selTracks_thxz_thyz->GetEntries()/
                                                    allTracks_thxz_thyz->GetEntries() << std::endl;
  std::cout << " tracks pass deflection:      " << defTracks_length->GetEntries() << std::endl;

}
