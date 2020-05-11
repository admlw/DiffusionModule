void make_diffusion_selection_validation_plots(){

  std::cout << "Starting script" << std::endl;
  TTree* t = (TTree*)_file0->Get("diffsel/difffiltertree");
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

  // save plots
  TCanvas* c1 = new TCanvas();
  c1->SetRightMargin(0.18);
  c1->SetLogz();

  // draw to histograms using TTree::Draw()
  t->Draw("trackThetaYZ:trackThetaXZ >> allTracks_thxz_thyz");
  t->Draw("trackThetaYZ:trackThetaXZ >> lengthTracks_thxz_thyz", "trackLength>50");
  t->Draw("trackThetaYZ:trackThetaXZ >> t0Tracks_thxz_thyz", "trackLength>50 && trackIsHasT0 == 1");
  t->Draw("trackThetaYZ:trackThetaXZ >> fidVolTracks_thxz_thyz", "trackLength>50 && trackIsPassVolumeCut == 1 && trackIsHasT0 == 1");
  t->Draw("trackThetaYZ:trackThetaXZ >> selTracks_thxz_thyz", "trackIsSelected == 1");
  t->Draw("trackLength >> allTracks_length");
  t->Draw("trackLength >> lengthTracks_length", "trackLength>50");
  t->Draw("trackLength >> t0Tracks_length", "trackLength>50 && trackIsHasT0 == 1");
  t->Draw("trackLength >> fidVolTracks_length", "trackLength>50 && trackIsPassVolumeCut == 1 && trackIsHasT0 == 1");
  t->Draw("trackLength >> selTracks_length", "trackIsSelected ==1");

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

  c1->SetLogy(1);

  //allTracks_length->SetLineColor(kBlack);
  //allTracks_length->SetLineWidth(2);
  //allTracks_length->Draw();
  gStyle->SetOptStat(0);
  lengthTracks_length->SetLineColor(kBlack);
  lengthTracks_length->SetLineWidth(2);
  lengthTracks_length->Draw();
  t0Tracks_length->SetLineColor(kAzure+1);
  t0Tracks_length->SetLineWidth(2);
  t0Tracks_length->Draw("same");
  fidVolTracks_length->SetLineColor(kRed+1);
  fidVolTracks_length->SetLineWidth(2);
  fidVolTracks_length->Draw("same");
  selTracks_length->SetLineColor(kGreen+1);
  selTracks_length->SetLineWidth(2);
  selTracks_length->Draw("same");
  TLegend *leg = new TLegend(0.55, 0.65, 0.8, 0.9);
  gStyle->SetLegendBorderSize(0);
  leg->AddEntry(lengthTracks_length, "Tracks > 50 cm");
  leg->AddEntry(t0Tracks_length, "T0-tagged Tracks");
  leg->AddEntry(fidVolTracks_length, "Tracks in FV");
  leg->AddEntry(selTracks_length, "Angle Cuts");
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

}
