// produces world data and theory summaries for
// - epsilon_l
// - D_L
//
// Theory:
//  - Atrazhev-Timoshkin
//  - BNL parametrization
//
// Datasets:
//  - MicroBooNE
//  - ICARUS

//.............................................................................
Double_t bnl_mu_param(double ef){
  float a0 = 551.6;
  float a1 = 7953.7;
  float a2 = 4440.43;
  float a3 = 4.29;
  float a4 = 43.63;
  float a5 = 0.2053;

  ef=ef/1000;

  float num = a0 + a1*ef + a2*std::pow(ef, 1.5) + a3*std::pow(ef, 2.5);
  float den = 1 + (a1/a0)*ef + a4*std::pow(ef, 2) + a5*std::pow(ef,3);

  return (num/den) * std::pow(89./87., -1.5);

}

//.............................................................................
Double_t bnl_el_param(double ef){
  float b0 = 0.0075; 
  float b1 = 742.9;
  float b2 = 3269.6;
  float b3 = 31678.2;
  float T  = 89;
  float T1 = 87;

  ef=ef/1000;

  return ((b0 + b1*ef + b2*ef*ef)   /
         (1+(b1/b0)*ef + b3*ef*ef)) *
         (T/T1);
}

//.............................................................................
// the atrazhev parameterisation comes from running the fit_atrazhev script
Double_t atrazhev_el_param(double ef){

  double offs = 89 * 0.000086173;
  double p0    = -0.000140438;
  double p1    = 1.38788e-05;
  double p2    = 2.19193e-09;
  double p3    = -2.69324e-13;
  double p4    = 1.14616e-17;

  return offs + p0 +
         p1 * ef   +
         p2 * std::pow(ef, 2) +
         p3 * std::pow(ef, 3) +
         p4 * std::pow(ef, 4);

}

//.............................................................................
Double_t bnl_dl_param(double ef){
  return bnl_mu_param(ef) * bnl_el_param(ef);
}

//.............................................................................
Double_t atrazhev_dl_param(double ef){
  return bnl_mu_param(ef) * atrazhev_el_param(ef);
}

//.............................................................................
// main
void make_theory_world_data_plots(){

  gStyle->SetOptStat(0);

  //.............................................
  // datasets
  //.............................................

  // uboone

  Double_t dl_ubx[1]  = {273.9};
  Double_t dl_ubxu[1] = {273.9*1.15-(273.9)};
  Double_t dl_ubxd[1] = {273.9*1.075-(273.9)};
  Double_t dl_uby[1]  = {3.74};
  Double_t dl_ubyu[1]  = {3.74*1.076 - 3.74};
  Double_t dl_ubyd[1]  = {3.74*1.077 - 3.74};

  TGraphAsymmErrors* dl_ub = new TGraphAsymmErrors(1, dl_ubx, dl_uby, dl_ubxd, dl_ubxu, dl_ubyd, dl_ubyu);
  dl_ub->SetMarkerStyle(20);
  dl_ub->SetMarkerColor(kAzure+1);
  dl_ub->SetLineWidth(2);
  dl_ub->SetLineColor(kAzure+1);

  Double_t el_ubx[1]   = {273.9};
  Double_t el_ubxu[1]  = {273.9*1.15-(273.9)};
  Double_t el_ubxd[1]  = {273.9*1.075-(273.9)};
  Double_t el_uby[1]   = {3.74/bnl_mu_param(273.9)};
  Double_t el_ubyu[1]  = {(3.74/bnl_mu_param(273.9))*1.076 - 3.74/bnl_mu_param(273.9)};
  Double_t el_ubyd[1]  = {(3.74/bnl_mu_param(273.9))*1.077 - 3.74/bnl_mu_param(273.9)};

  TGraphAsymmErrors* el_ub = new TGraphAsymmErrors(1, el_ubx, el_uby, el_ubxd, el_ubxu, el_ubyd, el_ubyu);
  el_ub->SetMarkerStyle(20);
  el_ub->SetMarkerColor(kAzure+1);
  el_ub->SetLineWidth(2);
  el_ub->SetLineColor(kAzure+1);

  // icarus
  Double_t dl_icarusx[1]  = {225};
  Double_t dl_icarusxu[1] = {350-225};
  Double_t dl_icarusxd[1] = {225-100};
  Double_t dl_icarusy[1]  = {4.8};
  Double_t dl_icarusyu[1]  = {0.2};
  Double_t dl_icarusyd[1]  = {0.2};

  TGraphAsymmErrors* dl_icarus = new TGraphAsymmErrors(1, dl_icarusx, dl_icarusy, dl_icarusxd, dl_icarusxu, dl_icarusyd, dl_icarusyu);
  dl_icarus->SetMarkerStyle(20);
  dl_icarus->SetMarkerColor(kPink+1);
  dl_icarus->SetLineWidth(2);
  dl_icarus->SetLineColor(kPink+1);

  Double_t el_icarusx[1]  = {225};
  Double_t el_icarusxu[1] = {350-225};
  Double_t el_icarusxd[1] = {225-100};
  Double_t el_icarusy[1]  = {4.8/bnl_mu_param(225)};
  Double_t el_icarusyu[1]  = {0.2/bnl_mu_param(225)};
  Double_t el_icarusyd[1]  = {0.2/bnl_mu_param(225)};

  TGraphAsymmErrors* el_icarus = new TGraphAsymmErrors(1, el_icarusx, el_icarusy, el_icarusxd, el_icarusxu, el_icarusyd, el_icarusyu);
  el_icarus->SetMarkerStyle(20);
  el_icarus->SetMarkerColor(kPink+1);
  el_icarus->SetLineWidth(2);
  el_icarus->SetLineColor(kPink+1);

  //.............................................
  // electron energy
  //.............................................

  TCanvas* c1 = new TCanvas("c1", "", 500, 500);
  c1->SetGridy();
  c1->SetGridx();
  c1->SetLogx();
  c1->SetLogy();
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.12);

  TH2D* el_bg = new TH2D("el_bg", ";E (V/cm);Electron Energy, #epsilon_{L} (eV)", 100, 80, 10000, 100, 5e-3, 3.1e-1);
  el_bg->GetXaxis()->CenterTitle();
  el_bg->GetYaxis()->CenterTitle();
  el_bg->Draw();

  TF1* bnl_el = new TF1("bnl_el", "bnl_el_param(x)", 10e-2, 1500);
  bnl_el->SetLineStyle(5);
  bnl_el->SetLineColor(kOrange+3);
  bnl_el->SetLineWidth(2);
  bnl_el->Draw("l same");
  
  TF1* atrazhev_el = new TF1("atrazhev_el", "atrazhev_el_param(x)", 10e-2, 10000);
  atrazhev_el->Draw("l same");
  
  el_icarus->Draw("p same");
  el_ub->Draw("p same");

  TLegend* leg = new TLegend(0.15, 0.70, 0.75, 0.85);
  leg->AddEntry(atrazhev_el, "Atrazhev-Timoshkin [10.1109/94.689434]", "l");
  leg->AddEntry(bnl_el     , "BNL [10.1016/j.nima.2016.01.094]", "l");
  leg->AddEntry(el_icarus  , "ICARUS [10.1016/0168-9002(94)90996-2]", "pe");
  leg->AddEntry(el_ub      , "MicroBooNE Data", "pe");
  leg->Draw("same");

  c1->SaveAs("el_summary.pdf");

  //.............................................
  // diffusion plots
  //.............................................

  TCanvas* c2 = new TCanvas("c2", "", 500, 500);
  c2->SetGridy();
  c2->SetGridx();
  c2->SetLogx();
  c2->SetLeftMargin(0.12);
  c2->SetBottomMargin(0.12);

  TH2D* dl_bg = new TH2D("dl_bg", ";E (V/cm);D_{L}", 100, 80, 10000, 100, 0, 10);
  dl_bg->GetXaxis()->CenterTitle();
  dl_bg->GetYaxis()->CenterTitle();
  dl_bg->Draw();

  TF1* bnl_dl = new TF1("bnl_dl", "bnl_dl_param(x)", 10e-2, 1500);
  bnl_dl->SetLineStyle(5);
  bnl_dl->SetLineColor(kOrange+3);
  bnl_dl->SetLineWidth(2);
  bnl_dl->Draw("l same");
  
  TF1* atrazhev_dl = new TF1("atrazhev_dl", "atrazhev_dl_param(x)", 10e-2, 10000);
  atrazhev_dl->Draw("l same");

  dl_icarus->Draw("same p");
  dl_ub->Draw("same p");
  leg->Draw("same");

  c2->SaveAs("dl_summary.pdf");
}
