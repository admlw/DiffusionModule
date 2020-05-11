void make_diffusion_coefficient_plot() {
  // Velocity from uBooNE measurement
  double v   = 0.1098*0.1098; // v^2, in cm^2/us^2
 
  // Diffusion coefficients from BNL at 89 K and 273 V/cm, 
  // https://lar.bnl.gov/properties/
  double D_L = 6.4e-6; // cm^2/microsecond
  double D_T = 9.7e-6;

  // D_T and D_L treated in same way; see slide 8 at
  // https://indico.hep.manchester.ac.uk/getFile.py/access?sessionId=18&resId=0&materialId=0&confId=5346
  
  // Transverse
  TF1 *f_dt = new TF1("f_dt", "sqrt(2*[0]*x/[1])", 0, 2300); // Note x here is drift time in us
  f_dt->SetParameter(0, D_T);
  f_dt->SetParameter(1, v);
  f_dt->SetParName  (0, "D_T");
  f_dt->SetParName  (1, "Velocity");
  f_dt->SetLineColor(kRed+1);
  f_dt->SetLineStyle(10);
  f_dt->GetXaxis()->SetTitle("Drift Time (#mus)");
  f_dt->GetXaxis()->SetNdivisions(5);
  f_dt->GetYaxis()->SetTitle("#sigma (#mus)");
  f_dt->GetYaxis()->SetTitleOffset(1.0);
  f_dt->SetTitle("");

  // Longitudinal
  TF1 *f_dl = new TF1("f_dl", "sqrt(2*[0]*x/[1])", 0, 2300); // Note x here is drift time in us
  f_dl->SetParameter(0, D_L);
  f_dl->SetParameter(1, v);
  f_dl->SetParName  (0, "D_L");
  f_dl->SetParName  (1, "Velocity");
  f_dl->SetLineColor(kAzure+2);
  f_dl->SetLineStyle(2);

  TCanvas *c = new TCanvas;
  f_dt->Draw();
  f_dl->Draw("same");

  TLegend *l = new TLegend(0.17, 0.75, 0.44, 0.9);
  l->AddEntry(f_dt, "D_{T} = 9.7 cm^{2}/s", "l");
  l->AddEntry(f_dl, "D_{L} = 6.4 cm^{2}/s", "l");
  l->Draw("same");

  TPaveText *pt = new TPaveText(0.63, 0.13, 0.83, 0.33, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(22); // Was 132
  pt->AddText("E = 273 V/cm");
  pt->Draw("same");

  //c->SaveAs("diffusion_v_drift_time.png", "PNG");
  c->SaveAs("sigma_v_drift_time.pdf", "PDF");


}



