void nSigmasPlot(){

  TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
  c1->SetLeftMargin(0.22);
  c1->SetRightMargin(0.01);
  c1->SetLogx();

  double nSigma[4] = {0.25, 0.5, 1, 100};
  double nSigmaErr[4] = {0, 0, 0, 0};
  double diff[4] = {6.294/6.4, 6.290/6.4, 6.299/6.4, 6.442/6.4};
  double diffErr[4] = {0.0891/6.4, 0.0882/6.4, 0.0880/6.4, 0.0898/6.4};

  TGraphErrors *gr = new TGraphErrors(4, nSigma, diff, nSigmaErr, diffErr);
  gr->GetXaxis()->SetRangeUser(0.01, 110);
  gr->GetYaxis()->SetRangeUser(0.9, 1.1);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.4);
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("#sigma Cut Value");
  gr->GetYaxis()->SetTitle("#frac{Measured D_{L}}{Input D_{L}}");
  gr->GetYaxis()->SetTitleOffset(1.8);
  gr->Draw("ap");

  TF1* f = new TF1("f", "1.0", 0.01, 110);
  f->SetLineColor(kGray);
  f->Draw("same");
  
  gr->Draw("psame");

}
