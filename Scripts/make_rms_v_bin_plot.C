void make_rms_v_bin_plot() {
  // Setup
  const int n = 25;
  _file0->cd();
  TString dir = "DiffusionModule/";
  TH1D *h;
  TString histName;
  double binVals[n];
  double rmsVals[n];
  // x-errors will be half the width of a drift time bin
  double binValErrs[n];
  double rmsValsErrs[n];
  // Not sure how to do errors on an RMS, so just give it some small value
  std::fill_n(binValErrs, n, 1e-7); // TODO: Update with actual values
  std::fill_n(rmsValsErrs, n, 1e-7);

  gStyle->SetOptStat(0);
  for (int i = 0; i < n; i++) {
    histName = Form("h_sigma_%i", i);
    h = (TH1D*)_file0->Get(dir+histName);
    binVals[i] = i;
    rmsVals[i] = h->GetRMS();
  }

  TCanvas *c = new TCanvas("c", "c", 750, 550);
  c->cd();
  TGraphErrors *g = new TGraphErrors(n, binVals, rmsVals, binValErrs, rmsValsErrs);
  g->SetTitle("");
  g->GetXaxis()->SetTitle("Bin No.");
  g->GetYaxis()->SetTitle("RMS of #sigma_{t}^{2} Dist.");
  g->Draw();
  c->SaveAs("rms_v_bin_plot.png", "PNG");

}
