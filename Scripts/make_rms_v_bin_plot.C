// Script for making a plot of average waveform rms
// in each drift bin. 
// Now with all three planes!
// Takes sigma_map.root file as input

void make_rms_v_bin_plot() {
  // Setup
  const int n = 25;
  _file0->cd();
  TString dir0 = "DiffusionModule/plane0/";
  TString dir1 = "DiffusionModule/plane1/";
  TString dir2 = "DiffusionModule/plane2/";
  TH1D *hPlane0;
  TH1D *hPlane1;
  TH1D *hPlane2;
  TString histNamePlane0;
  TString histNamePlane1;
  TString histNamePlane2;

  double binVals[n]; // Bin no. doesn't change between planes (duh)
  double rmsValsPlane0[n];
  double rmsValsPlane1[n];
  double rmsValsPlane2[n];

  // x-errors will be half the width of a drift time bin
  double binValErrs[n];
  double rmsValsErrsPlane0[n];
  double rmsValsErrsPlane1[n];
  double rmsValsErrsPlane2[n];

  // Not sure how to do errors on an RMS, so just give it some small value
  std::fill_n(binValErrs, n, 1e-7); // TODO: Update with actual values
  std::fill_n(rmsValsErrsPlane0, n, 1e-7);
  std::fill_n(rmsValsErrsPlane1, n, 1e-7);
  std::fill_n(rmsValsErrsPlane2, n, 1e-7);

  gStyle->SetOptStat(0);
  for (int i = 0; i < n; i++) {
    std::cout << "------------------------------" << std::endl;
    std::cout << "Starting bin " << i << std::endl;
    histNamePlane0 = Form("h_sigma_%i_plane0", i);
    histNamePlane1 = Form("h_sigma_%i_plane1", i);
    histNamePlane2 = Form("h_sigma_%i_plane2", i);
    hPlane0 = (TH1D*)_file0->Get(dir0+histNamePlane0);
    hPlane1 = (TH1D*)_file0->Get(dir1+histNamePlane1);
    hPlane2 = (TH1D*)_file0->Get(dir2+histNamePlane2);
    binVals[i] = i;
    rmsValsPlane0[i] = hPlane0->GetRMS();
    rmsValsPlane1[i] = hPlane1->GetRMS();
    rmsValsPlane2[i] = hPlane2->GetRMS();
    std::cout << "Plane 0 val " << rmsValsPlane0[i] << std::endl;
    std::cout << "Plane 1 val " << rmsValsPlane1[i] << std::endl;
    std::cout << "Plane 2 val " << rmsValsPlane2[i] << std::endl;
    std::cout << "------------------------------" << std::endl;
  }

  TCanvas *c = new TCanvas("c", "c", 750, 550);
  c->cd();

  TGraphErrors *g0 = new TGraphErrors(n, binVals, rmsValsPlane0, binValErrs, rmsValsErrsPlane0);
  g0->SetTitle("");
  g0->GetXaxis()->SetTitle("Bin No.");
  g0->GetYaxis()->SetTitle("RMS of #sigma_{t}^{2} Dist.");
  g0->SetLineColor(kRed+1);
  g0->SetLineWidth(2);
  g0->Draw();

  /*
  TGraphErrors *g1 = new TGraphErrors(n, binVals, rmsValsPlane1, binValErrs, rmsValsErrsPlane1);
  g1->SetTitle("");
  g1->SetLineColor(kGreen+1);
  g1->SetLineWidth(2);
  g1->Draw("same");

  TGraphErrors *g2 = new TGraphErrors(n, binVals, rmsValsPlane2, binValErrs, rmsValsErrsPlane2);
  g2->SetTitle("");
  g2->SetLineColor(kAzure+1);
  g2->SetLineWidth(2);
  g2->Draw("same");
  */

  TLegend *l = new TLegend(0.6, 0.7, 0.9, 0.9);
  l->AddEntry(g0, "Plane 0", "l");
  //l->AddEntry(g1, "Plane 1", "l");
  //l->AddEntry(g2, "Plane 2", "l");
  l->Draw("same");

  c->SaveAs("rms_v_bin_plot.png", "PNG");

}
