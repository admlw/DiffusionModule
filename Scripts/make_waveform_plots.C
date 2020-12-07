void make_waveform_plots() {
  TFile *fin_single = new TFile("WirePlot/mp.root", "READ");
  if (!fin_single) {
    std::cout << "ERROR Bad input file" << std::endl;
    return;
  }
  
  TH1D *h_single = (TH1D*)fin_single->Get("Event_1478/Event_1478_TimeWfm_channel6945");
  if (!h_single) {
    std::cout << "ERROR Could not read histogram" << std::endl;
    return;
  }

  TCanvas *c1 = new TCanvas;
  h_single->SetLineColor(kAzure+2);
  h_single->SetTitle("");
  h_single->GetXaxis()->SetRangeUser(3655, 3715);
  h_single->GetYaxis()->SetRangeUser(-0.5, 10);
  h_single->GetXaxis()->SetTitle("Time (ticks)");
  h_single->GetYaxis()->SetTitle("Arb. Units");
  h_single->Draw();

  TPaveText *pt1 = new TPaveText(0.15, 0.7, 0.45, 0.9, "NDC");
  TString rms_string_single = Form("RMS: %1.2f", h_single->GetRMS() );
  pt1->SetTextAlign(11);
  pt1->SetFillStyle(0);
  pt1->SetLineWidth(0);
  pt1->SetTextFont(43);
  pt1->SetTextSize(24);
  pt1->AddText(rms_string_single);
  pt1->Draw("same");

  c1->SaveAs("waveform_plot_single.pdf", "PDF");

  TString summed_dir = "/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionAna/";
  TFile *fin_summed = new TFile(summed_dir+"diffmod_run3_crt.root", "READ");

  TCanvas *c2 = new TCanvas;
  TH1D *h_summed = (TH1D*)fin_summed->Get("DiffusionModule/plane2/summed_waveform_bin_15_plane2");
  h_summed->SetLineColor(kAzure+2);
  h_summed->SetTitle("");
  h_summed->GetXaxis()->SetRangeUser(3620, 3680);
  h_summed->GetXaxis()->SetTitle("Time (ticks)");
  h_summed->GetYaxis()->SetTitle("Arb. Units");
  // Summed waveform from file already has a gaussian fit; don't draw that
  h_summed->GetFunction("gaus")->SetBit(TF1::kNotDraw);
  h_summed->Draw();

  TPaveText *pt2 = new TPaveText(0.15, 0.7, 0.45, 0.9, "NDC");
  TString rms_string_summed = Form("RMS: %1.2f", h_summed->GetRMS() );
  pt2->SetTextAlign(11);
  pt2->SetFillStyle(0);
  pt2->SetLineWidth(0);
  pt2->SetTextFont(43);
  pt2->SetTextSize(24);
  pt2->AddText(rms_string_summed);
  pt2->Draw("same");

  c2->SaveAs("waveform_plot_summed.pdf", "PDF");


}
