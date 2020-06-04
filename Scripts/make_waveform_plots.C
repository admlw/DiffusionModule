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
  c1->SaveAs("waveform_plot_single.pdf", "PDF");

  TString summed_dir = "/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionAna/";
  TFile *fin_summed = new TFile(summed_dir+"diffmod_run3_crt.root", "READ");

  TCanvas *c2 = new TCanvas;
  TH1D *h_summed = (TH1D*)fin_summed->Get("DiffusionModule/plane2/summed_waveform_bin_19_plane2");
  h_summed->SetLineColor(kAzure+2);
  h_summed->SetTitle("");
  h_summed->GetXaxis()->SetRangeUser(4355, 4415);
  h_summed->GetXaxis()->SetTitle("Time (ticks)");
  h_summed->GetYaxis()->SetTitle("Arb. Units");
  // Summed waveform from file already has a gaussian fit; don't draw that
  h_summed->GetFunction("gaus")->SetBit(TF1::kNotDraw);
  h_summed->Draw();
  c2->SaveAs("waveform_plot_summed.pdf", "PDF");


}
