void make_hit_peak_time_plots() {

  //TString dir = "/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionAna/";
  //TFile *fin = new TFile(dir+"diffmod_run3_crt_fidVol.root", "READ");
  TString dir = "/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionSigmamap/";
  TFile *fin = new TFile(dir+"sigma_map_run3_crt_throughgoing_SCE.root", "READ");
  //TFile *fin = new TFile(dir+"diffmod_run3_crt.root", "READ");
  TTree *t = (TTree*)fin->Get("DiffusionModule/difftree");
  if (!t || !fin) {
    std::cout << "ERROR" << std::endl;
  }

  TH1D *h1 = new TH1D("h1", "", 100, 0, 7000);
  TH1D *h2 = new TH1D("h2", "", 100, 0, 7000);

  TCanvas *c = new TCanvas();
  
  t->Draw("hit_peak_time        >> h1");
  t->Draw("hit_peak_time_t0corr >> h2");

  h1->SetLineColor(kAzure+1);
  h1->SetTitle("");

  h2->SetLineColor(kRed+1);
  h2->SetTitle("");
  h2->GetXaxis()->SetTitle("Hit Peak Time (Ticks)");
  //h2->GetYaxis()->SetRangeUser(0, 550000);
  //h2->GetYaxis()->SetRangeUser(0, 1000000);
  h2->Draw();
  h1->Draw("same");

  TLegend *leg = new TLegend(0.15, 0.7, 0.45, 0.9);
  leg->AddEntry(h1, "Non-t_{0}-Corrected", "l");
  leg->AddEntry(h2, "t_{0}-Corrected", "l");
  leg->Draw("same");

  c->SaveAs("hit_peak_time.pdf", "PDF"); 

}
