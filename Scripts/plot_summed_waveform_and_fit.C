void plot_summed_waveform_and_fit(){

  TFile* f = new TFile("/pnfs/uboone/persistent/users/alister1/diffusion_crt_data_forPaper/ana/diffusion_ana_crt_run3_data_max.25.root", "read");
  TH1D* h = (TH1D*)f->Get("DiffusionModule/plane2/summed_waveform_bin_1_plane2");
  TF1*  f2 = h->GetFunction("gaus");
  h->ls();
  h->Draw();
  //h->Fit("gaus");
  f2->Draw("same");

}
