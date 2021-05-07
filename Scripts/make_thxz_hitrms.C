void make_thxz_hitrms(){

  //TFile* f2 = new TFile("/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionAna/diffmod_run3_crt_Aug2020_newFV_bugFix.root", "read");
  TFile* f2 = new TFile("/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionAna/diffmod_prod_muminus_200k.root", "read");
  TTree* t2 = (TTree*)f2->Get("DiffusionModule/difftree");

  TH2D* h = new TH2D("h", ";Track #theta_{XZ}; Hit RMS", 10, 0, 10, 100, 0, 4);

  t2->Draw("hit_rms:track_theta_xz >> h");
  h->Draw("colz");
 

}
