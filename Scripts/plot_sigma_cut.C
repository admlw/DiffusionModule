void plot_sigma_cut() {
  TString dir = "/uboone/app/users/amogan/mcc9_diffusion/workdir/"
  TFile *fin = new TFile(dir+"sigma_map_CV.root", "READ");
  TH2D *h_precut = (TH2D*)fin->Get(
}
