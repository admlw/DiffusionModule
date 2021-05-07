#include "StylePlots.h"

void make_diffusion_t0tagging_deepdive(){

  TFile* f2 = new TFile("/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionAna/diffmod_run3_crt_Aug2020_newFV_bugFix.root", "read");
  TTree* t2 = (TTree*)f2->Get("DiffusionModule/difftree");

  TFile* outfile = new TFile("hit_width_dists.root", "recreate");

  TH1D* hu0  = new TH1D("hu0" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu1  = new TH1D("hu1" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu2  = new TH1D("hu2" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu3  = new TH1D("hu3" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu4  = new TH1D("hu4" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu5  = new TH1D("hu5" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu6  = new TH1D("hu6" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu7  = new TH1D("hu7" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu8  = new TH1D("hu8" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu9  = new TH1D("hu9" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu10  = new TH1D("hu10" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu11  = new TH1D("hu11" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu12  = new TH1D("hu12" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu13  = new TH1D("hu13" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu14  = new TH1D("hu14" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu15  = new TH1D("hu15" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu16  = new TH1D("hu16" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu17  = new TH1D("hu17" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu18  = new TH1D("hu18" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu19  = new TH1D("hu19" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu20  = new TH1D("hu20" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu21  = new TH1D("hu21" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu22  = new TH1D("hu22" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu23  = new TH1D("hu23" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu24  = new TH1D("hu24" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);

  TH1D* hv0  = new TH1D("hv0" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv1  = new TH1D("hv1" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv2  = new TH1D("hv2" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv3  = new TH1D("hv3" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv4  = new TH1D("hv4" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv5  = new TH1D("hv5" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv6  = new TH1D("hv6" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv7  = new TH1D("hv7" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv8  = new TH1D("hv8" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv9  = new TH1D("hv9" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv10  = new TH1D("hv10" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv11  = new TH1D("hv11" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv12  = new TH1D("hv12" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv13  = new TH1D("hv13" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv14  = new TH1D("hv14" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv15  = new TH1D("hv15" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv16  = new TH1D("hv16" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv17  = new TH1D("hv17" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv18  = new TH1D("hv18" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv19  = new TH1D("hv19" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv20  = new TH1D("hv20" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv21  = new TH1D("hv21" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv22  = new TH1D("hv22" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv23  = new TH1D("hv23" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv24  = new TH1D("hv24" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);

  TH1D* hy0  = new TH1D("hy0" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy1  = new TH1D("hy1" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy2  = new TH1D("hy2" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy3  = new TH1D("hy3" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy4  = new TH1D("hy4" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy5  = new TH1D("hy5" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy6  = new TH1D("hy6" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy7  = new TH1D("hy7" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy8  = new TH1D("hy8" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy9  = new TH1D("hy9" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy10  = new TH1D("hy10" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy11  = new TH1D("hy11" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy12  = new TH1D("hy12" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy13  = new TH1D("hy13" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy14  = new TH1D("hy14" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy15  = new TH1D("hy15" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy16  = new TH1D("hy16" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy17  = new TH1D("hy17" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy18  = new TH1D("hy18" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy19  = new TH1D("hy19" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy20  = new TH1D("hy20" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy21  = new TH1D("hy21" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy22  = new TH1D("hy22" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy23  = new TH1D("hy23" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy24  = new TH1D("hy24" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);

  t2->Draw("(hit_rms/2.) >> hu0"  , "hit_view == 0 && wvfm_bin_no == 0");
  t2->Draw("(hit_rms/2.) >> hu1"  , "hit_view == 0 && wvfm_bin_no == 1");
  t2->Draw("(hit_rms/2.) >> hu2"  , "hit_view == 0 && wvfm_bin_no == 2");
  t2->Draw("(hit_rms/2.) >> hu3"  , "hit_view == 0 && wvfm_bin_no == 3");
  t2->Draw("(hit_rms/2.) >> hu4"  , "hit_view == 0 && wvfm_bin_no == 4");
  t2->Draw("(hit_rms/2.) >> hu5"  , "hit_view == 0 && wvfm_bin_no == 5");
  t2->Draw("(hit_rms/2.) >> hu6"  , "hit_view == 0 && wvfm_bin_no == 6");
  t2->Draw("(hit_rms/2.) >> hu7"  , "hit_view == 0 && wvfm_bin_no == 7");
  t2->Draw("(hit_rms/2.) >> hu8"  , "hit_view == 0 && wvfm_bin_no == 8");
  t2->Draw("(hit_rms/2.) >> hu9"  , "hit_view == 0 && wvfm_bin_no == 9");
  t2->Draw("(hit_rms/2.) >> hu10" , "hit_view == 0 && wvfm_bin_no == 10");
  t2->Draw("(hit_rms/2.) >> hu11" , "hit_view == 0 && wvfm_bin_no == 11");
  t2->Draw("(hit_rms/2.) >> hu12" , "hit_view == 0 && wvfm_bin_no == 12");
  t2->Draw("(hit_rms/2.) >> hu13" , "hit_view == 0 && wvfm_bin_no == 13");
  t2->Draw("(hit_rms/2.) >> hu14" , "hit_view == 0 && wvfm_bin_no == 14");
  t2->Draw("(hit_rms/2.) >> hu15" , "hit_view == 0 && wvfm_bin_no == 15");
  t2->Draw("(hit_rms/2.) >> hu16" , "hit_view == 0 && wvfm_bin_no == 16");
  t2->Draw("(hit_rms/2.) >> hu17" , "hit_view == 0 && wvfm_bin_no == 17");
  t2->Draw("(hit_rms/2.) >> hu18" , "hit_view == 0 && wvfm_bin_no == 18");
  t2->Draw("(hit_rms/2.) >> hu19" , "hit_view == 0 && wvfm_bin_no == 19");
  t2->Draw("(hit_rms/2.) >> hu20" , "hit_view == 0 && wvfm_bin_no == 20");
  t2->Draw("(hit_rms/2.) >> hu21" , "hit_view == 0 && wvfm_bin_no == 21");
  t2->Draw("(hit_rms/2.) >> hu22" , "hit_view == 0 && wvfm_bin_no == 22");
  t2->Draw("(hit_rms/2.) >> hu23" , "hit_view == 0 && wvfm_bin_no == 23");
  t2->Draw("(hit_rms/2.) >> hu24" , "hit_view == 0 && wvfm_bin_no == 24");

  t2->Draw("(hit_rms/2.) >> hv0"  , "hit_view == 1 && wvfm_bin_no == 0");
  t2->Draw("(hit_rms/2.) >> hv1"  , "hit_view == 1 && wvfm_bin_no == 1");
  t2->Draw("(hit_rms/2.) >> hv2"  , "hit_view == 1 && wvfm_bin_no == 2");
  t2->Draw("(hit_rms/2.) >> hv3"  , "hit_view == 1 && wvfm_bin_no == 3");
  t2->Draw("(hit_rms/2.) >> hv4"  , "hit_view == 1 && wvfm_bin_no == 4");
  t2->Draw("(hit_rms/2.) >> hv5"  , "hit_view == 1 && wvfm_bin_no == 5");
  t2->Draw("(hit_rms/2.) >> hv6"  , "hit_view == 1 && wvfm_bin_no == 6");
  t2->Draw("(hit_rms/2.) >> hv7"  , "hit_view == 1 && wvfm_bin_no == 7");
  t2->Draw("(hit_rms/2.) >> hv8"  , "hit_view == 1 && wvfm_bin_no == 8");
  t2->Draw("(hit_rms/2.) >> hv9"  , "hit_view == 1 && wvfm_bin_no == 9");
  t2->Draw("(hit_rms/2.) >> hv10" , "hit_view == 1 && wvfm_bin_no == 10");
  t2->Draw("(hit_rms/2.) >> hv11" , "hit_view == 1 && wvfm_bin_no == 11");
  t2->Draw("(hit_rms/2.) >> hv12" , "hit_view == 1 && wvfm_bin_no == 12");
  t2->Draw("(hit_rms/2.) >> hv13" , "hit_view == 1 && wvfm_bin_no == 13");
  t2->Draw("(hit_rms/2.) >> hv14" , "hit_view == 1 && wvfm_bin_no == 14");
  t2->Draw("(hit_rms/2.) >> hv15" , "hit_view == 1 && wvfm_bin_no == 15");
  t2->Draw("(hit_rms/2.) >> hv16" , "hit_view == 1 && wvfm_bin_no == 16");
  t2->Draw("(hit_rms/2.) >> hv17" , "hit_view == 1 && wvfm_bin_no == 17");
  t2->Draw("(hit_rms/2.) >> hv18" , "hit_view == 1 && wvfm_bin_no == 18");
  t2->Draw("(hit_rms/2.) >> hv19" , "hit_view == 1 && wvfm_bin_no == 19");
  t2->Draw("(hit_rms/2.) >> hv20" , "hit_view == 1 && wvfm_bin_no == 20");
  t2->Draw("(hit_rms/2.) >> hv21" , "hit_view == 1 && wvfm_bin_no == 21");
  t2->Draw("(hit_rms/2.) >> hv22" , "hit_view == 1 && wvfm_bin_no == 22");
  t2->Draw("(hit_rms/2.) >> hv23" , "hit_view == 1 && wvfm_bin_no == 23");
  t2->Draw("(hit_rms/2.) >> hv24" , "hit_view == 1 && wvfm_bin_no == 24");

  t2->Draw("(hit_rms/2.) >> hy0"  , "hit_view == 2 && wvfm_bin_no == 0");
  t2->Draw("(hit_rms/2.) >> hy1"  , "hit_view == 2 && wvfm_bin_no == 1");
  t2->Draw("(hit_rms/2.) >> hy2"  , "hit_view == 2 && wvfm_bin_no == 2");
  t2->Draw("(hit_rms/2.) >> hy3"  , "hit_view == 2 && wvfm_bin_no == 3");
  t2->Draw("(hit_rms/2.) >> hy4"  , "hit_view == 2 && wvfm_bin_no == 4");
  t2->Draw("(hit_rms/2.) >> hy5"  , "hit_view == 2 && wvfm_bin_no == 5");
  t2->Draw("(hit_rms/2.) >> hy6"  , "hit_view == 2 && wvfm_bin_no == 6");
  t2->Draw("(hit_rms/2.) >> hy7"  , "hit_view == 2 && wvfm_bin_no == 7");
  t2->Draw("(hit_rms/2.) >> hy8"  , "hit_view == 2 && wvfm_bin_no == 8");
  t2->Draw("(hit_rms/2.) >> hy9"  , "hit_view == 2 && wvfm_bin_no == 9");
  t2->Draw("(hit_rms/2.) >> hy10" , "hit_view == 2 && wvfm_bin_no == 10");
  t2->Draw("(hit_rms/2.) >> hy11" , "hit_view == 2 && wvfm_bin_no == 11");
  t2->Draw("(hit_rms/2.) >> hy12" , "hit_view == 2 && wvfm_bin_no == 12");
  t2->Draw("(hit_rms/2.) >> hy13" , "hit_view == 2 && wvfm_bin_no == 13");
  t2->Draw("(hit_rms/2.) >> hy14" , "hit_view == 2 && wvfm_bin_no == 14");
  t2->Draw("(hit_rms/2.) >> hy15" , "hit_view == 2 && wvfm_bin_no == 15");
  t2->Draw("(hit_rms/2.) >> hy16" , "hit_view == 2 && wvfm_bin_no == 16");
  t2->Draw("(hit_rms/2.) >> hy17" , "hit_view == 2 && wvfm_bin_no == 17");
  t2->Draw("(hit_rms/2.) >> hy18" , "hit_view == 2 && wvfm_bin_no == 18");
  t2->Draw("(hit_rms/2.) >> hy19" , "hit_view == 2 && wvfm_bin_no == 19");
  t2->Draw("(hit_rms/2.) >> hy20" , "hit_view == 2 && wvfm_bin_no == 20");
  t2->Draw("(hit_rms/2.) >> hy21" , "hit_view == 2 && wvfm_bin_no == 21");
  t2->Draw("(hit_rms/2.) >> hy22" , "hit_view == 2 && wvfm_bin_no == 22");
  t2->Draw("(hit_rms/2.) >> hy23" , "hit_view == 2 && wvfm_bin_no == 23");
  t2->Draw("(hit_rms/2.) >> hy24" , "hit_view == 2 && wvfm_bin_no == 24");




  outfile->Write();

}
