#include "StylePlots.h"

void NormalizeAndDraw(TH1D* h, std::string opt){
  if (h->Integral() > 0)
    h->Scale(1./h->Integral());
  h->GetYaxis()->SetRangeUser(0, 0.159);
  h->Draw(opt.c_str());
}

void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
   if (!C){
     std::cout << "canvas doesn't exist" << std::endl;
     return;
   }

   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;

   Float_t hSpacing = 0.0;
   Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;

   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;

   for (Int_t i=0;i<Nx;i++) {

      if (i==0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr-hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx-1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr-hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }

      for (Int_t j=0;j<Ny;j++) {

         if (j==0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu-vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         } else if (j == Ny-1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = tMargin / (vposu-vposd);
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }

         C->cd(0);

         char name[16];
         sprintf(name,"pad_%i_%i",i,j);
         TPad *pad = (TPad*) gROOT->FindObject(name);
         if (pad) delete pad;
         pad = new TPad(name,"",hposl,vposd,hposr,vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);

         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);

         std::cout << "Created pad " << pad->GetName() << std::endl;

         pad->Draw();
      }
   }
}

void make_diffusion_t0tagging_plots(){

  TFile* f2 = new TFile("/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionAna/diffmod_run3_crt_Aug2020_newFV_bugFix.root", "read");
  TTree* t2 = (TTree*)f2->Get("DiffusionModule/difftree");

  TH1D* hu0  = new TH1D("hu0" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu12 = new TH1D("hu12", ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hu24 = new TH1D("hu24", ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv0  = new TH1D("hv0" , ";Hit RMS (#mus);Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv12 = new TH1D("hv12", ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hv24 = new TH1D("hv24", ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy0  = new TH1D("hy0" , ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy12 = new TH1D("hy12", ";;Number of Hits (area norm.)", 50, 0.81, 2.59);
  TH1D* hy24 = new TH1D("hy24", ";;Number of Hits (area norm.)", 50, 0.81, 2.59);

  t2->Draw("(hit_rms/2.) >> hu0"  , "hit_view == 0 && wvfm_bin_no == 0");
  t2->Draw("(hit_rms/2.) >> hu12" , "hit_view == 0 && wvfm_bin_no == 12");
  t2->Draw("(hit_rms/2.) >> hu24" , "hit_view == 0 && wvfm_bin_no == 24");
  t2->Draw("(hit_rms/2.) >> hv0"  , "hit_view == 1 && wvfm_bin_no == 0");
  t2->Draw("(hit_rms/2.) >> hv12" , "hit_view == 1 && wvfm_bin_no == 12");
  t2->Draw("(hit_rms/2.) >> hv24" , "hit_view == 1 && wvfm_bin_no == 24");
  t2->Draw("(hit_rms/2.) >> hy0"  , "hit_view == 2 && wvfm_bin_no == 0");
  t2->Draw("(hit_rms/2.) >> hy12" , "hit_view == 2 && wvfm_bin_no == 12");
  t2->Draw("(hit_rms/2.) >> hy24" , "hit_view == 2 && wvfm_bin_no == 24");

  SetGenericStyle();
  gROOT->ForceStyle();

  TCanvas* ct = new TCanvas("ct", "ct", 1200, 600);
  ct->UseCurrentStyle();
  gStyle->SetLabelFont(43, "xyz");
  gStyle->SetLabelSize(10);

  CanvasPartition(ct, 3, 1, 0.1, 0.05, 0.1, 0.1);


  TPad* pad_0_0 = ((TPad*)gROOT->FindObject("pad_0_0"));
  pad_0_0->SetRightMargin(0.05);
  pad_0_0->cd();
  hu0->SetLabelFont(43, "xyz");
  hu0->SetLabelSize(20, "xyz");
  hu0->SetTitleFont(43, "xyz");
  hu0->SetTitleSize(25, "xyz");
  hu0->GetYaxis()->CenterTitle();
  hu0->GetYaxis()->SetTitleOffset(2.4);
  hu0->GetYaxis()->SetRangeUser(0, 0.1);
  hu0 ->SetLineWidth(2); hu0 ->SetLineColor(kPTRed);
  hu12->SetLineWidth(2); hu12->SetLineColor(kPTOrange);
  hu24->SetLineWidth(2); hu24->SetLineColor(kPTLightBlue);
  NormalizeAndDraw(hu0, "hist");
  NormalizeAndDraw(hu12, "hist same");
  NormalizeAndDraw(hu24, "hist same");

  TLegend* leg = new TLegend(0.3, 0.73, 0.95, 0.97);
  leg->AddEntry(hu0 , "Drift time = 45 #mus", "l");
  leg->AddEntry(hu12, "Drift time = 1150 #mus", "l");
  leg->AddEntry(hu24, "Drift time = 2254 #mus", "l");
  leg->SetLineWidth(0);
  leg->SetFillStyle(0);
  leg->Draw();


  TPad* pad_1_0 = ((TPad*)gROOT->FindObject("pad_1_0"));
  pad_1_0->SetRightMargin(0.05);
  pad_1_0->cd();
  hv0->SetLabelFont(43, "xyz");
  hv0->SetLabelSize(20, "xyz");
  hv0->GetXaxis()->CenterTitle();
  hv0->SetTitleFont(43, "xyz");
  hv0->SetTitleSize(25, "xyz");
  hv0->GetYaxis()->SetRangeUser(0, 0.1);
  hv0 ->SetLineWidth(2); hv0 ->SetLineColor(kPTRed);
  hv12->SetLineWidth(2); hv12->SetLineColor(kPTOrange);
  hv24->SetLineWidth(2); hv24->SetLineColor(kPTLightBlue);
  NormalizeAndDraw(hv0, "hist");
  NormalizeAndDraw(hv12, "hist same");
  NormalizeAndDraw(hv24, "hist same");
  TPad* pad_2_0 = ((TPad*)gROOT->FindObject("pad_2_0"));
  pad_2_0->cd();
  hy0->SetLabelFont(43, "xyz");
  hy0->SetLabelSize(20, "xyz");
  hy0->GetYaxis()->SetRangeUser(0, 0.1);
  hy0 ->SetLineWidth(2); hy0 ->SetLineColor(kPTRed);
  hy12->SetLineWidth(2); hy12->SetLineColor(kPTOrange);
  hy24->SetLineWidth(2); hy24->SetLineColor(kPTLightBlue);
  NormalizeAndDraw(hy0, "hist");
  NormalizeAndDraw(hy12, "hist same");
  NormalizeAndDraw(hy24, "hist same");
  ApplyLabel(kData, 0.67, 0.95);
  
  ct->cd();
  TLatex *upl = new TLatex(0.185, 0.91, "U Plane");
  upl->SetNDC();
  upl->SetTextSize(0.05);
  upl->Draw();
  TLatex *vpl = new TLatex(0.47, 0.91, "V Plane");
  vpl->SetNDC();
  vpl->SetTextSize(0.05);
  vpl->Draw();
  TLatex *ypl = new TLatex(0.76, 0.91, "Y Plane");
  ypl->SetNDC();
  ypl->SetTextSize(0.05);
  ypl->Draw();

  ct->SaveAs("t0tagging.pdf");

}
