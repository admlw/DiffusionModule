#include "StylePlots.h"

void areaNormaliseX(TH2* h){

  for (int i = 0; i < h->GetNbinsX()+1; ++i){

    float content = 0.0;
    for (int j = 0; j < h->GetNbinsY()+1; ++j){
      content += h->GetBinContent(i,j);
    }

    for (int j = 0; j < h->GetNbinsY()+1; ++j){
      h->SetBinContent(i,j,h->GetBinContent(i,j)/content);
    }

  }

}


// partition canvas
void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
   if (!C) return;

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

         pad->Draw();
      }
   }
}


void GetBoundaries(std::vector<TH1D*> &h_alls, std::vector<TH1D*> &h_outs, float percentile){
  
  for (int i = 0; i < h_alls.size(); ++i){
    std::cout << "i: " << i << std::endl;

    float maxValBin = h_alls[i]->GetMaximumBin();
    float integral  = h_alls[i]->Integral();

    int lowBin  = 0;
    int highBin = h_alls[i]->GetNbinsX()+1;
    // count half of max Val bin as being on the left, and half on the right
    float counter = h_alls[i]->GetBinContent(maxValBin)/2.;
    std::cout << "counter: " << counter << std::endl;
    for (int j = maxValBin-1; j > 0; j--){
      counter += h_alls[i]->GetBinContent(j);
      std::cout << "counter: " << counter << " integral: " << integral << std::endl << " counter/integral: " << counter/integral << std::endl;
      if (counter/integral >= percentile/2.){
        lowBin=j;
        break;
      }
    }
    counter = h_alls[i]->GetBinContent(maxValBin)/2.;
    std::cout << "counter: " << counter << std::endl;
    for (int j = maxValBin+1; j < h_alls[i]->GetNbinsX()+1; j++){
      counter += h_alls[i]->GetBinContent(j);
      std::cout << "counter: " << counter << " integral: " << integral << std::endl << " counter/integral: " << counter/integral << std::endl;
      if (counter/integral >= percentile/2.){
        highBin=j+1;
        break;
      }
    }

    std::cout << "lowBin: " << lowBin << " highBin: " << highBin << std::endl;

    TString name = h_alls[i]->GetName();
    name.Append(std::to_string(percentile)+"perc");

    h_outs.push_back(new TH1D(name, "", h_alls[i]->GetNbinsX(), 0, 2300));
    
    for (int j = lowBin; j < highBin; ++j){
      h_outs[i]->SetBinContent(j, h_alls[i]->GetBinContent(j));
      h_outs[i]->SetLineWidth(0);
    }

  }

}

void make_diffusion_t0tagging_plots_inverse(){

  TFile* f2 = new TFile("/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionAna/diffmod_run3_crt_Aug2020_newFV_bugFix.root", "read");
  TTree* t2 = (TTree*)f2->Get("DiffusionModule/difftree");

  SetGenericStyle();
  gROOT->ForceStyle();

  TH2D* h = new TH2D("h", ";;", 50, 0., 2300, 10, 1.1, 2.1);
  t2->Draw("hit_rms/2:(hit_peak_time_t0corr-800)/2. >> h", "hit_view == 2");

  areaNormaliseX(h);
  h->Draw("colz");
  std::vector<TH1D*> hy_alls;
  std::vector<TH1D*> hy_1sig;
  std::vector<TH1D*> hy_2sig;
  std::vector<TH1D*> hy_3sig;

  for (int i = 0 ; i < 9; ++i){

    float init      = 1.1;
    float step_size = 0.1;
    float low_cut   = init + (i*step_size);
    float high_cut  = low_cut+step_size;

    std::string name  = "hy" + std::to_string(i);
    std::string title = std::to_string(low_cut) + " < Hit RMS < " + std::to_string(high_cut) + ";Hit Time (#mus); Number of Hits";

    hy_alls.push_back((TH1D*)h->ProjectionX(name.c_str(),i+1,i+2));
    hy_alls[i]->SetTitle(title.c_str());
  }

  for (int i = 0; i < hy_alls.size(); ++i){
    hy_alls[i]->Sumw2();
    hy_alls[i]->Scale(1./hy_alls[0]->Integral());
    hy_alls[i]->GetYaxis()->SetRangeUser(0, 1.0);
    hy_alls[i]->SetLineColor(kBlack);
  }

  GetBoundaries(hy_alls, hy_1sig, 0.682); 
  GetBoundaries(hy_alls, hy_2sig, 0.954); 
  GetBoundaries(hy_alls, hy_3sig, 0.996); 

  TCanvas* tp = new TCanvas("tp", "tp", 1000, 1000);
  CanvasPartition(tp, 3, 3, 0.1, 0.1, 0.15, 0.1);

  // pad_0_2
  TPad* pad_0_2 = ((TPad*)gROOT->FindObject("pad_0_2"));
  pad_0_2->SetRightMargin(0.04);
  pad_0_2->cd();
  hy_alls[0]->Draw("hist");
  hy_3sig[0]->SetFillColor(kAzure-9);
  hy_3sig[0]->Draw("same hist");
  hy_2sig[0]->SetFillColor(kAzure-8);
  hy_2sig[0]->Draw("same hist");
  hy_1sig[0]->SetFillColor(kAzure-7);
  hy_1sig[0]->Draw("same hist");
  hy_alls[0]->Draw("same hist");

  TLegend* leg = new TLegend(0.15, 0.6, 0.85, 0.9);
  leg->AddEntry(hy_1sig[0], "1 Sigma", "f");
  leg->AddEntry(hy_2sig[0], "2 Sigma", "f");
  leg->AddEntry(hy_3sig[0], "3 Sigma", "f");
  leg->SetLineWidth(0);
  leg->SetFillStyle(0);
  leg->Draw("same");

  //pad_1_2
  TPad* pad_1_2 = ((TPad*)gROOT->FindObject("pad_1_2"));
  pad_1_2->SetRightMargin(0.04);
  pad_1_2->cd();
  hy_alls[1]->Draw("hist");
  hy_3sig[1]->SetFillColor(kAzure-9);
  hy_3sig[1]->Draw("same hist");
  hy_2sig[1]->SetFillColor(kAzure-8);
  hy_2sig[1]->Draw("same hist");
  hy_1sig[1]->SetFillColor(kAzure-7);
  hy_1sig[1]->Draw("same hist");
  hy_alls[1]->Draw("same hist");

  //pad_2_2
  TPad* pad_2_2 = ((TPad*)gROOT->FindObject("pad_2_2"));
  //pad_2_2->SetRightMargin(0.04);
  pad_2_2->cd();
  hy_alls[2]->Draw("hist");
  hy_3sig[2]->SetFillColor(kAzure-9);
  hy_3sig[2]->Draw("same hist");
  hy_2sig[2]->SetFillColor(kAzure-8);
  hy_2sig[2]->Draw("same hist");
  hy_1sig[2]->SetFillColor(kAzure-7);
  hy_1sig[2]->Draw("same hist");
  hy_alls[2]->Draw("same hist");

  // pad_0_1
  TPad* pad_0_1 = ((TPad*)gROOT->FindObject("pad_0_1"));
  pad_0_1->SetRightMargin(0.04);
  pad_0_1->cd();
  hy_alls[3]->Draw("hist");
  hy_3sig[3]->SetFillColor(kAzure-9);
  hy_3sig[3]->Draw("same hist");
  hy_2sig[3]->SetFillColor(kAzure-8);
  hy_2sig[3]->Draw("same hist");
  hy_1sig[3]->SetFillColor(kAzure-7);
  hy_1sig[3]->Draw("same hist");
  hy_alls[3]->Draw("same hist");

  //pad_1_1
  TPad* pad_1_1 = ((TPad*)gROOT->FindObject("pad_1_1"));
  pad_1_1->SetRightMargin(0.04);
  pad_1_1->cd();
  hy_alls[4]->Draw("hist");
  hy_3sig[4]->SetFillColor(kAzure-9);
  hy_3sig[4]->Draw("same hist");
  hy_2sig[4]->SetFillColor(kAzure-8);
  hy_2sig[4]->Draw("same hist");
  hy_1sig[4]->SetFillColor(kAzure-7);
  hy_1sig[4]->Draw("same hist");
  hy_alls[4]->Draw("same hist");

  //pad_2_1
  TPad* pad_2_1 = ((TPad*)gROOT->FindObject("pad_2_1"));
  //pad_2_1->SetRightMargin(0.04);
  pad_2_1->cd();
  hy_alls[5]->Draw("hist");
  hy_3sig[5]->SetFillColor(kAzure-9);
  hy_3sig[5]->Draw("same hist");
  hy_2sig[5]->SetFillColor(kAzure-8);
  hy_2sig[5]->Draw("same hist");
  hy_1sig[5]->SetFillColor(kAzure-7);
  hy_1sig[5]->Draw("same hist");
  hy_alls[5]->Draw("same hist");

  // pad_0_0
  TPad* pad_0_0 = ((TPad*)gROOT->FindObject("pad_0_0"));
  pad_0_0->SetRightMargin(0.04);
  pad_0_0->cd();
  hy_alls[6]->Draw("hist");
  hy_3sig[6]->SetFillColor(kAzure-9);
  hy_3sig[6]->Draw("same hist");
  hy_2sig[6]->SetFillColor(kAzure-8);
  hy_2sig[6]->Draw("same hist");
  hy_1sig[6]->SetFillColor(kAzure-7);
  hy_1sig[6]->Draw("same hist");
  hy_alls[6]->Draw("same hist");

  //pad_1_0
  TPad* pad_1_0 = ((TPad*)gROOT->FindObject("pad_1_0"));
  pad_1_0->SetRightMargin(0.04);
  pad_1_0->cd();
  hy_alls[7]->Draw("hist");
  hy_3sig[7]->SetFillColor(kAzure-9);
  hy_3sig[7]->Draw("same hist");
  hy_2sig[7]->SetFillColor(kAzure-8);
  hy_2sig[7]->Draw("same hist");
  hy_1sig[7]->SetFillColor(kAzure-7);
  hy_1sig[7]->Draw("same hist");
  hy_alls[7]->Draw("same hist");

  //pad_2_0
  TPad* pad_2_0 = ((TPad*)gROOT->FindObject("pad_2_0"));
  //pad_2_0->SetRightMargin(0.04);
  pad_2_0->cd();
  hy_alls[8]->Draw("hist");
  hy_3sig[8]->SetFillColor(kAzure-9);
  hy_3sig[8]->Draw("same hist");
  hy_2sig[8]->SetFillColor(kAzure-8);
  hy_2sig[8]->Draw("same hist");
  hy_1sig[8]->SetFillColor(kAzure-7);
  hy_1sig[8]->Draw("same hist");
  hy_alls[8]->Draw("same hist");

  tp->SaveAs("tmp.pdf");

  TCanvas* c2 = new TCanvas("c2", "c2", 500, 500);
  TH1D* h_test = (TH1D*)hy_alls[4]->Clone("h_test");
  h_test->Draw();
  TF1* gausfit = new TF1("gausfit", "gaus");
  //f->SetParameter(0,1);
  //f->SetParameter(1, 1000);
  //f->SetParameter(2,500);
  h_test->Fit(gausfit, "q", "", 500, 1500);
}
