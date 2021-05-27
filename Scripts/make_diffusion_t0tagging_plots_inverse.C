#include "StylePlots.h"

void plotText(std::string text, float posx, float posy){ 
  TLatex* textlat = new TLatex(posx, posy, text.c_str());
  textlat->SetNDC();
  textlat->SetTextFont(43);
  textlat->SetTextSize(20);
  textlat->SetTextAlign(22);
  textlat->Draw();
}


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
    float counter       = h_alls[i]->GetBinContent(maxValBin)/2.;
    float integralLeft  = h_alls[i]->Integral(1, maxValBin-1)+counter;
    float integralRight = h_alls[i]->Integral(maxValBin+1, h_alls[i]->GetNbinsX())+counter;

    for (int j = maxValBin-1; j > 0; j--){
      counter += h_alls[i]->GetBinContent(j);
      if (counter/integralLeft >= percentile){
        lowBin=j;
        break;
      }
    }

    counter = h_alls[i]->GetBinContent(maxValBin)/2.;
    for (int j = maxValBin+1; j < h_alls[i]->GetNbinsX()+1; j++){
      counter += h_alls[i]->GetBinContent(j);
      std::cout << "counter: " << counter 
                << " integralRight: " << integralRight
                << " c/i" << counter/integralRight << std::endl;
      if (counter/integralRight >= percentile){
        highBin=j;
        break;
      }
    }

    std::cout << "lowBin: " << lowBin << " highBin: " << highBin << std::endl;

    TString name = h_alls[i]->GetName();
    name.Append(std::to_string(percentile)+"perc");

    h_outs.push_back(new TH1D(name, "", h_alls[i]->GetNbinsX(), 0, 2350));
    
    for (int j = lowBin; j < highBin; ++j){
      h_outs[i]->SetBinContent(j, h_alls[i]->GetBinContent(j));
      h_outs[i]->SetLineWidth(0);
    }

  }

}

void make_diffusion_t0tagging_plots_inverse(){

  TFile* f = new TFile("/pnfs/uboone/persistent/diffusion_analysis_final_files/diffmod_sim_DL_3.74_paper.root", "read");
  TTree* t = (TTree*)f->Get("DiffusionModule/difftree");
  TFile* f2 = new TFile("/pnfs/uboone/persistent/diffusion_analysis_final_files/diffmod_run3data_paper.root", "read");
  TTree* t2 = (TTree*)f2->Get("DiffusionModule/difftree");

  SetGenericStyle();
  gROOT->ForceStyle();

  TH2D* hu_dat = new TH2D("hu_dat", ";;", 50, 0, 2300, 10, 1.1, 2.1);
  TH2D* hv_dat = new TH2D("hv_dat", ";;", 50, 0, 2300, 10, 1.1, 2.1);
  TH2D* hy_dat = new TH2D("hy_dat", ";;", 50, 0, 2300, 10, 1.1, 2.1);
  TH2D* hy_sim = new TH2D("hy_sim", ";;", 50, 0, 2300, 10, 1.1, 2.1);
  
  t->Draw("hit_rms/2:(hit_peak_time_t0corr-800)/2. >> hu_dat", "hit_view == 0");
  t->Draw("hit_rms/2:(hit_peak_time_t0corr-800)/2. >> hv_dat", "hit_view == 1");
  t->Draw("hit_rms/2:(hit_peak_time_t0corr-800)/2. >> hy_dat", "hit_view == 2");
  t2->Draw("hit_rms/2:(hit_peak_time_t0corr-800)/2. >> hy_sim", "hit_view == 2");


  areaNormaliseX(hu_dat);
  areaNormaliseX(hv_dat);
  areaNormaliseX(hy_dat);
  areaNormaliseX(hy_sim);
  std::vector<TH1D*> hu_dat_alls;
  std::vector<TH1D*> hv_dat_alls;
  std::vector<TH1D*> hy_dat_alls;
  std::vector<TH1D*> hy_sim_alls;

  //std::vector<TH1D*> hy_dat_1sig;
  //std::vector<TH1D*> hy_dat_2sig;
  //std::vector<TH1D*> hy_dat_3sig;

  for (int i = 0 ; i < 9; ++i){

    float init      = 1.1;
    float step_size = 0.1;
    float low_cut   = init + (i*step_size);
    float high_cut  = low_cut+step_size;
     
    std::string nameu  = "hu_dat" + std::to_string(i);
    std::string namev  = "hv_dat" + std::to_string(i);
    std::string namey  = "hy_dat" + std::to_string(i);
    std::string namesimy  = "hy_sim" + std::to_string(i);

    std::string title = ";Hit Time (#mus); Number of Hits";

    hu_dat_alls.push_back((TH1D*)hu_dat->ProjectionX(nameu.c_str(),i+1,i+2));
    hv_dat_alls.push_back((TH1D*)hv_dat->ProjectionX(namev.c_str(),i+1,i+2));
    hy_dat_alls.push_back((TH1D*)hy_dat->ProjectionX(namey.c_str(),i+1,i+2));
    hy_sim_alls.push_back((TH1D*)hy_sim->ProjectionX(namesimy.c_str(),i+1,i+2));
    hu_dat_alls[i]->SetTitle(title.c_str());
    hv_dat_alls[i]->SetTitle(title.c_str());
    hy_dat_alls[i]->SetTitle(title.c_str());
    hy_sim_alls[i]->SetTitle(title.c_str());

    hy_dat_alls[i]->GetXaxis()->SetTitleFont(43);
    hy_dat_alls[i]->GetXaxis()->SetLabelFont(43);
    hy_dat_alls[i]->GetXaxis()->SetTitleSize(24);
    hy_dat_alls[i]->GetXaxis()->SetLabelSize(24);
    hy_dat_alls[i]->GetXaxis()->SetTitleOffset(2.5);
    hy_dat_alls[i]->GetXaxis()->CenterTitle();
    hy_dat_alls[i]->GetYaxis()->SetTitleFont(43);
    hy_dat_alls[i]->GetYaxis()->SetLabelFont(43);
    hy_dat_alls[i]->GetYaxis()->SetTitleSize(24);
    hy_dat_alls[i]->GetYaxis()->SetLabelSize(24);
    hy_dat_alls[i]->GetYaxis()->SetTitleOffset(3.2);
    hy_dat_alls[i]->GetYaxis()->CenterTitle();
    hy_dat_alls[i]->GetXaxis()->SetLimits(1, 2300);

    hy_sim_alls[i]->GetXaxis()->SetTitleFont(43);
    hy_sim_alls[i]->GetXaxis()->SetLabelFont(43);
    hy_sim_alls[i]->GetXaxis()->SetTitleSize(24);
    hy_sim_alls[i]->GetXaxis()->SetLabelSize(24);
    hy_sim_alls[i]->GetXaxis()->SetTitleOffset(2.5);
    hy_sim_alls[i]->GetXaxis()->CenterTitle();
    hy_sim_alls[i]->GetYaxis()->SetTitleFont(43);
    hy_sim_alls[i]->GetYaxis()->SetLabelFont(43);
    hy_sim_alls[i]->GetYaxis()->SetTitleSize(24);
    hy_sim_alls[i]->GetYaxis()->SetLabelSize(24);
    hy_sim_alls[i]->GetYaxis()->SetTitleOffset(3.2);
    hy_sim_alls[i]->GetYaxis()->CenterTitle();
    hy_sim_alls[i]->GetXaxis()->SetLimits(1, 2300);

  }
  for (int i = 0; i < hy_dat_alls.size(); ++i){
    hy_dat_alls[i]->Sumw2();
    hy_dat_alls[i]->Scale(1./hy_dat_alls[0]->Integral());
    hy_dat_alls[i]->GetYaxis()->SetRangeUser(0.001, 0.999);
    hy_dat_alls[i]->SetLineColor(kBlack);
    hy_dat_alls[i]->SetMarkerStyle(20);
    hy_dat_alls[i]->SetMarkerSize(0.5);
    hy_dat_alls[i]->SetMarkerColor(kBlack);

    hu_dat_alls[i]->Sumw2();
    hu_dat_alls[i]->Scale(1./hu_dat_alls[0]->Integral());
    hu_dat_alls[i]->GetYaxis()->SetRangeUser(0.001, 0.999);
    hu_dat_alls[i]->SetLineColor(kAzure+1);

    hv_dat_alls[i]->Sumw2();
    hv_dat_alls[i]->Scale(1./hv_dat_alls[0]->Integral());
    hv_dat_alls[i]->GetYaxis()->SetRangeUser(0.001, 0.999);
    hv_dat_alls[i]->SetLineColor(kGreen+1);

    hy_sim_alls[i]->Sumw2();
    hy_sim_alls[i]->Scale(1./hy_sim_alls[0]->Integral());
    hy_sim_alls[i]->Scale(hy_dat_alls[i]->Integral()/hy_sim_alls[i]->Integral());
    hy_sim_alls[i]->GetYaxis()->SetRangeUser(0.001, 1.79);
    hy_sim_alls[i]->SetLineColor(kGreen+1);

  }

  //GetBoundaries(hy_dat_alls, hy_dat_1sig, 0.682); 
  //GetBoundaries(hy_dat_alls, hy_dat_2sig, 0.954); 
  //GetBoundaries(hy_dat_alls, hy_dat_3sig, 0.996); 

  TCanvas* tp = new TCanvas("tp", "tp", 1000, 1000);
  CanvasPartition(tp, 3, 3, 0.1, 0.1, 0.15, 0.1);

  // pad_0_2
  TPad* pad_0_2 = ((TPad*)gROOT->FindObject("pad_0_2"));
  pad_0_2->SetRightMargin(0.0);
  pad_0_2->cd();
  hy_sim_alls[0]->Draw("hist");
  hy_dat_alls[0]->Draw("p");
  plotText("1.1 #mus < Hit RMS < 1.2 #mus", 0.64, 0.65);
  plotText("MicroBooNE Data", 0.64, 0.56);
  //hu_dat_alls[0]->Draw("hist same");
  //hv_dat_alls[0]->Draw("hist same");
  //hy_dat_3sig[0]->SetFillColor(kAzure-9);
  //hy_dat_3sig[0]->Draw("same hist");
  //hy_dat_2sig[0]->SetFillColor(kAzure-8);
  //hy_dat_2sig[0]->Draw("same hist");
  //hy_dat_1sig[0]->SetFillColor(kAzure-7);
  //hy_dat_1sig[0]->Draw("same hist");
  //hy_dat_alls[0]->Draw("same hist");

  //TLegend* leg = new TLegend(0.15, 0.6, 0.85, 0.9);
  //leg->AddEntry(hy_dat_1sig[0], "1 Sigma", "f");
  //leg->AddEntry(hy_dat_2sig[0], "2 Sigma", "f");
  //leg->AddEntry(hy_dat_3sig[0], "3 Sigma", "f");
  //leg->SetLineWidth(0);
  //leg->SetFillStyle(0);
  //leg->Draw("same");

  //TLegend* leg = new TLegend(0.3, 0.4, 0.95, 0.6);
  //leg->AddEntry(hy_sim_alls[0], "Sim. (D_{L} = 3.23 cm^{2}/s)", "l");
  //leg->AddEntry(hy_dat_alls[0], "MicroBooNE Data", "p");
  //leg->SetLineWidth(0);
  //leg->SetFillStyle(0);
  //leg->Draw("same");


  //pad_1_2
  TPad* pad_1_2 = ((TPad*)gROOT->FindObject("pad_1_2"));
  pad_1_2->SetRightMargin(0.0);
  pad_1_2->cd();
  hy_sim_alls[1]->Draw("hist");
  hy_dat_alls[1]->Draw("p");
  plotText("1.2 #mus < Hit RMS < 1.3 #mus", 0.5, 0.65);
  //hu_dat_alls[1]->Draw("hist same");
  //hv_dat_alls[1]->Draw("hist same");

  //hy_dat_3sig[1]->SetFillColor(kAzure-9);
  //hy_dat_3sig[1]->Draw("same hist");
  //hy_dat_2sig[1]->SetFillColor(kAzure-8);
  //hy_dat_2sig[1]->Draw("same hist");
  //hy_dat_1sig[1]->SetFillColor(kAzure-7);
  //hy_dat_1sig[1]->Draw("same hist");
  //hy_dat_alls[1]->Draw("same hist");

  //pad_2_2
  TPad* pad_2_2 = ((TPad*)gROOT->FindObject("pad_2_2"));
  //pad_2_2->SetRightMargin(0.0);
  pad_2_2->cd();
  hy_sim_alls[2]->Draw("hist");
  hy_dat_alls[2]->Draw("p");
  plotText("1.3 #mus < Hit RMS < 1.4 #mus", 0.38, 0.65);
  //hu_dat_alls[2]->Draw("hist same");
  //hv_dat_alls[2]->Draw("hist same");

  //hy_dat_3sig[2]->SetFillColor(kAzure-9);
  //hy_dat_3sig[2]->Draw("same hist");
  //hy_dat_2sig[2]->SetFillColor(kAzure-8);
  //hy_dat_2sig[2]->Draw("same hist");
  //hy_dat_1sig[2]->SetFillColor(kAzure-7);
  //hy_dat_1sig[2]->Draw("same hist");
  //hy_dat_alls[2]->Draw("same hist");

  // pad_0_1
  TPad* pad_0_1 = ((TPad*)gROOT->FindObject("pad_0_1"));
  pad_0_1->SetRightMargin(0.0);
  pad_0_1->cd();
  hy_sim_alls[3]->Draw("hist");
  hy_dat_alls[3]->Draw("p");
  plotText("1.4 #mus < Hit RMS < 1.5 #mus", 0.64, 0.89);
  //hu_dat_alls[3]->Draw("hist same");
  //hv_dat_alls[3]->Draw("hist same");

  //hy_dat_3sig[3]->SetFillColor(kAzure-9);
  //hy_dat_3sig[3]->Draw("same hist");
  //hy_dat_2sig[3]->SetFillColor(kAzure-8);
  //hy_dat_2sig[3]->Draw("same hist");
  //hy_dat_1sig[3]->SetFillColor(kAzure-7);
  //hy_dat_1sig[3]->Draw("same hist");
  //hy_dat_alls[3]->Draw("same hist");

  //pad_1_1
  TPad* pad_1_1 = ((TPad*)gROOT->FindObject("pad_1_1"));
  pad_1_1->SetRightMargin(0.0);
  pad_1_1->cd();
  hy_sim_alls[4]->Draw("hist");
  hy_dat_alls[4]->Draw("p");
  plotText("1.5 #mus < Hit RMS < 1.6 #mus", 0.5, 0.89);
  //hu_dat_alls[4]->Draw("hist same");
  //hv_dat_alls[4]->Draw("hist same");

  //hy_dat_3sig[4]->SetFillColor(kAzure-9);
  //hy_dat_3sig[4]->Draw("same hist");
  //hy_dat_2sig[4]->SetFillColor(kAzure-8);
  //hy_dat_2sig[4]->Draw("same hist");
  //hy_dat_1sig[4]->SetFillColor(kAzure-7);
  //hy_dat_1sig[4]->Draw("same hist");
  //hy_dat_alls[4]->Draw("same hist");

  //pad_2_1
  TPad* pad_2_1 = ((TPad*)gROOT->FindObject("pad_2_1"));
  //pad_2_1->SetRightMargin(0.0);
  pad_2_1->cd();
  hy_sim_alls[5]->Draw("hist");
  hy_dat_alls[5]->Draw("p");
  plotText("1.6 #mus < Hit RMS < 1.7 #mus", 0.38, 0.89);
  //hu_dat_alls[5]->Draw("hist same");
  //hv_dat_alls[5]->Draw("hist same");

  //hy_dat_3sig[5]->SetFillColor(kAzure-9);
  //hy_dat_3sig[5]->Draw("same hist");
  //hy_dat_2sig[5]->SetFillColor(kAzure-8);
  //hy_dat_2sig[5]->Draw("same hist");
  //hy_dat_1sig[5]->SetFillColor(kAzure-7);
  //hy_dat_1sig[5]->Draw("same hist");
  //hy_dat_alls[5]->Draw("same hist");

  // pad_0_0
  TPad* pad_0_0 = ((TPad*)gROOT->FindObject("pad_0_0"));
  pad_0_0->SetRightMargin(0.0);
  pad_0_0->cd();
  hy_sim_alls[6]->Draw("hist");
  hy_dat_alls[6]->Draw("p");
  plotText("1.7 #mus < Hit RMS < 1.8 #mus", 0.64, 0.94);
  //hu_dat_alls[6]->Draw("hist same");
  //hv_dat_alls[6]->Draw("hist same");

  //hy_dat_3sig[6]->SetFillColor(kAzure-9);
  //hy_dat_3sig[6]->Draw("same hist");
  //hy_dat_2sig[6]->SetFillColor(kAzure-8);
  //hy_dat_2sig[6]->Draw("same hist");
  //hy_dat_1sig[6]->SetFillColor(kAzure-7);
  //hy_dat_1sig[6]->Draw("same hist");
  //hy_dat_alls[6]->Draw("same hist");

  //pad_1_0
  TPad* pad_1_0 = ((TPad*)gROOT->FindObject("pad_1_0"));
  pad_1_0->SetRightMargin(0.0);
  pad_1_0->cd();
  hy_sim_alls[7]->Draw("hist");
  hy_dat_alls[7]->Draw("p");
  plotText("1.8 #mus < Hit RMS < 1.9 #mus", 0.5, 0.94);
  //hu_dat_alls[7]->Draw("hist same");
  //hv_dat_alls[7]->Draw("hist same");

  //hy_dat_3sig[7]->SetFillColor(kAzure-9);
  //hy_dat_3sig[7]->Draw("same hist");
  //hy_dat_2sig[7]->SetFillColor(kAzure-8);
  //hy_dat_2sig[7]->Draw("same hist");
  //hy_dat_1sig[7]->SetFillColor(kAzure-7);
  //hy_dat_1sig[7]->Draw("same hist");
  //hy_dat_alls[7]->Draw("same hist");

  //pad_2_0
  TPad* pad_2_0 = ((TPad*)gROOT->FindObject("pad_2_0"));
  //pad_2_0->SetRightMargin(0.0);
  pad_2_0->cd();
  hy_sim_alls[8]->Draw("hist");
  hy_dat_alls[8]->Draw("p");
  plotText("1.9 #mus < Hit RMS < 2.0 #mus", 0.38, 0.94);
  //hu_dat_alls[8]->Draw("hist same");
  //hv_dat_alls[8]->Draw("hist same");

  //hy_dat_3sig[8]->SetFillColor(kAzure-9);
  //hy_dat_3sig[8]->Draw("same hist");
  //hy_dat_2sig[8]->SetFillColor(kAzure-8);
  //hy_dat_2sig[8]->Draw("same hist");
  //hy_dat_1sig[8]->SetFillColor(kAzure-7);
  //hy_dat_1sig[8]->Draw("same hist");
  //hy_dat_alls[8]->Draw("same hist");

  tp->SaveAs("hit_time_per_hit_rms.pdf");
}
