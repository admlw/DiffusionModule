#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TF1.h"
#include "TPaveText.h"

void diffusionDemo()
{
  gStyle->SetOptStat(0);
  gStyle->SetTitleFont(132);
  gStyle->SetTextFont(22);
  gStyle->SetTextSize(0.04);
  //TCanvas *c1 = new TCanvas("c1","",550,1000);
  TCanvas *c1 = new TCanvas("c1","",650,1000);
  c1->SetLeftMargin(0.20);

  TH1D* h[24];
  TF1* lin[24];

  TFile *fin = new TFile("diffmod_CV.root", "READ");

  for (int i = 0; i < 25; i++) {

    //TString histoName = Form("calWireHistoBin_%i", i);
    TString histoName = Form("DiffusionModule/summed_waveform_bin_%i", i);
    if (fin->Get(histoName))
      h[i]=(TH1D*)fin->Get(histoName);
    else {
      std::cout << "No summed histogram in bin " << i << std::endl;
      continue;
    }

    int maxBin = h[i]->GetMaximumBin();
    int halfNBins = h[i]->GetNbinsX()/2.;

    if (maxBin >= halfNBins){
      for (int j = 0; j <= h[i]->GetNbinsX(); j++){
        double newVal = h[i]->GetBinContent(j);
        h[i]->SetBinContent(j + halfNBins - maxBin, newVal);
      }
    }

    if (maxBin < halfNBins){
      for (int j = h[i]->GetNbinsX(); j >= 0; j--){
        double newVal = h[i]->GetBinContent(j);
      
        h[i]->SetBinContent(j + halfNBins - maxBin, newVal);
      }
    }

    int nEntries = h[i]->GetEntries()/2.;
    double nStackedHistos = (double)nEntries/(((double)4600/(double)25)+2);
    double scaleFactor = 1.0 / (2.0*nStackedHistos);
    // Arbitrarily messed with this until the plot looked good
    h[i]->Scale(scaleFactor*20.);

    int maxTick = (h[i]->GetBinLowEdge(h[i]->GetMaximumBin())/2.); // Divide by two to convert to microseconds
    h[i]->GetXaxis()->SetLimits(-92, 92);
    h[i]->GetXaxis()->SetRangeUser(-20,20);

    TString lins = Form("lin%i", i);
    lin[i] = new TF1(lins, "[0]", -92, 92);
    lin[i]->SetParameter(0,maxTick);
    lin[i]->SetLineColor(kGray+2);
    lin[i]->SetLineWidth(1);
    lin[i]->SetLineStyle(2);

    for (int j = 1; j <= h[i]->GetNbinsX(); j++){
      h[i]->SetBinContent(j, h[i]->GetBinContent(j) + maxTick);
    }

    h[i]->SetLineColor(kBlue);
    h[i]->SetMinimum(400);
    h[i]->SetMaximum(3200);

    TAxis* xaxis = h[i]->GetXaxis();
    xaxis->SetTitle("Arb. time (ticks)");
    xaxis->SetTitleSize(0.05);

    TAxis* yaxis = h[i]->GetYaxis();
    yaxis->SetTitle("Drift time (#mus)");
    yaxis->SetTitleSize(0.05);
    yaxis->SetTitleOffset(1.75);

    h[i]->Draw("hist same");

    lin[i]->Draw("same");
    h[i]->Draw("same");

    //TPaveText* pt1 = new TPaveText(0.16, 0.86, 0.88, 0.92, "NDC");
    TPaveText* pt1 = new TPaveText(0.21, 0.86, 0.93, 0.91, "NDC");
    pt1->AddText("MicroBooNE Simulation");
    pt1->SetFillStyle(0);
    pt1->SetBorderSize(0);
    pt1->SetTextAlign(12);
    pt1->Draw("same");    

    //TPaveText* pt2 = new TPaveText(0.16, 0.82, 0.88, 0.86, "NDC");
    TPaveText* pt2 = new TPaveText(0.21, 0.82, 0.93, 0.87, "NDC");
    pt2->AddText("E = 273 V/cm");
    pt2->SetFillStyle(0);
    pt2->SetBorderSize(0);
    pt2->SetTextAlign(12);
    pt2->Draw("same");    


  }


  c1->SaveAs("diffusionDemo.png");
  c1->SaveAs("diffusionDemo.pdf");

}
