#include "../StylePlots.h"

double getRms2(TH1D* h){

    double mean = 0;
    double mean2 = 0;
    double totalCharge = 0;
    double rms2 = 0;
    double threshold = 10.0 * h->GetMaximum()/100.0;

    for (int j = 1; j < h->GetNbinsX(); j++) {

        h->SetBinContent(j, h->GetBinContent(j) - threshold);

        if (h->GetBinContent(j) >= 0){

            totalCharge = totalCharge + h->GetBinContent(j);
            mean = mean + (double)h->GetBinCenter(j) * 0.5 * (double)h->GetBinContent(j);
            mean2 = mean2 + (double)h->GetBinCenter(j) * 0.5 * (double)h->GetBinCenter(j) *0.5 * (double)h->GetBinContent(j);
        }

        h->SetBinContent(j, h->GetBinContent(j) + threshold);

    }

    if (totalCharge !=0){

        mean = mean/totalCharge;
        mean2 = mean2/totalCharge;
        rms2 = sqrt(mean2 - mean*mean) * sqrt(mean2-mean*mean);

    }

    return rms2;

}

void makePlot(){

  SetGenericStyle();
 TH1D* baseHisto =  (TH1D*)_file0->Get("Event_66/Event_66_TimeWfm_channel7000");
 TH1D* comparisonHisto = (TH1D*)_file0->Get("Event_66/Event_66_TimeWfm_channel7001");

 // open 3 tick window1 around peak
 TH1D* window1 = new TH1D("window1", "", 40, -20, 20);
 TH1D* window2 = new TH1D("window2", "", 40, -20, 20);

 for (int i = 0; i < window1->GetNbinsX(); i++){

    window1->SetBinContent(i, baseHisto->GetBinContent(baseHisto->GetMaximumBin()-20+i)+0.2);
    window2->SetBinContent(i, comparisonHisto->GetBinContent(comparisonHisto->GetMaximumBin()-20+i)+0.2);

 }

 window1->SetMaximum(baseHisto->GetMaximum()*2.2);
 window1->SetLineColor(kPTVibrantMagenta);
 window1->SetFillColor(kPTVibrantMagenta);
 window1->SetFillStyle(3345);
 window1->GetXaxis()->SetTitle("Time (ticks)");
 window1->GetYaxis()->SetTitle("Arb. Units");

 TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);

 for (int i = -5; i <= 5; i++){

   window1->GetXaxis()->CenterTitle();
   window1->GetYaxis()->CenterTitle();
   // make new histogram offset by a few
   TH1D* offset = new TH1D("offset", "",40, -20, 20);

   for (int j = 0; j < offset->GetNbinsX(); j++){
  
     offset->SetBinContent(j+i, window2->GetBinContent(j));

   }

   offset->SetLineColor(kPTVibrantCyan);
   offset->SetFillColor(kPTVibrantCyan);
   offset->SetFillStyle(3354);
 
   TH1D* sum = (TH1D*)offset->Clone("sum");
   sum->Add(window1);
   sum->SetLineColor(kBlack);
   sum->SetFillStyle(0);
   
   window1->GetYaxis()->SetRangeUser(0, sum->GetMaximum()*1.25);
   window1->Draw();
   offset->Draw("same");
   sum->Draw("same");

   double rms2 = getRms2(sum);

   TPaveText* pt = new TPaveText(0.6, 0.7, 0.85, 0.85, "NDC");
   pt->SetFillStyle(0);
   pt->SetBorderSize(0);
   TString strings = Form("RMS^{2}: %0.2f", rms2);
   pt->AddText(strings);
   
   pt->Draw("same");
 
   TString saveString = Form("offsetPlot_%i", i);
   ApplyLabel(kData);
   c1->SaveAs(saveString+".pdf");

 }

}
