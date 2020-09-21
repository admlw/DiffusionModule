// Plot two waveforms from different files on the same Canvas
// Purpose is to compare DT up waveforms with CV (or similar)
// Takes two files that were generated using the same gen and G4
// fcl files, then run through detsims with DT varied (everything
// else the same)
// Written by A. Mogan Sep. 2, 2020

void make_wvfm_comp_plot() {

  //TFile *fin_cv = new TFile("out_DT_default.root"); // DT = 9.8
  //TFile *fin_up = new TFile("out_DT_way_way_up.root"); // DT = 100.0
  TString dir = "/uboone/data/users/amogan/wireplot_outfiles_diffusion/";
  TFile *fin_cv = new TFile(dir+"out_DT_off_thetaXZ_20.root"); // DT = 0.0, high-angle
  TFile *fin_up  = new TFile(dir+"out_DT_way_way_up_thetaXZ_20.root"); // DT = 100.0, high-angle
  
  // Initialization
  // Specify event and channel number of desired waveform
  // This can be found by looking at the output of GetWaveform.cc
  // and searching for events with peak above threshold (~2050 for 
  // U and V planes, ~480 for Y plane)
  
  Int_t event = 3;
  Int_t channel = 7432;

  TString histname = Form("Event_%i/Event_%i_TimeWfm_channel%i", event, event, channel); 

  TCanvas *c = new TCanvas("c");

  fin_cv->cd();
  TH1D *h_cv = (TH1D*)fin_cv->Get(histname);
  int maxbin = h_cv->GetMaximumBin();
  double maxval = h_cv->GetBinCenter(maxbin);
  //h_cv->GetXaxis()->SetRangeUser(maxval-100., maxval+100.);
  h_cv->GetXaxis()->SetTitle("Time (Ticks)");
  h_cv->GetYaxis()->SetTitle("ADC");
  h_cv->SetLineColor(kAzure+2);
  //h_cv->Scale(1.0/h_cv->Integral() );
  h_cv->Draw("hist");

  fin_up->cd();
  TH1D *h_up = (TH1D*)fin_up->Get(histname);
  h_up->SetLineColor(kRed+1);
  //h_up->Scale(1.0/h_up->Integral() );
  h_up->Draw("hist same");

  TLegend *l = new TLegend(0.6, 0.7, 0.85, 0.9);
  l->AddEntry(h_cv, "D_{T} = 0.0 cm^{2}/s", "l");
  l->AddEntry(h_up, "D_{T} = 100.0 cm^{2}/s", "l");
  //l->AddEntry(h_cv, "D_{T} = 9.80 cm^{2}/s", "l");
  //l->AddEntry(h_up, "D_{T} = 100.0 cm^{2}/s", "l");
  //l->Draw("same");

  TPaveText *pt = new TPaveText(0.15, 0.75, 0.4, 0.9, "NDC");
  pt->SetTextAlign(11);
  pt->SetFillStyle(0);
  pt->SetLineWidth(0);
  //pt->SetTextFont(42);
  //pt->SetTextSize(18);
  TString event_text   = Form("Event %i"  , event);
  TString channel_text = Form("Channel %i", channel);
  pt->AddText("MicroBooNE Simulation");
  pt->AddText("#theta_{xz} = 20 deg.");
  pt->AddText(event_text);
  pt->AddText(channel_text);
  pt->Draw("same");

  TString outname = Form("Event_%i_TimeWfm_channel%i", event, channel); 
  c->SaveAs(outname+".pdf", "PDF");


}



