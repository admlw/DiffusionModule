#include "StylePlots.h"

void style(TGraph* g, std::string plane){

  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.6);
  g->SetLineWidth(2);
  g->GetXaxis()->CenterTitle();
  g->GetYaxis()->CenterTitle();

  if (plane == "u"){
    g->SetMarkerColor(kAzure+1);
    g->SetLineColor(kAzure+1);
  }
  else if (plane == "v"){
    g->SetMarkerColor(kGreen+1);
    g->SetLineColor(kGreen+1);
  }
  else{

  }
}

void area_normalize(TH1D* h){
  h->Sumw2();
  h->Scale(1./h->Integral());
}

std::vector<double> generate_times(){
  std::vector<double> times;
  for (int i = 0; i < 25; ++i){
    times.push_back(46+(i*92));
  }
  return times;
}

std::vector<double> get_histograms_and_fit(TFile* f, std::string plane){
  std::vector<double> widths;
  for (int i = 0; i < 25; i++){
    std::string name = "h"+plane+std::to_string(i);
    TH1D* h = (TH1D*)f->Get(name.c_str());
    area_normalize(h);

    float max = h->GetBinCenter(h->GetMaximumBin());
    float width = 0.1;

    TF1* gaus = new TF1("gaussfit", "gaus", 0, 5);
    gaus->SetParameter(0,1);
    gaus->SetParameter(1, max);
    gaus->SetParameter(2, width);

    TCanvas* c1 = new TCanvas();
    h->Fit(gaus, "", "", max-width, max+width);
    h->Draw();
    c1->SaveAs((name+".pdf").c_str());

    widths.push_back(h->GetFunction("gaussfit")->GetParameter(2)/h->GetFunction("gaussfit")->GetParameter(1));
  }
  return widths;
}

void ana_hit_dists_file(){

  SetGenericStyle();  

  TFile* f = new TFile("hit_width_dists.root", "read");


  std::vector<double> times = generate_times();
  std::vector<double> uwidths = get_histograms_and_fit(f, "u");
  std::vector<double> vwidths = get_histograms_and_fit(f, "v");
  std::vector<double> ywidths = get_histograms_and_fit(f, "y");

  TCanvas* c2 = new TCanvas();
  TGraph* drifttimewidthu = new TGraph(25, &times[0], &uwidths[0]);
  TGraph* drifttimewidthv = new TGraph(25, &times[0], &vwidths[0]);
  TGraph* drifttimewidthy = new TGraph(25, &times[0], &ywidths[0]);

  TH2D* bg = new TH2D("bg", ";Drift Time (#mus);Hit RMS Width/Mean", 1, 0, 2300, 1, 0.05,  0.15);
  bg->GetXaxis()->CenterTitle();
  bg->GetYaxis()->CenterTitle();
  bg->Draw();

  style(drifttimewidthu, "u");
  style(drifttimewidthv, "v");
  style(drifttimewidthy, "y");

  drifttimewidthu->Draw("pl same");
  drifttimewidthv->Draw("pl same");
  drifttimewidthy->Draw("pl same");

  TLatex *upl = new TLatex(0.15, 0.80, "MicroBooNE Data");
  upl->SetNDC();
  upl->SetTextSize(0.05);
  upl->Draw();

 // drifttimewidthu->Fit("pol1", "", "", 1000, 2300);
 // drifttimewidthv->Fit("pol1", "", "", 1000, 2300);
 // drifttimewidthy->Fit("pol1", "", "", 1000, 2300);

  TLegend* leg = new TLegend(0.6, 0.7, 0.85, 0.85);
  leg->AddEntry(drifttimewidthu, "U Plane");
  leg->AddEntry(drifttimewidthv, "V Plane");
  leg->AddEntry(drifttimewidthy, "Y Plane");
  leg->Draw("same");

  c2->SaveAs("graph.pdf");

}
