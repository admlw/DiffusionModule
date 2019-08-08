#include <iostream>
void makeSigmaDistPlots() {
  std::ifstream in("vals.txt");
  if (!in) {
    std::cout << "Bad input file" << std::endl;
    return;
  }

  // Get median, max and truncated mean values (from module output)
  int const nlines = 25;
  double medians[nlines];
  double maxs[nlines];
  double truncMeans[nlines];

  for (int i = 0; i < nlines; i++) {
    in >> medians[i];
    in >> maxs[i];
    in >> truncMeans[i];
    /*
    std::cout << "Medians at " << i << ": " << medians[i] << std::endl;
    std::cout << "Maxs at " << i << ": " << maxs[i] << std::endl;
    std::cout << "Truncs at " << i << ": " << truncMeans[i] << std::endl;
    */
  }

  TString dir = "/uboone/app/users/amogan/mcc9_diffusion/workdir/";
  TFile *fin = new TFile(dir+"sigma_map.root", "READ");
  //std::vector<TH1D*> h_sigma_hists(nlines);
  TH1D *h_sigma;
  TCanvas *c = new TCanvas("c", "c", 750, 550);
  gStyle->SetOptStat(0);
  for (int i = 0; i < nlines; i++) {
    TString sigmaHistName = Form("DiffusionModule/h_sigma_%i", i);
    //std::cout << "Hist name " << sigmaHistName << std::endl;
    //h_sigma_hists.at(i) = (TH1D*)fin->Get(sigmaHistName);
    h_sigma = (TH1D*)fin->Get(sigmaHistName);
    TString sigmaHistAxisName = Form("#sigma_{t}^{2} Bin %i", i);
    h_sigma->GetXaxis()->SetTitle(sigmaHistAxisName);
    h_sigma->GetXaxis()->SetTitleSize(0.04);
    h_sigma->GetXaxis()->SetRangeUser(2., 6.);
    h_sigma->Draw();

    TLine *l1 = new TLine(medians[i], 0, medians[i], h_sigma->GetMaximum() );
    TLine *l2 = new TLine(maxs[i], 0, maxs[i], h_sigma->GetMaximum() );
    TLine *l3 = new TLine(truncMeans[i], 0, truncMeans[i], h_sigma->GetMaximum() );
    l1->SetLineColor(kRed+1);
    l1->SetLineWidth(2);
    l1->SetLineStyle(1);
    l2->SetLineColor(kAzure+2);
    l2->SetLineWidth(2);
    l2->SetLineStyle(9);
    l3->SetLineColor(kGreen+1);
    l3->SetLineWidth(2);
    l3->SetLineStyle(2);
    l1->Draw("same NDC");
    l2->Draw("same NDC");
    l3->Draw("same NDC");

    TLegend *l = new TLegend(0.6, 0.7, 0.9, 0.9);
    l->AddEntry(l1, "Median", "l");
    l->AddEntry(l2, "Max", "l");
    l->AddEntry(l3, "Trunc. Mean", "l");
    l->Draw("same");

    TString plotName = Form("other_plots/plot_sigma_%i.pdf", i);
    c->SaveAs(plotName, "PDF");

  }


}



























