#include <iostream>
void make_sigma_dist_plots() {
  // Input txt file comes from output of ana stage in diffusion module
  // Just look for three columns of 25 numbers all in a row
  // Column structure is median \t max \t truncated mean
  //std::ifstream in("vals.txt");
  std::ifstream in("vals_single_muons_Aug2020_fwdgoing.txt");
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

  TString dir = "/pnfs/uboone/persistent/users/amogan/v08_00_00_25/diffusion_output_files/diffusionSigmamap/";
  //TFile *fin = new TFile(dir+"sigma_map_run3_crt_Aug2020_newFV.root", "READ");
  TFile *fin = new TFile(dir+"sigma_map_single_muons_Aug2020_newFV.root", "READ");
  if (!fin) {
    std::cout << "ERROR: Bad input file" << std::endl;
    return;
  }

  TH1D *h_sigma;
  TCanvas *c = new TCanvas("c", "c", 750, 550);
  gStyle->SetOptStat(0);
  for (int i = 0; i < nlines; i++) {
    TString sigmaHistName = Form("DiffusionModule/plane2/h_sigma_%i_plane2", i);
    //std::cout << "Hist name " << sigmaHistName << std::endl;
    //h_sigma_hists.at(i) = (TH1D*)fin->Get(sigmaHistName);
    h_sigma = (TH1D*)fin->Get(sigmaHistName);
    TString sigmaHistAxisName = Form("#sigma_{t}^{2} Bin %i", i);
    h_sigma->GetXaxis()->SetTitle(sigmaHistAxisName);
    h_sigma->GetXaxis()->SetTitleSize(0.04);
    h_sigma->GetXaxis()->SetRangeUser(2., 6.);
    h_sigma->Draw();

    TLine *l1 = new TLine(medians[i],    0, medians[i],    h_sigma->GetMaximum() );
    TLine *l2 = new TLine(maxs[i],       0, maxs[i],       h_sigma->GetMaximum() );
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



