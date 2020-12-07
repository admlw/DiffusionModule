// use the well defined atrazhev-timoshkin regions and interpolate between the 
// two, providing a functional form which is close enough

//.............................................................................
// returns e_l for well defined regions as per altrazhev-timoshkin paper
// values of 10 and 1200 are chosen to get a reasonable fit and for no
// other reason
Double_t atrazhev_el_from_paper(double ef){

  //1 degree kelvin = 8.621738 X10-5  eV
  float T_eV = 89 * 0.000086173;
  float Eh = 73*2; // V/cm

  if( ef < 10)
    return T_eV;
  else if (ef > 1200)
    return (0.5 * 0.8 * T_eV * (ef / Eh));
  else return 0;

}


//.............................................................................
// main
void fit_atrazhev(){

  TCanvas* c1 = new TCanvas("c1", "", 500, 500);
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.12);
  c1->SetGridy();
  c1->SetGridx();
  c1->SetLogx();
  c1->SetLogy();

  TF1* atrazhev_el = new TF1("atrazhev_el", "atrazhev_el_from_paper(x)", 10e-2, 10000);
  const int npoints = 100000;
  Double_t xs[npoints] = {};
  Double_t ys[npoints] = {};
  float min = 10e-2;
  float max = 10000;
  float step = (max-min)/npoints;
  for (int i = 0; i < npoints; ++i){
    xs[i] = min + (i*step);
    ys[i] = atrazhev_el_from_paper(xs[i]);
    if (ys[i] == 0){
      xs[i] = 0.1;
      ys[i] = 89*0.000086173;
    }
  }
  TGraph* bg = new TGraph(npoints, xs, ys);
  bg->SetTitle("");
  bg->GetXaxis()->SetTitle("E (V/cm)");
  bg->GetYaxis()->SetTitle("Electron Energy, #epsilon_{L} (eV)");
  bg->GetXaxis()->CenterTitle();
  bg->GetYaxis()->CenterTitle();

  bg->Draw("pa");
  TF1* fitter = new TF1("fitter", "(89*0.000086173) + pol4", 1, 10000);
  bg->Fit("fitter");
  c1->SaveAs("AtrazhevFit.pdf");
  c1->SaveAs("AtrazhevFit.png");

}
