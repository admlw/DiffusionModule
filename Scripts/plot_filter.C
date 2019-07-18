void plot_filter(){
  TF1 *filter_u = new TF1("filter_u","exp(-0.5*pow(x/[0],[1]))");
  double par[2]={1.43555e+01/200.*2.,4.95096e+00};
  filter_u->SetParameters(par);

  TF1 *filter_v = new TF1("filter_v","exp(-0.5*pow(x/[0],[1]))");
  double par1[2]={1.47404e+01/200.*2.,4.97667e+00};
  filter_v->SetParameters(par1);

  TF1 *filter_w = new TF1("filter_y","exp(-0.5*pow(x/[0],[1]))");
  double par2[2]={1.45874e+01/200.*2.,5.02219e+00};
  filter_w->SetParameters(par2);

  TF1 *filter_g = new TF1("filter_g","exp(-0.5*pow(x/[0],2))");
  double par3[1]={1.11408e-01};
  filter_g->SetParameters(par3);

  int nbin = 9600;
  TH1F *hfilter_u = new TH1F("hfilter_u","hfilter_u",nbin,0,nbin);
  TH1F *hfilter_v = new TH1F("hfilter_v","hfilter_v",nbin,0,nbin);
  TH1F *hfilter_w = new TH1F("hfilter_w","hfilter_w",nbin,0,nbin);
  TH1F *hfilter_g = new TH1F("hfilter_g","hfilter_g",nbin,0,nbin);
  
  for (Int_t i=0;i!=nbin;i++){
    
    Double_t frequency;
    frequency = 0;
    if (i!=0){
      frequency= (nbin/2.-fabs(i-nbin/2.))*2./nbin;
    }
    hfilter_u->SetBinContent(i+1,filter_u->Eval(frequency));
    hfilter_v->SetBinContent(i+1,filter_v->Eval(frequency));
    hfilter_w->SetBinContent(i+1,filter_w->Eval(frequency));
    hfilter_g->SetBinContent(i+1,filter_g->Eval(frequency));
  }

  TH1F *hfilter_time_u = new TH1F("hfilter_time_u","hfilter_time_u",nbin,-nbin/4.,nbin/4.);
  TH1F *hfilter_time_v = new TH1F("hfilter_time_v","hfilter_time_v",nbin,-nbin/4.,nbin/4.);
  TH1F *hfilter_time_w = new TH1F("hfilter_time_w","hfilter_time_w",nbin,-nbin/4.,nbin/4.);
  TH1F *hfilter_time_gaus = new TH1F("hfilter_time_gaus","hfilter_time_gaus",nbin,-nbin/4.,nbin/4.);

  double value_re[9600],value_im[9600];
  Int_t n = nbin;
  TVirtualFFT *ifft = TVirtualFFT::FFT(1,&n,"C2R M K");
  TH1 *fb;
  double baseline ;

  for (Int_t i=0;i!=nbin;i++){
    Double_t filter = hfilter_u->GetBinContent(i+1);
    value_re[i] = filter/nbin;
    value_im[i] = 0;
  }
  ifft->SetPointsComplex(value_re,value_im);
  ifft->Transform();  
  fb = 0;
  fb = TH1::TransformHisto(ifft,fb,"Re");
  
  //  fb->Draw();

  // baseline = -1./nbin;
  for (int i=0;i!=nbin;i++){
    Double_t content = fb->GetBinContent(i+1);
  //   content -=baseline;
    int binp = i+nbin/2+1;
    if (binp > nbin) binp -=nbin;
    hfilter_time_u->SetBinContent(binp,content);
  }
  // hfilter_time_u->Draw();

  // V-plane
  for (Int_t i=0;i!=nbin;i++){
    Double_t filter = hfilter_v->GetBinContent(i+1);
    value_re[i] = filter/nbin;
    value_im[i] = 0;
  }
  ifft->SetPointsComplex(value_re,value_im);
  ifft->Transform();  
  fb = 0;
  fb = TH1::TransformHisto(ifft,fb,"Re");
  
  baseline = -1./nbin;
  for (int i=0;i!=nbin;i++){
    Double_t content = fb->GetBinContent(i+1);
    //content -=baseline;
    int binp = i+nbin/2+1;
    if (binp > nbin) binp -=nbin;
    hfilter_time_v->SetBinContent(binp,content);
    
  }
  // W-plane
  for (Int_t i=0;i!=nbin;i++){
    Double_t filter = hfilter_w->GetBinContent(i+1);
    value_re[i] = filter/nbin;
    value_im[i] = 0;
  }
  ifft->SetPointsComplex(value_re,value_im);
  ifft->Transform();  
  fb = 0;
  fb = TH1::TransformHisto(ifft,fb,"Re");
  
  baseline = -1./nbin;
  for (int i=0;i!=nbin;i++){
    Double_t content = fb->GetBinContent(i+1);
    //content -=baseline;
    int binp = i+nbin/2+1;
    if (binp > nbin) binp -=nbin;
    hfilter_time_w->SetBinContent(binp,content);
  }
  
  // Gaussian one 
  for (Int_t i=0;i!=nbin;i++){
    Double_t filter = hfilter_g->GetBinContent(i+1);
    value_re[i] = filter/nbin;
    value_im[i] = 0;
  }
  ifft->SetPointsComplex(value_re,value_im);
  ifft->Transform();  
  fb = 0;
  fb = TH1::TransformHisto(ifft,fb,"Re");
  
  baseline = -1./nbin;
  for (int i=0;i!=nbin;i++){
    Double_t content = fb->GetBinContent(i+1);
    //content -=baseline;
    int binp = i+nbin/2+1;
    if (binp > nbin) binp -=nbin;
    hfilter_time_gaus->SetBinContent(binp,content);
  }


  TCanvas *c1 = new TCanvas("c1","c1", 600,600);
  hfilter_time_u->SetTitle("Time Domain");
  hfilter_time_u->Draw();
  hfilter_time_u->SetXTitle("t (#mus)");
  hfilter_time_u->SetYTitle("Filter(t)");
  hfilter_time_v->Draw("same");
  hfilter_time_w->Draw("same");
  hfilter_time_u->SetLineColor(2);
  hfilter_time_v->SetLineColor(4);
  hfilter_time_w->SetLineColor(6);
  //hfilter_time_u->GetYaxis()->SetRangeUser(-0.05,1.0);
  hfilter_time_u->GetYaxis()->SetRangeUser(-0.05,0.2);
  hfilter_time_u->GetXaxis()->SetRangeUser(-20,20);
  std::cout << "StdDev of w time filter: " << hfilter_time_w->GetStdDev() << std::endl;
  std::cout << "StdDev of w time filter: " << hfilter_time_w->GetStdDev() << std::endl;
 
  
  

  // TF1 *f1 = new TF1("f1","1./sqrt(2.*3.1415926)/[0]*exp(-x*x/2./[0]/[0])");
  // Double_t par3[1] = {3./2.1}; //
  // f1->SetParameters(par3);
  
  // TH1F *hfilter_time_gaus = (TH1F*)hfilter_time_u->Clone("hfilter_time_gaus");
  // hfilter_time_gaus->Reset();
  // for (Int_t i=0;i!=hfilter_time_gaus->GetNbinsX();i++){
  //   hfilter_time_gaus->SetBinContent(i+1,f1->Eval(hfilter_time_gaus->GetBinCenter(i+1)));
  // }
  // hfilter_time_gaus->Scale(1./hfilter_time_gaus->GetSum());
  hfilter_time_gaus->Draw("same");
  std::cout << "StdDev of gaus time filter: " << hfilter_time_gaus->GetStdDev() << std::endl;
  hfilter_time_gaus->SetLineColor(1);

  hfilter_time_u->SetLineWidth(2);
  hfilter_time_v->SetLineWidth(2);
  hfilter_time_w->SetLineWidth(2);
  hfilter_time_gaus->SetLineWidth(2);
  
  hfilter_time_u->SetLineStyle(2);
  hfilter_time_v->SetLineStyle(3);
  hfilter_time_w->SetLineStyle(4);
  hfilter_time_gaus->SetLineStyle(1);


  // hfilter_time_u->Rebin(4);
  // hfilter_time_v->Rebin(4);
  // hfilter_time_w->Rebin(4);
  // hfilter_time_gaus->Rebin(4);

  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->SetFillColor(10);
  le1->AddEntry(hfilter_time_u,"U-plane","l");
  le1->AddEntry(hfilter_time_v,"V-plane","l");
  le1->AddEntry(hfilter_time_w,"W-plane","l");
  le1->AddEntry(hfilter_time_gaus,"Gaus","l");
  le1->Draw();

  //std::cout << hfilter_time_u->GetSum() << std::endl;

  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  TH1 *hmu = 0;
  hmu = hfilter_time_u->FFT(hmu,"MAG");

  TH1 *hmv = 0;
  hmv = hfilter_time_v->FFT(hmv,"MAG");
  
  TH1 *hmw = 0;
  hmw = hfilter_time_w->FFT(hmw,"MAG");
  
  TH1 *hmg = 0;
  hmg = hfilter_time_gaus->FFT(hmg,"MAG");
  

  TH1F *hmu1 = new TH1F("hmu1","hmu1",4800,0,1);
  TH1F *hmv1 = new TH1F("hmv1","hmv1",4800,0,1);
  TH1F *hmw1 = new TH1F("hmw1","hmw1",4800,0,1);
  TH1F *hmg1 = new TH1F("hmg1","hmg1",4800,0,1);

  for (Int_t i=0;i!=4800;i++){
    hmu1->SetBinContent(i+1,hmu->GetBinContent(i+1));
    hmv1->SetBinContent(i+1,hmv->GetBinContent(i+1));
    hmw1->SetBinContent(i+1,hmw->GetBinContent(i+1));
    hmg1->SetBinContent(i+1,hmg->GetBinContent(i+1));
  }


  hmu1->Draw();
  hmu1->SetXTitle("#omega (MHz)");
  hmu1->SetYTitle("Filter(#omega)");
  hmv1->Draw("same");
  hmw1->Draw("same");
  hmg1->Draw("same");
  
  hmu1->SetLineColor(2);
  hmv1->SetLineColor(4);
  hmw1->SetLineColor(6);
  hmg1->SetLineColor(1);

  hmu1->SetLineWidth(2);
  hmv1->SetLineWidth(2);
  hmw1->SetLineWidth(2);
  hmg1->SetLineWidth(2);
  
  hmu1->SetLineStyle(2);
  hmv1->SetLineStyle(3);
  hmw1->SetLineStyle(4);
  hmg1->SetLineStyle(1);

  hmu1->SetTitle("Frequency Domain");
  
  
}
