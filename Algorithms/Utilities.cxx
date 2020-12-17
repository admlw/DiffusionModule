#include "Utilities.h"

namespace diffmod {

  double Utilities::getReturnerAngle(double bareAngle){

    double returnerAngle;
    if (bareAngle < 90)
      returnerAngle = bareAngle;
    else 
      returnerAngle = 180-bareAngle;

    return returnerAngle;

  }

  double Utilities::getThetaXZ(art::Ptr<recob::Track> thisTrack){

    recob::Track::Vector_t trkDir = thisTrack->StartDirection();

    double bareThetaXZ = std::abs(std::atan2(trkDir.X(), trkDir.Z())) *
                         180. / TMath::Pi();

    // We don't actually want the bare theta_xz, we're really interested
    // in the angle with respect to the beam direction, so
    double returnerThetaXZ = this->getReturnerAngle(bareThetaXZ);

    return returnerThetaXZ;

  }

  double Utilities::getThetaYZ(art::Ptr<recob::Track> thisTrack){
    recob::Track::Vector_t trkDir = thisTrack->StartDirection();

    double bareThetaYZ = std::abs(std::atan2(trkDir.Y(), trkDir.Z())) *
                         180. / TMath::Pi();

    // We don't actually want the bare theta_xz, we're really interested
    // in the angle with respect to the beam direction, so
    double returnerThetaYZ = this->getReturnerAngle(bareThetaYZ);

    return returnerThetaYZ;
  }

  /*
  double Utilities::FindMaximumValue(TH1D* h){
  
    double maxBinValue = h->GetBinCenter((h->GetMaximumBin()));
    
    TF1* gaus = new TF1("f1", "gaus", 0, 100);
    gaus->SetParameter(0, 1);
    gaus->SetParameter(1,maxBinValue);
    gaus->SetParameter(2, 0.25);
    h->Fit(gaus,"","", maxBinValue-0.25, maxBinValue+0.25);

    double value = (double)h->GetFunction("f1")->GetParameter("Mean");

    return value;

  }
  */

  // Michelle's version
  double Utilities::FindMaximumValue(TH1D* h){

    float qq[25]; float eq[25];
    float   dt[25],edt[25];
    float   kk[25],ek[25];

    int iwlow = 0; int iwhigh=0;
    int npts = 0;
    TH1F dq[25];  TH1F dtime[25];

    gStyle->SetOptFit(1111);
    //TCanvas* c1 = new TCanvas("c1");

      //TString label = Form("DiffusionModule/plane%1d/h_sigma_%1d_plane%1d",iplane,i,iplane);
      //dq[i] = *((TH1F*)f->Get(label)); 
      dq[i] = h; 
      //std::cout << dq[i].GetEntries() << std::endl;
      //  Could be 50000  . . . see if uncertainty reflects this
      if (dq[i].GetEntries()>2000) {
        int mbin = dq[i].GetMaximumBin();
        float xmid = dq[i].GetBinCenter(mbin);
        float maxh=dq[i].GetBinContent(mbin);
        // dynamic width
        float target=0.5*maxh;
        int thisbin=mbin-1;  float thish=dq[i].GetBinContent(thisbin);
        while (thish>target) {
          thisbin--;
          thish=dq[i].GetBinContent(thisbin);
        }
        iwlow = thisbin;
        std::cout << "[MAXFIT] lowedge " << iwlow << " center " << mbin << std::endl;
        thisbin=mbin+1;  thish=dq[i].GetBinContent(thisbin);
        while (thish>target) {
          thisbin++;
          thish=dq[i].GetBinContent(thisbin);
        }
        iwhigh = thisbin;
        std::cout << "[MAXFIT] highedge " << iwhigh << " center " << mbin << std::endl;
        /// end dynamic width
        float lowedge = dq[i].GetBinCenter(iwlow);
        float highedge = dq[i].GetBinCenter(iwhigh);
        TF1 *f1 = new TF1("f1", "gaus",lowedge,highedge);
        f1->SetParameter(0,xmid);
        f1->SetParameter(1,xmid);
        float temp = xmid-lowedge;
        f1->SetParameter(2,temp);
        dq[i].Fit("f1","R");
        qq[npts] = f1->GetParameter(1);
        eq[npts] = f1->GetParError(1);
        dt[npts]=i; edt[npts]=0;
        kk[npts] = f1->GetParameter(2);
        ek[npts] =f1->GetParError(2);
        //      edt[npts]= 69./2/1.414;  // uncertainty = full width of a strip in drift time divided by 2*sqrt(2)
        npts++;
      }
      else {
        // qq[i]=1500.0;0
        // eq[i]=1500.0;
        // dt[i]=(i+0.5)*1250.0/18.;
        // edt[i]=1000.0;
      }
      
      TString name = outdir+Form("fit_plane%1d_bin%2d",iplane,i)+".png";
      c1->Print(name);
      c1->Update();
    //    std::cout << "dq " << i << " " << qq[i] << " " << dt[i] << std::endl;


      std::cout << "npts " << npts << std::endl;
      TCanvas* c2 = new TCanvas("c2");
      TGraphErrors *gr = new TGraphErrors(npts,dt,qq,edt,eq);
      gr->SetMarkerColor(4);
      gr->SetMarkerStyle(21);
      gr->GetXaxis()->SetTitle("drift time bin");
      gr->GetYaxis()->SetTitle("peak position");
      gr->Draw("AP");
      c2->Print(outdir+Form("peak_plane%1d",iplane)+".png");
      c2->Update();
      

      TGraphErrors *gr2 = new TGraphErrors(npts,dt,kk,edt,ek);
      gr2->SetMarkerColor(4);
      gr2->SetMarkerStyle(21);
      gr2->GetXaxis()->SetTitle("drift time bin");
      gr2->GetYaxis()->SetTitle("fit width");
      gr2->Draw("AP");

      c2->Print(outdir+Form("width_plane%1d",iplane)+".png");
      c2->Update();
      

      f->Close();

  }

}
