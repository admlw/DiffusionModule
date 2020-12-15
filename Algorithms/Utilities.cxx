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

}
