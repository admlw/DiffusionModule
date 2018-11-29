#include "VectorFunctions.h"

namespace diffmod {
    TVector3 VectorFunctions::getMomentumVector(const recob::Track& a) {
      auto const& smv = a.StartMomentumVector();
      TVector3 aMom(smv.X(), smv.Y(), smv.Z() );

      return aMom;
    }

    double VectorFunctions::getAngle(TVector3 a, TVector3 b, diffmod::VectorFunctions vecFuncs, std::string proj) {

      const double rad2deg = 57.295827; // Convert to degrees
      TVector3 aMom = a;
      TVector3 bMom = b;

      /*
      if (proj == "no"){
        aMom = (a.X(), a.Y(), a.Z() );
        bMom = (b.X(), b.Y(), b.Z() );
      }
      */

      if ((proj == "xz") | (proj == "zx")){
        aMom.SetY(0);
        bMom.SetY(0);
      }

      else if ((proj == "xy") | (proj == "yx")){
        aMom.SetZ(0);
        bMom.SetZ(0);
      }

      else if ((proj == "yz") | (proj == "zy")){
        aMom.SetX(0);
        bMom.SetX(0);
      }
      else
        throw std::invalid_argument("Valid arguments are 'no', 'xy', 'xz' and 'yz'");

      TVector3 aMomUnit = vecFuncs.getUnitVector(aMom);
      TVector3 bMomUnit = vecFuncs.getUnitVector(bMom);

      //float angle = std::acos(vecFuncs.getDotProduct(aMomUnit, bMomUnit)) * 180 / 3.14159;
      float angle = std::acos(aMomUnit.Dot(bMomUnit) ) * rad2deg;

      if (angle > 90) return 180-angle;
      else return angle;
    }

    TVector3 VectorFunctions::getUnitVector(TVector3 a) {
      return a.Unit();
    }

    double VectorFunctions::getDotProduct(TVector3 a, TVector3 b) {
      return a.Dot(b);
    }

}
