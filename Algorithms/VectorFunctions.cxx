#include "VectorFunctions.h"

namespace diffmod {

    TVector3 VectorFunctions::getMomentumVector(const recob::Track& a) {
      auto const& smv = a.StartMomentumVector();
      TVector3 aMom(smv.X(), smv.Y(), smv.Z() );

      return aMom;
    }

    double VectorFunctions::getAngle(const recob::Track &track, diffmod::VectorFunctions vecFuncs, std::string proj) {

      const double rad2deg = 57.295827; // Convert to degrees

      //TVector3 const& dir_start = track.Direction();
      recob::Track::Vector_t dir_start = track.StartDirection();
      double angle;

      if ((proj == "xz") | (proj == "zx")){
        angle = std::atan2(dir_start.X(), dir_start.Z() ) * rad2deg;
      }

      else if ((proj == "xy") | (proj == "yx")){
        angle = std::atan2(dir_start.Y(), dir_start.Z() ) * rad2deg;
      }
      else {
        throw std::invalid_argument("Valid arguments are 'no', 'xy', 'xz' and 'yz'");
      }

      //if (angle > 90) return 180-angle;
      // else
      return angle;
    }

    TVector3 VectorFunctions::getUnitVector(TVector3 a) {
      return a.Unit();
    }

    double VectorFunctions::getDotProduct(TVector3 a, TVector3 b) {
      return a.Dot(b);
    }

}
