#include "VectorFunctions.h"

namespace diffmod {

    // Because StartMomentumVector() is a weird type by default
    TVector3 VectorFunctions::getMomentumVector(const recob::Track& a) {
      auto const& smv = a.StartMomentumVector();
      TVector3 aMom(smv.X(), smv.Y(), smv.Z() );

      return aMom;
    }

    // Get theta_xz or theta_yz, depending on passed string parameter
    double VectorFunctions::getAngle(const recob::Track &track, diffmod::VectorFunctions vecFuncs, std::string proj) {

      const double rad2deg = 57.295827; // Convert to degrees

      recob::Track::Vector_t dir_start = track.StartDirection();
      double angle;

      if ((proj == "xz") | (proj == "zx")){
        angle = std::atan2(dir_start.X(), dir_start.Z() ) * rad2deg;
      }

      else if ((proj == "xy") | (proj == "yx")){
        angle = std::atan2(dir_start.Y(), dir_start.Z() ) * rad2deg;
      }
      else {
        throw std::invalid_argument("Valid arguments are 'xz' and 'yz'");
      }

      //if (angle > 90) return 180-angle;
      // else
      return angle;
    }

    // You'll never guess what this function does
    TVector3 VectorFunctions::getUnitVector(TVector3 a) {
      return a.Unit();
    }

    // Or this one
    double VectorFunctions::getDotProduct(TVector3 a, TVector3 b) {
      return a.Dot(b);
    }

}
