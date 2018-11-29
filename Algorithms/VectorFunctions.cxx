#include "VectorFunctions.h"

TVector3 getMomentumVector(const recob::Track& a) {
  auto const& smv = a.StartMomentumVector();
  TVector3 aMom(smv.X(), smv.Y(), smv.Z() );

  return aMom;
}

double getAngle(TVector3 a, TVector3 b) {
  return a.Angle(b);
}

TVector3 getUnitVector(TVector3 a) {
  return a.Unit();
}

double getDotProduct(TVector3 a, TVector3 b) {
  return a.Dot(b);
}
