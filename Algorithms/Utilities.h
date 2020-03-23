#ifndef DIFFMODUTILITIES_H
#define DIFFMODUTILITIES_H

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"

// larsoft includes
#include "lardataobj/RecoBase/Track.h"

// root includes
#include "TMath.h"

namespace diffmod{

  class Utilities {

    public:

      /// in order to make things easier to interperet, we 
      /// return any angle > 90 degrees as 180-angle
      /// this means forward and backward angles are treated
      /// in the same way
      double getReturnerAngle(double bareAngle);

      /// method gets theta xz from track
      double getThetaXZ(art::Ptr<recob::Track> thisTrack);

      /// method gets theta yz from track
      double getThetaYZ(art::Ptr<recob::Track> thisTrack);

  };

}

#endif
