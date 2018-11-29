#ifndef VECTORUTILS_H
#define VECTORUTILS_H

// larsoft includes
#include "lardataobj/RecoBase/Track.h"

// root includes
#include "TVector3.h"

namespace diffmod {

    class VectorFunctions {
        public: 

            TVector3 getMomentumVector(const recob::Track& a);

            double getAngle(TVector3 a, TVector3 b);

            TVector3 getUnitVector(TVector3 a, TVector3 b);

            double getDotProduct(TVector3 a, TVector3 b);
    };

}

#endif
