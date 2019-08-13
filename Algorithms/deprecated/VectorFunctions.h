#ifndef VECTORFUNCTIONS_H
#define VECTORFUNCTIONS_H

// larsoft includes
#include "lardataobj/RecoBase/Track.h"

// root includes
#include "TVector3.h"

namespace diffmod {

    class VectorFunctions {
        public: 

            TVector3 getMomentumVector(const recob::Track& a);

            double getAngle(const recob::Track& track, VectorFunctions vecFuncs, std::string proj);

            TVector3 getUnitVector(TVector3 a);

            double getDotProduct(TVector3 a, TVector3 b);
    };

}

#endif
