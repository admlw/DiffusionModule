#ifndef TRACKFUNCTIONS_H
#define TRACKFUNCTIONS_H

// larsoft includes
#include "lardataobj/RecoBase/Track.h"

// root includes
#include "TVector3.h"

// Local includes
#include "VectorFunctions.h"

namespace diffmod {

    class TrackFunctions {
        public:
            bool passesPreSelection(recob::Track const& track,
                 VectorFunctions vecFuncs,
                 TrackFunctions trackFuncs,
                 double trackAngleXZLowBound,
                 double trackAngleXZHighBound,
                 double trackAngleYZLowBound,
                 double trackAngleYZHighBound
            );

            double findTrackStraightness(recob::Track const& track,
                   TrackFunctions trackFuncs
            );
    };

}

#endif
