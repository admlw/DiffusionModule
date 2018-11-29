#include "TrackFunctions.h"

namespace diffmod {
    bool passesPreSelection(recob::Track const& track,
                            diffmod::VectorFunctions vecFuncs,
                            diffmod::TrackFunctions trackFuncs,
                            double trackAngleXZLowBound,
                            double trackAngleXZHighBound,
                            double trackAngleYZLowBound,
                            double trackAngleYZHighBound )
    {
        double trackLength = track.Length();
        if (trackLength < 25) {
            return false;
        }

        //double trackBeginningEndAngle = trackFuncs.findTrackStraightness(track, vecFuncs);

        float trackAngleXZ = vecFuncs.getAngle(track, vecFuncs, "xz");
        float trackAngleYZ = vecFuncs.getAngle(track, vecFuncs, "yz");

        if ( trackAngleXZ < trackAngleXZLowBound || trackAngleXZ > trackAngleXZHighBound) {
          //std::cout << "Angle " << trackAngleXZ << " out of XZ bound" << std::endl;   
          return false;
        }
              
        else if (trackAngleYZ < trackAngleYZLowBound || trackAngleYZ > trackAngleYZHighBound) {
          //std::cout << "Angle " << trackAngleYZ << " out of YZ bound" << std::endl;   
          return false;
        }
        else if (trackLength < 25){
          //std::cout << "Track too short" << std::endl;
          return false;
        }
        else {
          //std::cout << "Track passed" << std::endl;    
          return true;
        }

    }

}
