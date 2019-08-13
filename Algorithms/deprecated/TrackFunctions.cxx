#include "TrackFunctions.h"

namespace diffmod {
    bool passesPreSelection(recob::Track const& track,
                            diffmod::VectorFunctions vecFuncs,
                            diffmod::TrackFunctions trackFuncs,
                            double trackAngleXZLowBound,
                            double trackAngleXZHighBound,
                            double trackAngleYZLowBound,
                            double trackAngleYZHighBound,
                            double trackLengthLowBound)
    {

        //double trackBeginningEndAngle = trackFuncs.findTrackStraightness(track, vecFuncs);

        double trackLength = track.Length();
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
        else if (trackLength < trackLengthLowBound){
          //std::cout << "Track too short" << std::endl;
          return false;
        }
        else {
          //std::cout << "Track passed" << std::endl;    
          return true;
        }

    }

}
