#include "TrackFunctions.h"

bool passesPreSelection(recob::Track const& track,
                        VectorFunctions vecFuncs,
                        TrackFunctions trackFuncs,
                        double trackAngleXZLowBound,
                        double trackAngleXZHighBound,
                        double trackAngleYZLowBound,
                        double trackAngleYZHighBound )
{
    double trackLength = track.Length();
    double trackBeginningEndAngle = trackFuncs.findTrackStraightness(track, vecFuncs);

    TVector3 trackMomentum = vecFuncs.getMomentumVector(track);
    TVector3 zDir(0, 0, 1);

    float trackAngleXZ = 

}
