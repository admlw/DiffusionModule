#ifndef DIFFUSIONFUNCTIONS_H
#define DIFFUSIONFUNCTIONS_H

// larsoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"

// root includes
#include "TH1.h"
#include "TF1.h"
#include "TVector.h"

// cpp includes
#include <numeric>
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>

namespace diffmod {
    class diffusionFunctions{

      public:

        std::vector<double> getMomentumVector(const recob::Track& a);

        float getAngle(const std::vector<double> a, const std::vector<double> b, diffusionFunctions diffutil, std::string proj);

        std::vector<float> getUnitVector(std::vector<double> a);

        float getDotProduct(std::vector<float>, std::vector<float> b);

        bool passesPreSelection(recob::Track const& track, diffusionFunctions _diffusionFunctions_instance, double trackAngleXZLowBound, double trackAngleXZHighBound, double trackAngleYZLowBound, double trackAngleYZHighBound);

        bool passesHitSelection(recob::Hit const* hit, double HIT_GOODNESSOFFIT_CUT);

        std::vector<double> getUnitVector(recob::tracking::Vector_t vec3);

        double getDotProduct(std::vector<double> a, std::vector<double> b);

        double findTrackStraightness(recob::Track const& track, diffusionFunctions _diffusionFunctions_instance);

        double convertXToTicks(double xPosition, int WAVEFORM_DRIFT_START_TICK, int WAVEFORM_DRIFT_SIZE, double X_WIDTH);

        TH1D* applyGlobalBaselineCorrection(TH1D* h, TH1D* h_corrected);

        double findXCorrection(diffusionFunctions _diffusionFunctions_instance, TH1D* summedWaveform, TH1D* h, int NUMBER_TICKS_PER_BIN, double mean);

        std::vector<double> getSigma(TH1D* h);

        double convertTicksToX(int tick, int WAVEFORM_DRIFT_START_TICK, int WAVEFORM_DRIFT_SIZE, double X_WIDTH);

        double getMedian(TH1D* h);

        double getRms2(TH1D* h);

    };
}

#endif
