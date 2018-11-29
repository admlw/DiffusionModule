#ifndef WAVEFORMFUNCTIONS_H
#define WAVEFORMFUNCTIONS_H

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

    class WaveformFunctions {
        public: 
          bool passesHitSelection(recob::Hit const* hit, double HIT_GOODNESSOFFIT_CUT); 

          double convertXToTicks(double xPosition, int WAVEFORM_DRIFT_START_TICK, int WAVEFORM_DRIFT_SIZE, double X_WIDTH);

          TH1D* applyGlobalBaselineCorrection(TH1D *h_rawD, TH1D *h_rawDCorrected);

          double findXCorrection(WaveformFunctions waveFuncs, TH1D *summedWaveform, TH1D *h, int NUMBER_TICKS_PER_BIN, double mean);

          std::vector<double> getSigma(TH1D *h_rawDCorrected);

          double convertTicksToX(int tick, int WAVEFORM_DRIFT_START_TICK, int WAVEFORM_DRIFT_SIZE, double X_WIDTH);

          double getMedian(TH1D *h);

          double getRms2(TH1D *h);

    };

}

#endif
