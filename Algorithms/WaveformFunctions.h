#ifndef WAVEFORMFUNCTIONS_H
#define WAVEFORMFUNCTIONS_H

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"

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
          bool passesHitSelection(art::Ptr< recob::Hit > hit, double HIT_GOODNESSOFFIT_CUT, 
              int HIT_MULTIPILCITY); 

          double convertXToTicks(double xPosition, int WAVEFORM_DRIFT_START_TICK, int WAVEFORM_DRIFT_SIZE, double X_WIDTH);

          double convertTicksToX(int tick, int WAVEFORM_DRIFT_START_TICK, int WAVEFORM_DRIFT_SIZE, double X_WIDTH);

          TH1D* applyGlobalBaselineCorrection(TH1D *h_rawD, TH1D *h_rawDCorrected);

          double findXCorrection(TH1D *summedWaveform, TH1D *h, int NUMBER_TICKS_PER_BIN, double mean);

          std::vector<double> getSigma(TH1D *h_rawDCorrected);

          double getMedian(TH1D *h);

          double getRms2(TH1D *h);

    };

}

#endif
