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

      /// does hit pass selection?
      bool passesHitSelection(art::Ptr< recob::Hit > hit, double hitGOFCut, int hitMultiplicity); 

      /// convert x position to number of ticks
      double convertXToTicks(double xPosition, int wvfmDriftStartTick, int wvfmDriftSize, double xWidth);

      /// convert number of ticks to an x position
      double convertTicksToX(int tick, int wvfmDriftStartTick, int wvfmDriftSize, double xWidth);

      /// apply a global baseline correction to the waveforms
      TH1D* applyGlobalBaselineCorrection(TH1D *h_rawD, TH1D *h_rawDCorrected);

      /// find x correction factor
      double findXCorrection(TH1D *summedWaveform, TH1D *h, int nTicksPerBin, double mean);

      /// get sigma value
      std::vector<double> getSigma(TH1D *h_rawDCorrected);

      /// get median value
      double getMedian(TH1D *h);

      /// get rms2 value
      double getRms2(TH1D *h);

  };

}

#endif
