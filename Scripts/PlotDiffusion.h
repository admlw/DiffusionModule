#include <iostream>
#include <utility>

#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TVector.h"
#include "TPad.h"
#include "TGraphErrors.h"

// cpp includes
#include <numeric>
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>

void increaseError(TH1D *h);
std::pair<int, int> getBinXError(TH1D *driftHisto);

namespace diffmod {

    class WaveformFunctions {
        public: 
          double convertXToTicks(double xPosition);

          double convertTicksToX(int tickVal);

    };

}


