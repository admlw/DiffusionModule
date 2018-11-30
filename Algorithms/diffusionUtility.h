#ifndef DIFFUSIONUTILITY_H
#define DIFFUSIONUTILITY_H

// larsoft includes
#include "lardataobj/RecoBase/Track.h"

// cc includes
#include <numeric>
#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>

namespace diffUtil{

  class diffusionUtility{

    public:

      std::vector<double> getMomentumVector(const recob::Track& a);

      float getAngle(const std::vector<double> a, const std::vector<double> b, diffusionUtility diffutil, std::string proj);

      std::vector<float> getUnitVector(std::vector<double> a); 

      float getDotProduct(std::vector<float> a, std::vector<float> b); 
  };

}

#endif
