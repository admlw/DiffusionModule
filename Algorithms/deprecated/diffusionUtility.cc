#include "diffusionUtility.h"

namespace diffUtil{

  std::vector<double> diffusionUtility::getMomentumVector(const recob::Track& a){

    std::vector<double> aMom;
    auto const& smv = a.StartMomentumVector();
    aMom = {smv.X(), smv.Y(), smv.Z()};

    return aMom;

  }

  float diffusionUtility::getAngle(const std::vector<double> a, const std::vector<double> b, diffusionUtility diffutil, std::string proj){

    std::vector<double> aMom;
    std::vector<double> bMom;

    if (proj == "no"){

      aMom = {a.at(0), a.at(1), a.at(2)};
      bMom = {b.at(0), b.at(1), b.at(2)};

    }

    else if ((proj == "xz") | (proj == "zx")){

      aMom = {a.at(0), a.at(2)};
      bMom = {b.at(0), b.at(2)};

    }

    else if ((proj == "xy") | (proj == "yx")){

      aMom = {a.at(0), a.at(1)};
      bMom = {b.at(0), b.at(1)};

    }

    else if ((proj == "yz") | (proj == "zy")){

      aMom = {a.at(1), a.at(2)};
      bMom = {b.at(1), b.at(2)};

    }
    else
      throw std::invalid_argument("Valid arguments are 'no', 'xy', 'xz' and 'yz'");


    std::vector<float> aMomUnit = diffutil.getUnitVector(aMom);
    std::vector<float> bMomUnit = diffutil.getUnitVector(bMom);

    float angle = std::acos(diffutil.getDotProduct(aMomUnit, bMomUnit)) * 180 / 3.14159;

    return angle;
  }


  std::vector<float> diffusionUtility::getUnitVector(std::vector<double> a){
    /*
     *     float ax = a.at(0);
     *         float ay = a.at(1);
     *             float az = a.at(2);
     *                 float aMag = std::sqrt(ax*ax + ay*ay + az*az);
     *                     float axNorm = ax / aMag;
     *                         float ayNorm = ay / aMag;
     *                             float azNorm = az / aMag;
     *                             */

    std::vector<float> aNorm;

    double magnitude = 0;
    for (size_t i = 0; i < a.size(); i++){

      magnitude += (a.at(i) * a.at(i));

    }

    magnitude = std::sqrt(magnitude);

    for (size_t i = 0; i < a.size(); i++){

      aNorm.push_back(a.at(i)/magnitude);

    }

    return aNorm;

  }

  float diffusionUtility::getDotProduct(std::vector<float> a, std::vector<float> b){

    if (a.size() != b.size())
      throw std::invalid_argument("Cannot dot product vectors of different sizes");

    float dotProduct = 0;
    for (size_t i = 0; i < a.size(); i++){

      dotProduct += a.at(i)*b.at(i);

    }

    return dotProduct;

  }


}
