#include "diffusionFunctions.h"

/*
namespace diffmod{
    std::vector<double> diffusionFunctions::getMomentumVector(const recob::Track& a){

      std::vector<double> aMom;
      auto const& smv = a.StartMomentumVector();
      aMom = {smv.X(), smv.Y(), smv.Z()};

      return aMom;

    }
}

float diffusionFunctions::getAngle(const std::vector<double> a, const std::vector<double> b, diffusionFunctions diffutil, std::string proj){

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

  if (angle > 90)
    return 180-angle;
  else return angle;
}
*/

/*
std::vector<float> diffusionFunctions::getUnitVector(std::vector<double> a){
   //
   //      *     float ax = a.at(0);
   //           *         float ay = a.at(1);
   //                *             float az = a.at(2);
   //                     *                 float aMag = std::sqrt(ax*ax + ay*ay + az*az);
   //                          *                     float axNorm = ax / aMag;
   //                               *                         float ayNorm = ay / aMag;
   //                                    *                             float azNorm = az / aMag;
   //                                         *                             

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

float diffusionFunctions::getDotProduct(std::vector<float> a, std::vector<float> b){

  if (a.size() != b.size())
    throw std::invalid_argument("Cannot dot product vectors of different sizes");

  float dotProduct = 0;
  for (size_t i = 0; i < a.size(); i++){

    dotProduct += a.at(i)*b.at(i);

  }

  return dotProduct;

}
*/

/*
bool diffusionFunctions::passesPreSelection(recob::Track const& track, diffusionFunctions _diffusionFunctions_instance, double trackAngleXZLowBound, double trackAngleXZHighBound, double trackAngleYZLowBound, double trackAngleYZHighBound){

  double trackLength = track.Length();

  double trackBeginningEndAngle = findTrackStraightness(track, _diffusionFunctions_instance);
 
  std::vector<double> trackMomentum = _diffusionFunctions_instance.getMomentumVector(track);
  std::vector<double> zDir = {0,0,1};

  float trackAngleXZ = _diffusionFunctions_instance.getAngle(trackMomentum, zDir, _diffusionFunctions_instance, "xz");
  float trackAngleYZ = _diffusionFunctions_instance.getAngle(trackMomentum, zDir, _diffusionFunctions_instance, "yz");

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

bool diffusionFunctions::passesHitSelection(recob::Hit const* hit, double HIT_GOODNESSOFFIT_CUT){

  if (hit->Multiplicity() == 1 && hit->GoodnessOfFit() < HIT_GOODNESSOFFIT_CUT && hit->View() ==2 && hit->Channel() > 6150) {
    //std::cout << "Hit passed" << std::endl;
    return true;
  }
  else {
    if (hit->Multiplicity() != 1) {
        //std::cout << "Hit multiplicity is " << hit->Multiplicity() << std::endl;
    }
    if (hit->GoodnessOfFit() < HIT_GOODNESSOFFIT_CUT) {
        //std::cout << "Hit goodness-of-fit " << hit->GoodnessOfFit() << std::endl;
    }
    if (hit->View() != 2) {
        //std::cout << "Hit view " << hit->View() << std::endl;
    }
    if (hit->Channel() < 6150) {
        //std::cout << "Hit in U-shorted region" << std::endl;
    }
    return false;
  }
}

double diffusionFunctions::findTrackStraightness(recob::Track const& track, diffusionFunctions _diffusionFunctions_instance){


  using Vector_t = recob::tracking::Vector_t;
  //using recob::Track::Vector_t = tracking::Vector_t;
  using TrajectoryPoint_t = recob::Track::TrajectoryPoint_t;

  size_t firstPoint = track.FirstPoint();

  // get reference vector
  TrajectoryPoint_t firstTrajPoint = track.TrajectoryPoint(firstPoint);
  auto firstPointMom = firstTrajPoint.momentum;

  std::vector<double> firstVec = _diffusionFunctions_instance.getUnitVector(firstPointMom); 

  size_t lastPoint = track.LastPoint();
  TrajectoryPoint_t lastTrajPoint = track.TrajectoryPoint(lastPoint); 
  auto lastPointMom = lastTrajPoint.momentum;

  std::vector<double> lastVec = _diffusionFunctions_instance.getUnitVector(lastPointMom);

  double dotProd = _diffusionFunctions_instance.getDotProduct(firstVec, lastVec);

  double angle = std::acos(1-dotProd) * 180/3.14159;

  return angle;
}

double diffusionFunctions::getDotProduct(std::vector<double> a, std::vector<double> b){

  double dotProduct = 0.0;
  for (int i = 0; i < a.size(); i++){

    dotProduct += a.at(i) * b.at(i);

  }

  return dotProduct;

}
*/

/*
std::vector<double> diffusionFunctions::getUnitVector(recob::tracking::Vector_t vec3){

  std::vector<double> returnVec;


  double x = vec3.x();
  double y = vec3.y();
  double z = vec3.z();


  double vMag = std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2));

  double xNorm = x/vMag;
  double yNorm = y/vMag;
  double zNorm = z/vMag;

  returnVec.push_back(xNorm);
  returnVec.push_back(yNorm);
  returnVec.push_back(yNorm);


  return returnVec;
}

double diffusionFunctions::convertXToTicks(double xPosition, int WAVEFORM_DRIFT_START_TICK, int WAVEFORM_DRIFT_SIZE, double X_WIDTH){

  //    double tick = WAVEFORM_START_TICK + 17.9427*xPosition;
  double tick = WAVEFORM_DRIFT_START_TICK + (((double)WAVEFORM_DRIFT_SIZE/X_WIDTH)*xPosition);
  return tick;

}

TH1D* diffusionFunctions::applyGlobalBaselineCorrection(TH1D* h_rawD, TH1D* h_rawDCorrected){

  // first loop over the bins and find the beginning and end of the ROI
  std::pair<int, int> firstLastBinROI = std::make_pair(-1, -1);
  std::pair<int, int> clippedSignalBins = std::make_pair(-1,-1);
  double firstBinContent = h_rawD->GetBinContent(1);
  double lastBinContent = h_rawD->GetBinContent(h_rawD->GetNbinsX()); 
  double binContent;

  for (int i = 1; i <= h_rawD->GetNbinsX(); i++) {

    binContent = h_rawD->GetBinContent(i);
    if (binContent != firstBinContent){

      firstLastBinROI.first = i+5;
      break;
    }

  }

  for (int i = h_rawD->GetNbinsX(); i >=1; i--){ 

    binContent = h_rawD->GetBinContent(i);
    if (binContent != lastBinContent) {

      firstLastBinROI.second = i-5;
      break;
    }
  }

  // clip out signal region;
  clippedSignalBins.first = h_rawD->GetMaximumBin()-20;
  clippedSignalBins.second = h_rawD->GetMaximumBin()+20;

  // loop over bins and find baseline within ROI (outside of signal)
  double cumulativeSum = 0;
  int tickCounter = 0;
  for (int i = firstLastBinROI.first; i <= firstLastBinROI.second; i++){

    if (i >= clippedSignalBins.first && i <= clippedSignalBins.second) continue;

    cumulativeSum = cumulativeSum + h_rawD->GetBinContent(i);
    tickCounter++;

  }

  double baseline = cumulativeSum/tickCounter;


  // correct baseline and return histogram
  for (int i = 0; i <= h_rawDCorrected->GetNbinsX(); i++){ 
    h_rawDCorrected->SetBinContent(i, h_rawD->GetBinContent(i) - baseline);
  }

  return h_rawDCorrected;

}
*/

/*
TH1D* diffusionFunctions::rebin(TH1D* h_rawD, TH1D* h_rebinned, int NUMBER_BINS_PER_BIN, int WAVEFORM_DRIFT_SIZE, int NUMBER_DRIFT_BINS){

  //
  // rebin function
  //

  // no standard way to split single bin in to multiple bins in root
  // used to change every 1-tick bin in to a 1/NUMBER_BINS_PER_BIN-tick bin
  // note that the ADC value is not scaled by 1/NUMBER_BINS_PER_BIN

  h_rebinned->SetBins(WAVEFORM_DRIFT_SIZE/NUMBER_DRIFT_BINS * NUMBER_BINS_PER_BIN, h_rawD->GetXaxis()->GetXmin(), h_rawD->GetXaxis()->GetXmax());

  for (int i=1; i <= h_rawD->GetNbinsX(); i++){

    for (int j = 1; j <= NUMBER_BINS_PER_BIN; j++){

      h_rebinned->SetBinContent(((i-1)*NUMBER_BINS_PER_BIN)+j, h_rawD->GetBinContent(i));

    }

  }

  return h_rebinned;

}

double diffusionFunctions::findXCorrection(diffusionFunctions _diffusionFunctions_instance, TH1D* summedWaveform, TH1D* h, int NUMBER_TICKS_PER_BIN, double mean){

  int centerBin = -1;
  if (summedWaveform->GetMaximum() == 0){
    centerBin = NUMBER_TICKS_PER_BIN/2;
  }
  else {

    centerBin = summedWaveform->FindBin(_diffusionFunctions_instance.getSigma(summedWaveform).at(0));

  }


  int distanceToBinCenter = 0;
  int testXCorr;
  double rms2 = 1000;
  for (int i = -5; i <= 5; i++){
    TH1D* h_summedClone = (TH1D*)summedWaveform->Clone("h_summedClone");
    TH1D* h_clone = (TH1D*)h->Clone("h_clone");


    testXCorr = centerBin - h->FindBin(mean) + i;
    for (int ntick = 1; ntick <= h->GetNbinsX(); ntick++){

      h_clone->SetBinContent(ntick, h->GetBinContent(ntick + testXCorr));

    }

    h_summedClone->Add(h_clone);


    double rms2Test = getRms2(h_summedClone);

    if (rms2Test < rms2){

      rms2 = rms2Test;
      distanceToBinCenter = testXCorr;

    }

    h_clone->Delete();
    h_summedClone->Delete();

  }

  return distanceToBinCenter;
}
*/

/*
std::vector<double> diffusionFunctions::getSigma(TH1D* h_rawDCorrected){

  std::vector<double> returnVector;
  if (h_rawDCorrected->Integral() == 0 || h_rawDCorrected->GetMaximum() <= 0 ) {
    returnVector = {-1.0, -1.0, -1.0};
    return returnVector;
  }

  // loop over histogram, find first point below 0.1 ADC, fit only between these regions
  int maxBin = h_rawDCorrected->GetMaximumBin();
  double maxVal = h_rawDCorrected->GetMaximum();
  double fitLowerLimit = 0;
  double fitHigherLimit = 0;
  double percentageOfHeight = 2;
  double cutOff = maxVal* percentageOfHeight/100;

  for (int i = 0; i < 200; i++){

    double binValueLow = h_rawDCorrected->GetBinContent(maxBin-i);
    if (binValueLow < cutOff){

      fitLowerLimit = h_rawDCorrected->GetBinCenter(maxBin-i); 
      break;
    }

  }

  for (int i = 0; i < 200; i++){

    double binValueHigh = h_rawDCorrected->GetBinContent(maxBin+i);
    if (binValueHigh < cutOff){

      fitHigherLimit = h_rawDCorrected->GetBinCenter(maxBin+i); 
      break;
    }

  }


  h_rawDCorrected->Fit("gaus", "q", "", fitLowerLimit, fitHigherLimit);

  double sigma;
  double chisq;
  double mean;
  double mean_err;

  if (h_rawDCorrected->GetFunction("gaus")){

    mean = h_rawDCorrected->GetFunction("gaus")->GetParameter(1);
    mean_err = h_rawDCorrected->GetFunction("gaus")->GetParError(1);
    sigma = h_rawDCorrected->GetFunction("gaus")->GetParameter(2);
    chisq = h_rawDCorrected->GetFunction("gaus")->GetChisquare();

  }
  else {
    sigma = -1;
    chisq = -1;
  }

  returnVector = {mean, sigma, chisq};

  return returnVector;
}

double diffusionFunctions::convertTicksToX(int tick, int WAVEFORM_DRIFT_START_TICK, int WAVEFORM_DRIFT_SIZE, double X_WIDTH){

  //double xPos = (tick - 812)/17.9427;
  double xPos = (tick - WAVEFORM_DRIFT_START_TICK)/(WAVEFORM_DRIFT_SIZE/X_WIDTH);
  return xPos;

}

double diffusionFunctions::getMedian(TH1D* h){

  double quantile = 0.5;
  double median; 

  h->GetQuantiles(1,&median,&quantile);                                

  return median;

}
*/

/*
double diffusionFunctions::getRms2(TH1D* h){

  double mean = 0;
  double mean2 = 0;
  double totalCharge = 0;
  double rms2 = 0;
  double threshold = 10.0 * h->GetMaximum()/100.0;

  for (int j = 1; j < h->GetNbinsX(); j++) {

    h->SetBinContent(j, h->GetBinContent(j) - threshold);

    if (h->GetBinContent(j) >= 0){

      totalCharge = totalCharge + h->GetBinContent(j);
      mean = mean + (double)h->GetBinCenter(j) * 0.5 * (double)h->GetBinContent(j);
      mean2 = mean2 + (double)h->GetBinCenter(j) * 0.5 * (double)h->GetBinCenter(j) *0.5 * (double)h->GetBinContent(j);
    }

    h->SetBinContent(j, h->GetBinContent(j) + threshold);

  }

  if (totalCharge !=0){

    mean = mean/totalCharge;
    mean2 = mean2/totalCharge;
    rms2 = sqrt(mean2 - mean*mean) * sqrt(mean2-mean*mean);

  }

  return rms2;

}

*/

