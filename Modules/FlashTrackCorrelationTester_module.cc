////////////////////////////////////////////////////////////////////////
// Class:       FlashTrackCorrelationTester
// Plugin Type: analyzer (art v2_05_00)
// File:        FlashTrackCorrelationTester_module.cc
//
// Generated at Mon Oct 23 14:00:49 2017 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
//
// Intended to be used with a single muon sample, else this won't make sense
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "TTree.h"

class FlashTrackCorrelationTester;

class FlashTrackCorrelationTester : public art::EDAnalyzer {
public:
  explicit FlashTrackCorrelationTester(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FlashTrackCorrelationTester(FlashTrackCorrelationTester const &) = delete;
  FlashTrackCorrelationTester(FlashTrackCorrelationTester &&) = delete;
  FlashTrackCorrelationTester & operator = (FlashTrackCorrelationTester const &) = delete;
  FlashTrackCorrelationTester & operator = (FlashTrackCorrelationTester &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  art::ServiceHandle< art::TFileService > tfs;
  TTree* flashChecker;

  double trackCenterY = -9999;
  double trackCenterZ = -9999;
  double trackCenterX = -9999;
  double trackWidthY = -9999;
  double trackWidthZ = -9999;
  double trackWidthX = -9999;

  double flashCenterY = -9999;
  double flashCenterZ = -9999;
  double flashWidthY = -9999;
  double flashWidthZ = -9999;

  double trMinusFlWidthZ = -9999;
  double trMinusFlWidthY = -9999;
  double trMinusFlCenterZ = -9999;
  double trMinusFlCenterY = -9999;

};


FlashTrackCorrelationTester::FlashTrackCorrelationTester(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
 
{}

void FlashTrackCorrelationTester::analyze(art::Event const & e)
{
  // Implementation of required member function here.

  art::Handle< std::vector<recob::Track> > trackHandle;
  e.getByLabel("pandoraCosmic", trackHandle);

  art::Handle< std::vector<recob::OpFlash> > opFlashHandle;
  e.getByLabel("simpleFlashBeam", opFlashHandle);

  int trackCounter = 0;
  for (auto const& track : (*trackHandle)){

    if (track.Length() < 30) continue;

    auto trackStart = track.Vertex();
    double trackStartY = trackStart.Y();
    double trackStartZ = trackStart.Z();
    double trackStartX = trackStart.X();

    auto trackEnd = track.End();
    double trackEndY = trackEnd.Y();
    double trackEndZ = trackEnd.Z();
    double trackEndX = trackEnd.X();

    // construct yCenter, yWidth, zCenter and zWidth for track

    trackWidthY = std::abs(trackStartY-trackEndY)/2.0;
    trackCenterY = std::min(trackStartY, trackEndY) + trackWidthY;

    trackWidthZ = std::abs(trackStartZ-trackEndZ)/2.0;
    trackCenterZ = std::min(trackStartZ, trackEndZ) + trackWidthZ;

    trackWidthX = std::abs(trackStartX-trackEndX)/2.0;
    trackCenterX = std::min(trackStartX, trackEndX) + trackWidthX;

    trackCounter++;
  }

  std::cout << "NUMBER OF TRACKS: " << trackCounter << std::endl;

  if (trackCounter != 1) return;

  std::vector<recob::OpFlash> opFlashVector;

  double largestFlashChecker = 0;
  for (auto const& opFlash : (*opFlashHandle)){

    double opFlashPeSum = 0;
    for (int i = 0; i < 32; i++){

      opFlashPeSum = opFlashPeSum + opFlash.PE(i);

    }

     // just fill vector with opFlashes. Ensure largest flash is first in vector;
     if (opFlashPeSum > largestFlashChecker){
       opFlashVector.insert(opFlashVector.begin(), opFlash);
       largestFlashChecker = opFlashPeSum;  
     }
     else opFlashVector.push_back(opFlash);

  }

  std::cout << "NUMBER OF FLASHES: " << opFlashVector.size() << std::endl;

  // get largest opFlash
  if (opFlashVector.size() == 0) return;
  auto const& largestOpFlash = opFlashVector.at(0);

  flashWidthY = largestOpFlash.YWidth();
  flashCenterY = largestOpFlash.YCenter();
  flashWidthZ = largestOpFlash.ZWidth();
  flashCenterZ = largestOpFlash.ZCenter();

  trMinusFlWidthZ = std::abs(trackWidthZ - flashWidthZ);
  trMinusFlWidthY = std::abs(trackWidthY - flashWidthY);
  trMinusFlCenterZ = std::abs(trackCenterZ - flashCenterZ);
  trMinusFlCenterY = std::abs(trackCenterY - flashCenterY);

  std::cout << ">> PRINTING TRACK MINUS FLASH INFORMATION"
    << "\n-->> WidthZ: " << trMinusFlWidthZ
    << "\n-->> WidthY: " << trMinusFlWidthY
    << "\n-->> CenterZ: " << trMinusFlCenterZ
    << "\n-->> CenterY: " << trMinusFlCenterY
    << std::endl;

  flashChecker->Fill();

}

void FlashTrackCorrelationTester::beginJob()
{
  // Implementation of optional member function here.

  flashChecker = tfs->make<TTree>("flashChecker", "flashChecker");

  flashChecker->Branch("flashWidthY", &flashWidthY, "flashWidthY/D");
  flashChecker->Branch("flashWidthZ", &flashWidthZ, "flashWidthZ/D");
  flashChecker->Branch("flashCenterY", &flashCenterY, "flashCenterY/D");
  flashChecker->Branch("flashCenterZ", &flashCenterZ, "flashCenterZ/D");
  flashChecker->Branch("trackWidthY", &trackWidthY, "trackWidthY/D");
  flashChecker->Branch("trackWidthZ", &trackWidthZ, "trackWidthZ/D");
  flashChecker->Branch("trackWidthX", &trackWidthX, "trackWidthX/D");
  flashChecker->Branch("trackCenterY", &trackCenterY, "trackCenterY/D");
  flashChecker->Branch("trackCenterZ", &trackCenterZ, "trackCenterZ/D");
  flashChecker->Branch("trackCenterX", &trackCenterX, "trackCenterX/D");
  flashChecker->Branch("trMinusFlCenterY", &trMinusFlCenterY, "trMinusFlCenterY/D");
  flashChecker->Branch("trMinusFlCenterZ", &trMinusFlCenterZ, "trMinusFlCenterZ/D");
  flashChecker->Branch("trMinusFlWidthY", &trMinusFlWidthY, "trMinusFlWidthY/D");
  flashChecker->Branch("trMinusFlWidthZ", &trMinusFlWidthZ, "trMinusFlWidthZ/D");

}

void FlashTrackCorrelationTester::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(FlashTrackCorrelationTester)
