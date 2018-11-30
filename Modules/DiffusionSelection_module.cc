////////////////////////////////////////////////////////////////////////
// Class:       DiffusionSelection
// Plugin Type: analyzer (art v2_05_00)
// File:        DiffusionSelection_module.cc
//
// Generated at Fri Oct 13 16:46:35 2017 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

// base includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// additional includes
#include "lardataobj/RecoBase/Track.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpFlash.h"

// ROOT includes
#include "TH2.h"

// local includes
#include "ubana/DiffusionSelection/Algos/diffusionUtility.h"

class DiffusionSelection;


class DiffusionSelection : public art::EDFilter {
  public:
    explicit DiffusionSelection(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    DiffusionSelection(DiffusionSelection const &) = delete;
    DiffusionSelection(DiffusionSelection &&) = delete;
    DiffusionSelection & operator = (DiffusionSelection const &) = delete;
    DiffusionSelection & operator = (DiffusionSelection &&) = delete;

    // Required functions.
    bool filter(art::Event & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

  private:

    bool isTrueT0 = true;

    art::ServiceHandle< art::TFileService > tfs;

    diffUtil::diffusionUtility _diffutilInstance;

    std::string fTrackLabel;
    std::string fTrackTruthMatcherLabel;
    double fTrackLengthCut;
    double fTrackAngleCutXZLow;
    double fTrackAngleCutXZHigh;
    double fTrackAngleCutYZLow;
    double fTrackAngleCutYZHigh;

    // track position stuff
    double recoY;
    double recoZ;

    // flash stuff
    std::vector<recob::OpFlash> opFlashVec;

    // track angle stuff
    std::vector<double> trackMomentum;
    std::vector<double> zDirection;

    double trackAngleXZ;
    double trackAngleYZ;

    double recoT0;

    // truth information
    double trueTrackT0;
    int trueTrackID;
    int trueTriggerType;

    // other stuff
    int channel;
    int startTick;
    int endTick;
    float peakTime;
    float sigmaPeakTime;
    float rms;
    float peakAmplitude;
    float sigmaPeakAmplitude;
    float summedAdc;
    float hitIntegral;
    float hitSigmaIntegral;
    short int multiplicity;
    short int local_index;
    float goodnessOfFit;
    int dof;

    TH2D* trackAngles;
    TH2D* trackAnglesT0Tagged;
};


DiffusionSelection::DiffusionSelection(fhicl::ParameterSet const & p)
  // :
  // More initializers here.
{

  fTrackLabel = p.get<std::string> ("TrackLabel");
  fTrackTruthMatcherLabel = p.get<std::string> ("TrackTruthMatcherLabel");
  fTrackLengthCut = p.get<double> ("TrackLengthCut");
  fTrackAngleCutXZLow  = p.get<double> ("TrackAngleCutXZLow");
  fTrackAngleCutXZHigh = p.get<double> ("TrackAngleCutXZHigh");
  fTrackAngleCutYZLow  = p.get<double> ("TrackAngleCutYZLow");
  fTrackAngleCutYZHigh = p.get<double> ("TrackAngleCutYZHigh");

  produces< std::vector<recob::Track> >();
  produces< std::vector<anab::T0> >();
  produces< art::Assns<recob::Track, anab::T0> >();
  produces< std::vector<recob::Hit> >();
  produces< art::Assns<recob::Track, recob::Hit> >();
}

void DiffusionSelection::beginJob()
{

  trackAngles = tfs->make<TH2D>("trackAngles", "", 36, 0, 180, 36, 0, 180);
  trackAnglesT0Tagged = tfs->make<TH2D>("trackAnglesT0Tagged", "", 36, 0, 180, 36, 0, 180); 

}

bool DiffusionSelection::filter(art::Event & e)
{
  bool isPass = false;

  art::Handle< std::vector<recob::Track> > trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);

  //  art::Handle< std::vector<recob::OpFlash> > opFlashHandle;
  //  e.getByLabel("simpleFlashCosmic", opFlashHandle);

  //art::FindManyP<simb::MCParticle, anab::BackTrackerMatchingData> mcpsFromTracks(trackHandle, e, fTrackTruthMatcherLabel);
  art::FindManyP<recob::Hit> hitsFromTracks(trackHandle, e, fTrackLabel);
  art::FindManyP<anab::T0>   t0FromTracks(trackHandle, e, "pandoraCosmicT0RecoLoose");

  // produces a new trackCollection for passing tracks result is that we have
  // only events in which at least one track passes the selection, and only 
  // those tracks are reconstructed
  std::unique_ptr< std::vector<recob::Track> > trackCollection( new std::vector<recob::Track>);
  std::unique_ptr< std::vector<anab::T0> > t0Collection( new std::vector<anab::T0>);
  std::unique_ptr< art::Assns<recob::Track, anab::T0> > trackT0Assn( new art::Assns<recob::Track, anab::T0>);
  std::unique_ptr< std::vector<recob::Hit> > hitCollection( new std::vector<recob::Hit> );
  std::unique_ptr< art::Assns<recob::Track, recob::Hit> > trackHitAssn( new art::Assns<recob::Track, recob::Hit>);

  for (auto const& thisTrack : (*trackHandle)){

    if (thisTrack.Length() < fTrackLengthCut) continue;
        
    std::vector< art::Ptr<anab::T0> > t0s = t0FromTracks.at(thisTrack.ID());
    
    trackMomentum = _diffutilInstance.getMomentumVector(thisTrack);
    zDirection = {0,0,1};

    trackAngleXZ = _diffutilInstance.getAngle(trackMomentum, zDirection, _diffutilInstance, "xz");
    trackAngleYZ = _diffutilInstance.getAngle(trackMomentum, zDirection, _diffutilInstance, "yz");

    trackAngles->Fill(trackAngleXZ, trackAngleYZ);
    
    if (t0s.size() == 1){

      trackAnglesT0Tagged->Fill(trackAngleXZ, trackAngleYZ);

    }

    if ((trackAngleXZ <= fTrackAngleCutXZHigh && trackAngleXZ >= fTrackAngleCutXZLow
          && trackAngleYZ <= fTrackAngleCutYZHigh && trackAngleYZ >= fTrackAngleCutYZLow) 
        || (trackAngleXZ >= (180 - fTrackAngleCutXZHigh) && (trackAngleXZ <= (180 - fTrackAngleCutXZLow)) 
          && trackAngleYZ >= (180 - fTrackAngleCutYZHigh) && trackAngleYZ <= (180 - fTrackAngleCutYZLow))){

      trackCollection->push_back(recob::Track(thisTrack.Trajectory(), 
            thisTrack.ParticleId(),
            thisTrack.Chi2(),
            thisTrack.Ndof(),
            thisTrack.VertexCovarianceLocal5D(),
            thisTrack.EndCovarianceLocal5D(),
            thisTrack.ID()));

      //std::vector< art::Ptr<simb::MCParticle> > mcps = mcpsFromTracks.at(thisTrack.ID());

      std::cout << ">> RECO TRACK INFORMATION " 
        << "\n >> Track Start x: " << thisTrack.Vertex().X() << " End x  : " << thisTrack.End().X() 
        << "\n >> Track Start y: " << thisTrack.Vertex().Y() << " End y  : " << thisTrack.End().Y()
        << "\n >> Track Start z: " << thisTrack.Vertex().Z() << " End z  : " << thisTrack.End().Z()
        << std::endl;

      recoY = std::min(thisTrack.Vertex().Y(), thisTrack.End().Y()) + (std::abs(thisTrack.Vertex().Y() - thisTrack.End().Y())/2.0);
      recoZ = std::min(thisTrack.Vertex().Z(), thisTrack.End().Z()) + (std::abs(thisTrack.Vertex().Z() - thisTrack.End().Z())/2.0);


      if (t0s.size() != 1) {

        std::cout << "nT0: " << t0s.size() << std::endl; 

        continue;
      }
      else {

        t0Collection->push_back(*(t0s.at(0).get()));

        std::cout << ">> T0 INFORMATION "
          << "\n >> T0 tagged time: " << t0s.at(0)->Time() << std::endl;

      }



      //      if (isTrueT0 == true){
      //        for (auto const& thisMcp : mcps){

      //          trueTrackT0 = thisMcp->T();

      //          std::cout << ">> TRUE TRACK INFORMATION "
      //            << "\n >> Track Start T: " << trueTrackT0 << " EndT: " << thisMcp->EndT()
      //            << "\n >> track Start X: " << thisMcp->Vx() << "EndX: "<< thisMcp->EndX() << std::endl;

      //trueTrackID = thisMcp->TrackId();
      //trueTriggerType = 2;
      //t0Collection->push_back(anab::T0(trueTrackT0,
      //      trueTriggerType,
      //      -1,
      //      trueTrackID,
      //      -1
      //      ));

      //        }
      //      }
      //      else {
      /*        for (auto const& thisFlash : (*opFlashHandle)){

                double yDistance = std::abs(thisFlash.YCenter() - recoY);
                double zDistance = std::abs(thisFlash.ZCenter() - recoZ);

                if (yDistance > 10.0 || zDistance > 10.0) continue;

                std::cout << "Y Distance: " << yDistance << std::endl;
                std::cout << "Z Distance: " << zDistance << std::endl;

                opFlashVec.push_back(thisFlash);

                }

                if (opFlashVec.size() !=1){
                std::cout << "opFlashVecSize: " << opFlashVec.size() << std::endl;
                recoT0 = -9999;
                }
                else
                recoT0 = opFlashVec.at(0).Time();

                t0Collection->push_back(anab::T0(recoT0, 2, -1, -1, -1));
                std::cout << ">> RECO FLASH INFORMATION " 
                << "\n >> recoT0: " << recoT0 << std::endl;

      //      }
      */
      std::vector< art::Ptr<recob::Hit> > hits = hitsFromTracks.at(thisTrack.ID());

      for (auto const& thisHit : hits){

        channel = thisHit->Channel();     
        startTick = thisHit->StartTick();
        endTick = thisHit->EndTick();
        peakTime = thisHit->PeakTime();
        sigmaPeakTime = thisHit->SigmaPeakTime();
        rms = thisHit->RMS();
        peakAmplitude = thisHit->PeakAmplitude();
        sigmaPeakAmplitude = thisHit->SigmaPeakAmplitude();
        summedAdc = thisHit->SummedADC();
        hitIntegral = thisHit->Integral();
        hitSigmaIntegral = thisHit->SigmaIntegral();
        multiplicity = thisHit->Multiplicity();
        local_index = thisHit->LocalIndex();
        goodnessOfFit = thisHit->GoodnessOfFit();
        dof = thisHit->DegreesOfFreedom();
        auto view = thisHit->View();
        auto signalType = thisHit->SignalType();
        auto wireID = thisHit->WireID();

        hitCollection->push_back(recob::Hit(channel, 
              startTick,
              endTick,
              peakTime,
              sigmaPeakTime,
              rms,
              peakAmplitude,
              sigmaPeakAmplitude,
              summedAdc,
              hitIntegral,
              hitSigmaIntegral,
              multiplicity,
              local_index,
              goodnessOfFit,
              dof,
              view,
              signalType,
              wireID));

      }

      std::cout << "!!Hit collection size: " << hitCollection->size() << std::endl;

      art::Ptr<recob::Track> trkPtr(trackHandle, (size_t)thisTrack.ID());

      util::CreateAssn(*this, e, *trackCollection, *hitCollection, *trackHitAssn, 0, hitCollection->size());
      util::CreateAssn(*this, e, *trackCollection, *t0Collection, *trackT0Assn, 0, t0Collection->size());

      isPass = true;
      break;

    }

  }

  e.put(std::move(trackCollection));
  e.put(std::move(t0Collection));
  e.put(std::move(hitCollection));
  e.put(std::move(trackT0Assn));
  e.put(std::move(trackHitAssn));
  return isPass;

}

void DiffusionSelection::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(DiffusionSelection)
