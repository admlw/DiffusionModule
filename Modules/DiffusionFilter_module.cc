////////////////////////////////////////////////////////////////////////
// Class:       DiffusionFilter
// Plugin Type: analyzer (art v2_05_00)
// File:        DiffusionFilter_module.cc
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
#include "TTree.h"
#include "TMath.h"

// local includes
#include "ubana/DiffusionModule/Algorithms/diffusionUtility.h"

class DiffusionFilter;


class DiffusionFilter : public art::EDFilter {
  public:
    explicit DiffusionFilter(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    DiffusionFilter(DiffusionFilter const &) = delete;
    DiffusionFilter(DiffusionFilter &&) = delete;
    DiffusionFilter & operator = (DiffusionFilter const &) = delete;
    DiffusionFilter & operator = (DiffusionFilter &&) = delete;

    // Required functions.
    bool filter(art::Event & e) override;

    // Selected optional functions.
    void beginJob() override;

  private:

    bool isTrueT0 = true;

    art::ServiceHandle< art::TFileService > tfs;

    diffUtil::diffusionUtility _diffutilInstance;

    TTree* tree;
    std::string fTrackLabel;
    std::string fTrackTruthMatcherLabel;
    std::string fT0Label;
    double fTrackLengthCut;
    double fTrackAngleCutXZLow;
    double fTrackAngleCutXZHigh;
    double fTrackAngleCutYZLow;
    double fTrackAngleCutYZHigh;

    // flash stuff
    std::vector<recob::OpFlash> opFlashVec;

    double track_length;
    double track_angle_xz;
    double track_angle_yz;
    bool track_has_t0;

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

};


DiffusionFilter::DiffusionFilter(fhicl::ParameterSet const & p)
  // :
  // More initializers here.
{

  fTrackLabel = p.get<std::string> ("TrackLabel");
  fTrackTruthMatcherLabel = p.get<std::string> ("TrackTruthMatcherLabel");
  fT0Label = p.get<std::string> ("T0Label");
  fTrackLengthCut = p.get<double> ("TrackLengthCut");
  fTrackAngleCutXZLow  = p.get<double> ("TrackAngleCutXZLow");
  fTrackAngleCutXZHigh = p.get<double> ("TrackAngleCutXZHigh");
  fTrackAngleCutYZLow  = p.get<double> ("TrackAngleCutYZLow");
  fTrackAngleCutYZHigh = p.get<double> ("TrackAngleCutYZHigh");

  produces< std::vector< recob::Track > >();
  produces< std::vector< anab::T0> >();
  produces< art::Assns< recob::Track, anab::T0 > >();
  produces< std::vector< recob::Hit > >();
  produces< art::Assns< recob::Track, recob::Hit > >();
}

void DiffusionFilter::beginJob()
{

    tree = tfs->make< TTree >("difffiltertree", "diffusion filter tree");
    tree->Branch("track_length", &track_length);
    tree->Branch("track_angle_xz", &track_angle_xz);
    tree->Branch("track_angle_yz", &track_angle_yz);
    tree->Branch("track_has_t0", &track_has_t0);

}

bool DiffusionFilter::filter(art::Event & e)
{
  bool isPass = false;

  art::Handle< std::vector<recob::Track> > trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);

  //  art::Handle< std::vector<recob::OpFlash> > opFlashHandle;
  //  e.getByLabel("simpleFlashCosmic", opFlashHandle);

  //art::FindManyP<simb::MCParticle, anab::BackTrackerMatchingData> mcpsFromTracks(trackHandle, e, fTrackTruthMatcherLabel);
  art::FindManyP<recob::Hit> hitsFromTracks(trackHandle, e, fTrackLabel);
  //art::FindManyP<anab::T0>   t0FromTracks(trackHandle, e, "pandoraCosmicT0RecoLoose");
  art::FindManyP<anab::T0>   t0FromTracks(trackHandle, e, fT0Label);

  // produces a new trackCollection for passing tracks result is that we have
  // only events in which at least one track passes the selection, and only 
  // those tracks are reconstructed
  std::unique_ptr< std::vector<recob::Track> > trackCollection( new std::vector<recob::Track>);
  std::unique_ptr< std::vector<anab::T0> > t0Collection( new std::vector<anab::T0>);
  std::unique_ptr< art::Assns<recob::Track, anab::T0> > trackT0Assn( new art::Assns<recob::Track, anab::T0>);
  std::unique_ptr< std::vector<recob::Hit> > hitCollection( new std::vector<recob::Hit> );
  std::unique_ptr< art::Assns<recob::Track, recob::Hit> > trackHitAssn( new art::Assns<recob::Track, recob::Hit>);

  for (auto const& thisTrack : (*trackHandle)){

    track_length = thisTrack.Length();

    // track angles
    // here we've defined it such that fwd/bwd going tracks are xz = 0, +/- pi
    // downward going tracks are yz = -pi/2
    //
    // the tracks we want are xz = 0, +/- pi. yz = 0, +/- pi
    track_angle_xz = std::atan2(thisTrack.StartDirection().Z(), thisTrack.StartDirection().X()) * 180/TMath::Pi();
    track_angle_yz = std::atan2(thisTrack.StartDirection().Y(), thisTrack.StartDirection().Z()) * 180/TMath::Pi();

    std::vector< art::Ptr<anab::T0> > t0s = t0FromTracks.at(thisTrack.ID());
    
    track_has_t0 = false;
    
    if (t0s.size() == 1)
        track_has_t0 = true;
   
    tree->Fill();

    if (track_length < fTrackLengthCut) continue;
    
    if ((track_angle_xz <= fTrackAngleCutXZHigh && track_angle_xz >= fTrackAngleCutXZLow
          && track_angle_yz <= fTrackAngleCutYZHigh && track_angle_yz >= fTrackAngleCutYZLow) 
        || (track_angle_xz >= (180 - fTrackAngleCutXZHigh) && (track_angle_xz <= (180 - fTrackAngleCutXZLow)) 
          && track_angle_yz >= (180 - fTrackAngleCutYZHigh) && track_angle_yz <= (180 - fTrackAngleCutYZLow))){

      trackCollection->push_back(recob::Track(thisTrack.Trajectory(), 
            thisTrack.ParticleId(),
            thisTrack.Chi2(),
            thisTrack.Ndof(),
            thisTrack.VertexCovarianceLocal5D(),
            thisTrack.EndCovarianceLocal5D(),
            thisTrack.ID()));

      if (t0s.size() != 1) {

        std::cout << "nT0: " << t0s.size() << std::endl; 

        continue;
      }
      else {

        t0Collection->push_back(*(t0s.at(0).get()));

        std::cout << ">> T0 INFORMATION "
          << "\n >> T0 tagged time: " << t0s.at(0)->Time() << std::endl;

      }

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

DEFINE_ART_MODULE(DiffusionFilter)
