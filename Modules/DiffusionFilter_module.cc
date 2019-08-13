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
#include "art/Persistency/Common/PtrMaker.h"

// ROOT includes
#include "TTree.h"
#include "TMath.h"

// cpp
#include <exception>

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

    // clear vectors
    void emptyVectors();

  private:

    // service and class instances
    art::ServiceHandle< art::TFileService > tfs;

    // fhicl
    std::string fTrackLabel;
    std::string fT0Label;
    double fTrackLengthCut;
    double fTrackAngleCutXZLow;
    double fTrackAngleCutXZHigh;
    double fTrackAngleCutYZLow;
    double fTrackAngleCutYZHigh;

    double thisTrackLength;
    double thisTrackThetaXZ;
    double thisTrackThetaYZ;
    double thisTrackTheta;
    double thisTrackPhi;
    bool thisTrackIsPassLengthCut;
    bool thisTrackIsPassAngularCut;
    bool thisTrackIsHasT0;
    bool thisTrackIsSelected;

    TTree* tree;
    std::vector<double>* trackLength = nullptr;
    std::vector<double>* trackThetaXZ = nullptr;
    std::vector<double>* trackThetaYZ = nullptr;
    std::vector<double>* trackTheta = nullptr;
    std::vector<double>* trackPhi = nullptr;
    std::vector<bool>* trackIsPassLengthCut = nullptr;
    std::vector<bool>* trackIsPassAngularCut = nullptr;
    std::vector<bool>* trackIsHasT0 = nullptr;
    std::vector<bool>* trackIsSelected = nullptr;

};


DiffusionFilter::DiffusionFilter(fhicl::ParameterSet const & p)
  // :
  // More initializers here.
{

  fTrackLabel = p.get<std::string> ("TrackLabel");
  fT0Label = p.get<std::string> ("T0Label");
  fTrackLengthCut = p.get<double> ("TrackLengthCut");
  fTrackAngleCutXZLow  = p.get<double> ("TrackAngleCutXZLow");
  fTrackAngleCutXZHigh = p.get<double> ("TrackAngleCutXZHigh");
  fTrackAngleCutYZLow  = p.get<double> ("TrackAngleCutYZLow");
  fTrackAngleCutYZHigh = p.get<double> ("TrackAngleCutYZHigh");

  std::cout << "--- Printing Configuration: " << std::endl;
  std::cout << "TrackLabel          : " << fTrackLabel << std::endl;
  std::cout << "fT0Label            : " << fT0Label << std::endl;
  std::cout << "TrackLengthCut      : " << fTrackLengthCut << std::endl;
  std::cout << "TrackAngleCutXZLow  : " << fTrackAngleCutXZLow << std::endl;
  std::cout << "TrackAngleCutXZHigh : " << fTrackAngleCutXZHigh << std::endl;
  std::cout << "TrackAngleCutYZLow  : " << fTrackAngleCutYZLow << std::endl;
  std::cout << "TrackAngleCutYZHigh : " << fTrackAngleCutYZHigh << std::endl;
  std::cout << "---------------------------" << std::endl;

  produces< std::vector< recob::Track > >();
  produces< std::vector< anab::T0> >();
  produces< art::Assns< recob::Track, anab::T0 > >();
  produces< std::vector< recob::Hit > >();
  produces< art::Assns< recob::Track, recob::Hit > >();

}

void DiffusionFilter::beginJob()
{

  tree = tfs->make< TTree >("difffiltertree", "diffusion filter tree");

  tree->Branch("trackLength", &trackLength);
  tree->Branch("trackTheta", &trackTheta);
  tree->Branch("trackPhi", &trackPhi);
  tree->Branch("trackThetaXZ", &trackThetaXZ);
  tree->Branch("trackThetaYZ", &trackThetaYZ);
  tree->Branch("trackIsPassLengthCut", &trackIsPassLengthCut);
  tree->Branch("trackIsPassAngularCut", &trackIsPassAngularCut);
  tree->Branch("trackIsHasT0", &trackIsHasT0);
  tree->Branch("trackIsSelected", &trackIsSelected);

}

bool DiffusionFilter::filter(art::Event & e)
{

  // initialise variables
  this->emptyVectors();
  bool isPass = false;

  art::Handle< std::vector<recob::Track> > trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);
  std::vector< art::Ptr< recob::Track > > trackPtrVector;
  art::fill_ptr_vector(trackPtrVector, trackHandle);

  art::FindManyP<recob::Hit> hitsFromTracks(trackHandle, e, fTrackLabel);
  art::FindManyP<anab::T0>   t0FromTracks(trackHandle, e, fT0Label);

  // produces a new trackCollection for passing tracks result is that we have
  // only events in which at least one track passes the selection, and only 
  // those tracks are reconstructed
  std::unique_ptr< std::vector<recob::Track> >            trackCollection( new std::vector<recob::Track>);
  std::unique_ptr< std::vector<anab::T0> >                t0Collection( new std::vector<anab::T0>);
  std::unique_ptr< art::Assns<recob::Track, anab::T0> >   trackT0Assn( new art::Assns<recob::Track, anab::T0>);
  std::unique_ptr< std::vector<recob::Hit> >              hitCollection( new std::vector<recob::Hit> );
  std::unique_ptr< art::Assns<recob::Track, recob::Hit> > trackHitAssn( new art::Assns<recob::Track, recob::Hit>);

  art::PtrMaker< recob::Track > makeTrackPtr(e);
  art::PtrMaker< recob::Hit >   makeHitPtr(e);
  art::PtrMaker< anab::T0 >     makeT0Ptr(e);

  for (size_t iTrack = 0; iTrack < trackPtrVector.size(); iTrack++){
    art::Ptr< recob::Track > thisTrack = trackPtrVector.at(iTrack);

    thisTrackLength = thisTrack->Length();
    thisTrackTheta = thisTrack->Theta();
    thisTrackPhi = thisTrack->Phi();
    thisTrackThetaXZ 
      = std::atan2(
          thisTrack->StartDirection().X(), 
          thisTrack->StartDirection().Z()) 
      * 180/TMath::Pi();
    thisTrackThetaYZ 
      = std::atan2(
          thisTrack->StartDirection().Y(), 
          thisTrack->StartDirection().Z()) 
      * 180/TMath::Pi();

    std::vector< art::Ptr<anab::T0> > t0s = t0FromTracks.at(thisTrack.key());

    thisTrackIsPassLengthCut = (thisTrackLength > fTrackLengthCut);

    thisTrackIsPassAngularCut = 
      (( thisTrackThetaXZ <= fTrackAngleCutXZHigh 
         && thisTrackThetaXZ >= fTrackAngleCutXZLow
         && thisTrackThetaYZ <= fTrackAngleCutYZHigh 
         && thisTrackThetaYZ >= fTrackAngleCutYZLow)) 
      ||
      (( thisTrackThetaXZ >= (180 - fTrackAngleCutXZHigh) 
         && thisTrackThetaXZ <= (180 - fTrackAngleCutXZLow)) 
       && thisTrackThetaYZ >= (180 - fTrackAngleCutYZHigh) 
       && thisTrackThetaYZ <= (180 - fTrackAngleCutYZLow));

    // make sure that the track has at least one associated t0
    if (t0s.size() == 1){
      thisTrackIsHasT0 = true;

      anab::T0 t0ForCollection = *((t0s.at(0)).get());
      t0Collection->push_back(t0ForCollection);
    }
    else if (t0s.size() == 0)
      thisTrackIsHasT0 = false;
    else {
      std::string errMsg(
          "Track "
          + std::to_string(thisTrack->ID())
          + " has "
          + std::to_string(t0s.size())
          + " associated t0s. That can't be right");
      throw std::logic_error(errMsg);
    }

    trackLength          ->push_back(thisTrackLength);
    trackTheta           ->push_back(thisTrackTheta);
    trackPhi             ->push_back(thisTrackPhi);
    trackThetaXZ         ->push_back(thisTrackThetaXZ);
    trackThetaYZ         ->push_back(thisTrackThetaYZ);
    trackIsPassLengthCut ->push_back(thisTrackIsPassLengthCut);
    trackIsPassAngularCut->push_back(thisTrackIsPassAngularCut);
    trackIsHasT0         ->push_back(thisTrackIsHasT0);
    trackIsSelected      ->push_back((thisTrackIsPassLengthCut && thisTrackIsPassAngularCut && thisTrackIsHasT0));

    if (trackIsSelected->at(iTrack)){

      isPass = true;

      // now create collections for the event
      recob::Track trackForCollection = *(thisTrack.get());

      trackCollection->push_back(trackForCollection); 

      art::Ptr< recob::Track > trackForCollectionPtr 
        = makeTrackPtr(trackCollection->size()-1);


      std::vector< art::Ptr<recob::Hit> > hits;
      std::vector< art::Ptr< recob::Hit > > hitPtrCollection;
      if ((int)hitsFromTracks.at(thisTrack->ID()).size() > 0){
        hits  = hitsFromTracks.at(thisTrack->ID());

        for (art::Ptr<recob::Hit>& thisHit : hits){

          recob::Hit hitForCollection = *(thisHit.get());
          hitCollection->push_back(hitForCollection);

          art::Ptr< recob::Hit > hitForCollectionPtr 
            = makeHitPtr(hitCollection->size()-1);

          hitPtrCollection.push_back(hitForCollectionPtr);
        }
      }


      util::CreateAssn(
          *this, 
          e, 
          trackForCollectionPtr, 
          hitPtrCollection,
          *trackHitAssn);

      if (trackIsHasT0){
        art::Ptr< anab::T0 > t0ForCollectionPtr 
          = makeT0Ptr(t0Collection->size()-1);

        util::CreateAssn(
            *this,
            e,
            t0ForCollectionPtr,
            trackForCollectionPtr,
            *trackT0Assn);
      }
    }
  }

  tree->Fill();

  e.put(std::move(trackCollection));
  e.put(std::move(t0Collection));
  e.put(std::move(hitCollection));
  e.put(std::move(trackT0Assn));
  e.put(std::move(trackHitAssn));

  return isPass;

}

void DiffusionFilter::emptyVectors(){
  trackLength->resize(0);
  trackThetaXZ->resize(0);
  trackThetaYZ->resize(0);
  trackTheta->resize(0);
  trackPhi->resize(0);
  trackIsPassLengthCut->resize(0);
  trackIsPassAngularCut->resize(0);
  trackIsHasT0->resize(0);
  trackIsSelected->resize(0);
}


DEFINE_ART_MODULE(DiffusionFilter)
