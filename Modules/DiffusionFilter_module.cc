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
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "art/Persistency/Common/PtrMaker.h"

// ROOT includes
#include "TTree.h"
#include "TMath.h"

// cpp
#include <exception>

// local includes
#include "ubana/DiffusionModule/Algorithms/Utilities.h"

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
    diffmod::Utilities _util;

    // fhicl
    std::string fTrackLabel;
    std::string fT0Label;
    std::string fHitLabel;
    std::string fHitSPAssnLabel;
    double fTrackLengthCut;
    double fTrackAngleCutXZLow;
    double fTrackAngleCutXZHigh;
    double fTrackAngleCutYZLow;
    double fTrackAngleCutYZHigh;

    double thisTrackLength;
    double thisTrackStartX;
    double thisTrackStartX_t0Corr;
    double thisTrackStartY;
    double thisTrackStartZ;
    double thisTrackEndX;
    double thisTrackEndX_t0Corr;
    double thisTrackEndY;
    double thisTrackEndZ;
    double thisTrackT0;
    double thisTrackThetaXZ;
    double thisTrackThetaYZ;
    double thisTrackTheta;
    double thisTrackPhi;
    bool thisTrackIsPassLengthCut;
    bool thisTrackIsPassThetaXZ;
    bool thisTrackIsPassThetaYZ;
    bool thisTrackIsPassAngularCut;
    bool thisTrackIsHasT0;
    bool thisTrackIsSelected;
    bool thisTrackIsAnodeCrosser;
    bool thisTrackIsCathodeCrosser;

    TTree* tree;
    std::vector<double>* trackLength         = nullptr;
    std::vector<double>* trackStartX         = nullptr;
    std::vector<double>* trackStartX_t0Corr  = nullptr;
    std::vector<double>* trackStartY         = nullptr;
    std::vector<double>* trackStartZ         = nullptr;
    std::vector<double>* trackEndX           = nullptr;
    std::vector<double>* trackEndX_t0Corr    = nullptr;
    std::vector<double>* trackEndY           = nullptr;
    std::vector<double>* trackEndZ           = nullptr;
    std::vector<double>* trackT0             = nullptr;
    std::vector<double>* trackThetaXZ        = nullptr;
    std::vector<double>* trackThetaYZ        = nullptr;
    std::vector<double>* trackTheta          = nullptr;
    std::vector<double>* trackPhi            = nullptr;
    std::vector<bool>* trackIsPassLengthCut  = nullptr;
    std::vector<bool>* trackIsPassThetaXZ    = nullptr;
    std::vector<bool>* trackIsPassThetaYZ    = nullptr;
    std::vector<bool>* trackIsPassAngularCut = nullptr;
    std::vector<bool>* trackIsHasT0          = nullptr;
    std::vector<bool>* trackIsSelected       = nullptr;
    std::vector<bool>* trackIsCathodeCrosser = nullptr;
    std::vector<bool>* trackIsAnodeCrosser   = nullptr;

};


DiffusionFilter::DiffusionFilter(fhicl::ParameterSet const & p)
  // :
  // More initializers here.
{

  fTrackLabel          = p.get<std::string> ("TrackLabel");
  fT0Label             = p.get<std::string> ("T0Label");
  fHitLabel            = p.get<std::string> ("HitLabel");
  fHitSPAssnLabel      = p.get<std::string> ("HitSPAssnLabel");
  fTrackLengthCut      = p.get<double> ("TrackLengthCut");
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
  produces< std::vector< recob::SpacePoint > >();
  produces< art::Assns< recob::Hit, recob::SpacePoint > >();

}

void DiffusionFilter::beginJob()
{


  std::cout << "DiffusionFilter::beginJob() begin" << std::endl;

  tree = tfs->make< TTree >("difffiltertree", "diffusion filter tree");

  tree->Branch("trackLength"           , &trackLength);
  tree->Branch("trackStartX"           , &trackStartX);
  tree->Branch("trackStartX_t0Corr"    , &trackStartX_t0Corr);
  tree->Branch("trackStartY"           , &trackStartY);
  tree->Branch("trackStartZ"           , &trackStartZ);
  tree->Branch("trackEndX"             , &trackEndX);
  tree->Branch("trackEndX_t0Corr"      , &trackEndX_t0Corr);
  tree->Branch("trackEndY"             , &trackEndY);
  tree->Branch("trackEndZ"             , &trackEndZ);
  tree->Branch("trackT0"               , &trackT0);
  tree->Branch("trackTheta"            , &trackTheta);
  tree->Branch("trackPhi"              , &trackPhi);
  tree->Branch("trackThetaXZ"          , &trackThetaXZ);
  tree->Branch("trackThetaYZ"          , &trackThetaYZ);
  tree->Branch("trackIsPassLengthCut"  , &trackIsPassLengthCut);
  tree->Branch("trackIsPassThetaXZ"    , &trackIsPassThetaXZ);
  tree->Branch("trackIsPassThetaYZ"    , &trackIsPassThetaYZ);
  tree->Branch("trackIsPassAngularCut" , &trackIsPassAngularCut);
  tree->Branch("trackIsHasT0"          , &trackIsHasT0);
  tree->Branch("trackIsSelected"       , &trackIsSelected);
  tree->Branch("trackIsAnodeCrosser"   , &trackIsAnodeCrosser);
  tree->Branch("trackIsCathodeCrosser" , &trackIsCathodeCrosser);

  std::cout << "DiffusionFilter::beginJob() end" << std::endl;

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

  art::Handle< std::vector<recob::Hit> > hitHandle;
  e.getByLabel(fHitLabel, hitHandle);
  std::vector< art::Ptr< recob::Hit > > hitPtrVector;
  art::fill_ptr_vector(hitPtrVector, hitHandle);

  art::FindManyP< recob::Hit >        hitsFromTracks(trackHandle, e, fTrackLabel);
  art::FindManyP< anab::T0 >          t0FromTracks  (trackHandle, e, fT0Label);
  art::FindManyP< recob::SpacePoint > spFromHits    (hitHandle  , e, fHitSPAssnLabel); 

  // produces a new trackCollection for passing tracks result is that we have
  // only events in which at least one track passes the selection, and only 
  // those tracks are reconstructed
  std::unique_ptr< std::vector<recob::Track> >                 trackCollection     ( new std::vector<recob::Track>);
  std::unique_ptr< std::vector<anab::T0> >                     t0Collection        ( new std::vector<anab::T0>);
  std::unique_ptr< art::Assns<recob::Track, anab::T0> >        trackT0Assn         ( new art::Assns<recob::Track, anab::T0>);
  std::unique_ptr< std::vector<recob::Hit> >                   hitCollection       ( new std::vector<recob::Hit> );
  std::unique_ptr< art::Assns<recob::Track, recob::Hit> >      trackHitAssn        ( new art::Assns<recob::Track, recob::Hit>);
  std::unique_ptr< std::vector<recob::SpacePoint> >            spacePointCollection( new std::vector<recob::SpacePoint>);
  std::unique_ptr< art::Assns<recob::Hit, recob::SpacePoint> > hitSpacePointAssn   ( new art::Assns<recob::Hit, recob::SpacePoint>); 

  art::PtrMaker< recob::Track >      makeTrackPtr(e);
  art::PtrMaker< recob::Hit >        makeHitPtr(e);
  art::PtrMaker< anab::T0 >          makeT0Ptr(e);
  art::PtrMaker< recob::SpacePoint > makeSpacePointPtr(e); 

  for (size_t iTrack = 0; iTrack < trackPtrVector.size(); iTrack++){
    art::Ptr< recob::Track > thisTrack = trackPtrVector.at(iTrack);
    std::vector< art::Ptr<anab::T0> > t0s = t0FromTracks.at(thisTrack.key());

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
          + std::to_string(thisTrack.key())
          + " has "
          + std::to_string(t0s.size())
          + " associated t0s. That can't be right");
      throw std::logic_error(errMsg);
    }

    thisTrackLength  = thisTrack->Length();
    thisTrackStartX  = thisTrack->Start().X();
    thisTrackStartY  = thisTrack->Start().Y();
    thisTrackStartZ  = thisTrack->Start().Z();
    thisTrackEndX    = thisTrack->End().X();
    thisTrackEndY    = thisTrack->End().Y();
    thisTrackEndZ    = thisTrack->End().Z();
    thisTrackTheta   = thisTrack->Theta();
    thisTrackPhi     = thisTrack->Phi();
    thisTrackThetaXZ = _util.getThetaXZ(thisTrack);
    thisTrackThetaYZ = _util.getThetaYZ(thisTrack);

    if (thisTrackIsHasT0){
      thisTrackT0 = t0s.at(0)->Time();
      thisTrackStartX_t0Corr = thisTrackStartX - thisTrackT0*1.098;
      thisTrackEndX_t0Corr   = thisTrackEndX   - thisTrackT0*1.098;
    }
    else{
      thisTrackT0 = -1e-9;
      thisTrackStartX_t0Corr = thisTrackStartX;
      thisTrackEndX_t0Corr = thisTrackEndX;
    }
    
    thisTrackIsPassLengthCut = (thisTrackLength > fTrackLengthCut);

    thisTrackIsPassThetaXZ = (thisTrackThetaXZ > fTrackAngleCutXZLow &&
                              thisTrackThetaXZ < fTrackAngleCutXZHigh);

    thisTrackIsPassThetaYZ = (thisTrackThetaYZ > fTrackAngleCutYZLow &&
                              thisTrackThetaYZ < fTrackAngleCutYZHigh);

    thisTrackIsPassAngularCut = thisTrackIsPassThetaXZ && thisTrackIsPassThetaYZ;

    // if the track starts or ends near the anode, and neither end is located near to the cathode
    if ((thisTrackStartX_t0Corr < 5 || thisTrackEndX_t0Corr < 5) && (thisTrackStartX_t0Corr < 250.0 && thisTrackEndX_t0Corr < 250.0)){
      thisTrackIsAnodeCrosser   = 1;
      thisTrackIsCathodeCrosser = 0;
    }
    // if the track starts or ends near the cathode, and neither end is located near to the anode
    else if ((thisTrackStartX_t0Corr > 250 || thisTrackEndX_t0Corr > 250) && (thisTrackStartX_t0Corr > 5.0 && thisTrackEndX_t0Corr > 5.0)){
      thisTrackIsAnodeCrosser   = 0;
      thisTrackIsCathodeCrosser = 1;
    }
    // catch the case where the t0_tagged track travels nearly the full width of the detector --- we don't care about these things
    else{
      thisTrackIsAnodeCrosser   = 0;
      thisTrackIsCathodeCrosser = 0;
    }

    trackLength          ->push_back(thisTrackLength);
    trackStartX          ->push_back(thisTrackStartX);
    trackStartX_t0Corr   ->push_back(thisTrackStartX_t0Corr);
    trackStartY          ->push_back(thisTrackStartY);
    trackStartZ          ->push_back(thisTrackStartZ);
    trackEndX            ->push_back(thisTrackEndX);
    trackEndX_t0Corr     ->push_back(thisTrackEndX_t0Corr);
    trackEndY            ->push_back(thisTrackEndY);
    trackEndZ            ->push_back(thisTrackEndZ);
    trackT0              ->push_back(thisTrackT0);
    trackTheta           ->push_back(thisTrackTheta);
    trackPhi             ->push_back(thisTrackPhi);
    trackThetaXZ         ->push_back(thisTrackThetaXZ);
    trackThetaYZ         ->push_back(thisTrackThetaYZ);
    trackIsPassLengthCut ->push_back(thisTrackIsPassLengthCut);
    trackIsPassThetaXZ   ->push_back(thisTrackIsPassThetaXZ);
    trackIsPassThetaYZ   ->push_back(thisTrackIsPassThetaYZ);
    trackIsPassAngularCut->push_back(thisTrackIsPassAngularCut);
    trackIsHasT0         ->push_back(thisTrackIsHasT0);
    trackIsSelected      ->push_back((thisTrackIsPassLengthCut && 
                                      thisTrackIsPassAngularCut && 
                                      thisTrackIsHasT0));
    trackIsCathodeCrosser->push_back(thisTrackIsCathodeCrosser);
    trackIsAnodeCrosser  ->push_back(thisTrackIsAnodeCrosser);

    if (trackIsSelected->at(iTrack)){

      isPass = true;

      // now create collections for the event
      recob::Track trackForCollection = *(thisTrack.get());

      trackCollection->push_back(trackForCollection); 

      art::Ptr< recob::Track > trackForCollectionPtr 
        = makeTrackPtr(trackCollection->size()-1);


      std::vector< art::Ptr<recob::Hit> > hits;
      std::vector< art::Ptr< recob::Hit > > hitPtrCollection;
      if ((int)hitsFromTracks.at(thisTrack.key()).size() > 0){
        hits  = hitsFromTracks.at(thisTrack.key());

        for (art::Ptr<recob::Hit>& thisHit : hits){

          recob::Hit hitForCollection = *(thisHit.get());
          hitCollection->push_back(hitForCollection);

          art::Ptr< recob::Hit > hitForCollectionPtr 
            = makeHitPtr(hitCollection->size()-1);

          hitPtrCollection.push_back(hitForCollectionPtr);

          // get spacepoints
          std::vector< art::Ptr< recob::SpacePoint > > spacePoints;
          std::vector< art::Ptr< recob::SpacePoint > > spacePointPtrCollection;
          if ((int)spFromHits.at(thisHit.key()).size() > 0){
            spacePoints = spFromHits.at(thisHit.key());

            for (art::Ptr< recob::SpacePoint >& thisSpacePoint : spacePoints){
               recob::SpacePoint spacePointForCollection = *(thisSpacePoint.get());
               spacePointCollection->push_back(spacePointForCollection);

               art::Ptr< recob::SpacePoint > spacePointForCollectionPtr
                 = makeSpacePointPtr(spacePointCollection->size()-1);

               spacePointPtrCollection.push_back(spacePointForCollectionPtr);
            }

            util::CreateAssn(
                *this,
                e,
                hitForCollectionPtr,
                spacePointPtrCollection,
                *hitSpacePointAssn);

          }


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


  std::cout << "DiffusionFilter::filter() --- filling tree..." << std::endl;
  tree->Fill();
  std::cout << "DiffusionFilter::filter() --- filled." << std::endl;

  e.put(std::move(trackCollection));
  e.put(std::move(t0Collection));
  e.put(std::move(hitCollection));
  e.put(std::move(spacePointCollection));
  e.put(std::move(trackT0Assn));
  e.put(std::move(trackHitAssn));
  e.put(std::move(hitSpacePointAssn));

  std::cout << "DiffusionFilter::filter() --- moved dataProducts to the event " << std::endl;

  return isPass;

}

void DiffusionFilter::emptyVectors(){
  trackLength           -> resize(0);
  trackStartX           -> resize(0);
  trackStartX_t0Corr    -> resize(0);
  trackStartY           -> resize(0);
  trackStartZ           -> resize(0);
  trackEndX             -> resize(0);
  trackEndX_t0Corr      -> resize(0);
  trackEndY             -> resize(0);
  trackEndZ             -> resize(0);
  trackT0               -> resize(0);
  trackThetaXZ          -> resize(0);
  trackThetaYZ          -> resize(0);
  trackTheta            -> resize(0);
  trackPhi              -> resize(0);
  trackIsPassLengthCut  -> resize(0);
  trackIsPassThetaXZ    -> resize(0);
  trackIsPassThetaYZ    -> resize(0);
  trackIsPassAngularCut -> resize(0);
  trackIsHasT0          -> resize(0);
  trackIsSelected       -> resize(0);
  trackIsCathodeCrosser -> resize(0);
  trackIsAnodeCrosser   -> resize(0);
}


DEFINE_ART_MODULE(DiffusionFilter)
