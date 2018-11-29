////////////////////////////////////////////////////////////////////////
// Class:       LArDiffusion
// Plugin Type: analyzer (art v2_11_03)
// File:        LArDiffusion_module.cc
//
// Generated at Thu Nov 29 09:47:03 2018 by Adam Lister using cetskelgen
// from cetlib version v3_03_01.
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

// LArSoft 
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/AnalysisBase/T0.h"

// art
#include "canvas/Persistency/Common/FindManyP.h"

namespace diffmod {
  class LArDiffusion;
}


class diffmod::LArDiffusion : public art::EDAnalyzer {
public:
  explicit LArDiffusion(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LArDiffusion(LArDiffusion const &) = delete;
  LArDiffusion(LArDiffusion &&) = delete;
  LArDiffusion & operator = (LArDiffusion const &) = delete;
  LArDiffusion & operator = (LArDiffusion &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // variables
  int run;
  int sub_run;
  int event;
  int is_real_data;

  bool use_t0tagged_tracks;

  // fhicl
  std::string track_label;
  std::string wire_label;

  std::string track_hit_assn;
  std::string track_t0_assn;

  // Declare member data here.

};


diffmod::LArDiffusion::LArDiffusion(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

    // defaults set to be MCC9 defaults
    track_label = p.get< std::string >("TrackLabel", "pandoraCosmic");
    wire_label  = p.get< std::string >("WireLabel", "butcher");

    track_hit_assn = p.get< std::string >("TrackHitAssn", "pandoraCosmic");
    track_t0_assn  = p.get< std::string >("TrackT0Assn", "pandoraCosmicT0RecoLoose");

    use_t0tagged_tracks = p.get< bool >("UseT0TaggedTracks", true); 
}

void diffmod::LArDiffusion::analyze(art::Event const & e)
{
    run = e.run();
    sub_run = e.subRun();
    event = e.event();
    is_real_data = e.isRealData();

    std::cout << "[DIFFMOD] --- Processing event " 
        << run << "." << sub_run << "." << event << std::endl;

    // get track information
    art::Handle< std::vector<recob::Track> > track_handle;
    e.getByLabel(track_label, track_handle);
    std::vector< art::Ptr<recob::Track> > track_ptr_vector;
    art::fill_ptr_vector(track_ptr_vector, track_handle);

    // get wire information
    art::Handle< std::vector<recob::Wire> > wire_handle;
    e.getByLabel(wire_label, wire_handle);

    // get associations
    art::FindManyP< recob::Hit > hits_from_tracks(track_handle, e, track_hit_assn);
    art::FindManyP< anab::T0 > t0_from_tracks(track_handle, e, track_t0_assn);

    // loop tracks, get associated hits
    for (size_t i_tr = 0; i_tr < track_ptr_vector.size(); i_tr++){

        art::Ptr< recob::Track > thisTrack = track_ptr_vector.at(i_tr);

        std::vector< art::Ptr< anab::T0 > > t0_from_track = t0_from_tracks.at(thisTrack.key());
        std::vector< art::Ptr< recob::Hit > > hit_from_track = hits_from_tracks.at(thisTrack.key());

    }

}

void diffmod::LArDiffusion::beginJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(diffmod::LArDiffusion)