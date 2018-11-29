////////////////////////////////////////////////////////////////////////
// Class:       LArDiffusion
// Plugin Type: analyzer (art v2_11_03)
// File:        LArDiffusion_module.cc
//
// Generated at Thu Nov 29 09:47:03 2018 by Adam Lister using cetskelgen
// from cetlib version v3_03_01.
// 
// authors: A. Lister (a.lister1@lancaster.ac.uk)
//          A. Mogan (a.mogan@vols.utk.edu)
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

// LArSoft 
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// art
#include "canvas/Persistency/Common/FindManyP.h"

// ROOT
#include "TH1D.h"

// local
#include "ubana/DiffusionModule/Algorithms/WaveformFunctions.h"

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

  float hit_peak_time;
  int number_ticks_per_bin;
  int tick_window_size;
  int tick_window_left;
  int tick_window_right;

  // fhicl
  std::string track_label;
  std::string wire_label;
  std::string hit_label;
  std::string track_hit_assn;
  std::string track_t0_assn;
  std::string hit_wire_assn;
  bool use_t0tagged_tracks;
  float drift_velocity;
  double hit_GOF_cut;
  int waveform_drift_size;
  int number_drift_bins;

  // manipulation histograms
  TH1D* h_wire_in_window = new TH1D("h_wire_in_window", "", 100, 0, 100);

  // detector properties
  ::detinfo::DetectorProperties const* _detprop;

  // other classes
  diffmod::WaveformFunctions _waveform_func;

};


diffmod::LArDiffusion::LArDiffusion(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{

    // defaults set to be MCC9 defaults
    track_label = p.get< std::string >("TrackLabel", "pandoraCosmic");
    wire_label  = p.get< std::string >("WireLabel", "butcher");
    hit_label   = p.get< std::string >("HitLabel", "gaushit");

    track_hit_assn = p.get< std::string >("TrackHitAssn", "pandoraCosmic");
    track_t0_assn  = p.get< std::string >("TrackT0Assn", "pandoraCosmicT0RecoLoose");
    hit_wire_assn = p.get< std::string >("HitWireAssn", "gaushit");

    use_t0tagged_tracks = p.get< bool >("UseT0TaggedTracks", true); 
    drift_velocity      = p.get< float >("DriftVelocity");
    hit_GOF_cut         = p.get< double >("HitGOFCut", 1.1);
    waveform_drift_size = p.get< int >("WindowSize", 6400); 
    number_drift_bins   = p.get< int > ("NumberDriftBins", 25); 

    _detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    number_ticks_per_bin = waveform_drift_size/number_drift_bins;
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

    // get hit information
    art::Handle< std::vector<recob::Hit> > hit_handle;
    e.getByLabel(hit_label, hit_handle);

    // get associations
    art::FindManyP< recob::Hit > hits_from_tracks(track_handle, e, track_hit_assn);
    art::FindManyP< anab::T0 > t0_from_tracks(track_handle, e, track_t0_assn);
    art::FindManyP< recob::Wire > wire_from_hits(hit_handle, e, hit_wire_assn);

    // loop tracks, get associated hits
    for (size_t i_tr = 0; i_tr < track_ptr_vector.size(); i_tr++){

        art::Ptr< recob::Track > thisTrack = track_ptr_vector.at(i_tr);

        std::vector< art::Ptr< anab::T0 > > t0_from_track = t0_from_tracks.at(thisTrack.key());
        std::vector< art::Ptr< recob::Hit > > hits_from_track = hits_from_tracks.at(thisTrack.key());

        // if a t0 exists (should, after diffusion filtering), then
        // go grab it and get the tick correction value
        if (t0_from_track.size() == 1){

            art::Ptr< anab::T0 > thisT0 = t0_from_track.at(0);

            // correction time is the t0 of the track 
            // + the number of ticks we drop
            // + the tick at which we begin to be in time
            // TODO: modify this to use fhicl where possible
            float correction_time = thisT0->Time()+2400+800;

            std::cout << "[DIFFMOD] Correction Time: " 
                << correction_time << std::endl;

            // loop hits
            for (size_t i_hit = 0; i_hit < hits_from_track.size(); i_hit++){

                art::Ptr< recob::Hit > thisHit = hits_from_track.at(i_hit);

                // if hit selection is not passed then ignore the hit
                if (!_waveform_func.passesHitSelection(thisHit, hit_GOF_cut)) 
                    continue;
             

                // get wire information for hit
                art::Ptr< recob::Wire > wire_from_hit = wire_from_hits.at(thisHit.key()).at(0);            
                
                hit_peak_time = thisHit->PeakTime(); 
                tick_window_size = number_ticks_per_bin;
                tick_window_left  = hit_peak_time - tick_window_size/2;
                tick_window_right = hit_peak_time + tick_window_size/2;

                // make sure that the window stops at the edge of the waveform
                if (tick_window_left < 0){
                    tick_window_left = 0;
                    tick_window_size = tick_window_right; 
                }
                if (tick_window_left > waveform_drift_size){
                    tick_window_right = waveform_drift_size;
                    tick_window_size = tick_window_right - tick_window_left;
                }

                // using the peak time, go to the recob::Wire, and grab the 
                // information within some number of ticks of the peak
                for (int i_tick = tick_window_left; i_tick < tick_window_right; i_tick++){

                    std::cout << wire_from_hit->Signal().at(i_tick) << std::endl;

                }

                // for later
                //detprop->ConvertTicksToX(thisHit->PeakTime(), geo::PlaneID(0,0,2));

                

            }

        }

    }

}

void diffmod::LArDiffusion::beginJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(diffmod::LArDiffusion)
