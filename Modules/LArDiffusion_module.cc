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
  int waveform_drift_size;

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
  int waveform_size;
  int waveform_intime_start;
  int waveform_intime_end;
  int number_drift_bins;
  float peak_finder_threshold;
  int number_dropped_ticks;

  // manipulation histograms
  // binnings are placeholders and will be modified later

  // histogram opened around the hit peak value 
  TH1D* h_wire_in_window = new TH1D("h_wire_in_window", "", 100, 0, 100);

  // after baseline correcting
  TH1D* h_wire_baseline_corrected = new TH1D("h_wire_baseline_corrected", "", 100, 0, 100);

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
    waveform_size = p.get< int >("WaveformSize", 6400); 
    waveform_intime_start = p.get< int >("WaveformIntimeStart", 800);
    waveform_intime_end = p.get< int >("WaveformIntimeEnd", 5600);
    number_drift_bins   = p.get< int >("NumberDriftBins", 25); 
    peak_finder_threshold = p.get< float >("PeakFinderThreshold", 3.0);
    number_dropped_ticks = p.get< int >("NumberDroppedTicks", 2400);

    waveform_drift_size = waveform_intime_end - waveform_intime_start;
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
            float t0 = thisT0->Time();
            float t0_tick = t0 * 2;
            float correction_x = t0 * drift_velocity;

            std::cout << "t0: " << t0 << std::endl;
            std::cout << "drift velocity: " << drift_velocity << std::endl;
            std::cout << "correction_x: " << correction_x << std::endl;
            std::cout << "uncorr: track start x: " << thisTrack->Start().X() << " end x: " << thisTrack->End().X() << std::endl;
            std::cout << "corr+: track start x: " << thisTrack->Start().X()+correction_x << " end x: " << thisTrack->End().X()+correction_x << std::endl;
            std::cout << "corr-: track start x: " << thisTrack->Start().X()-correction_x << " end x: " << thisTrack->End().X()-correction_x << std::endl;

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
                if (tick_window_left > waveform_size){
                    tick_window_right = waveform_size;
                    tick_window_size = tick_window_right - tick_window_left;
                }

                h_wire_in_window->SetBins(tick_window_size, tick_window_left, tick_window_right);

                // using the peak time, go to the recob::Wire, and grab the 
                // information within some number of ticks of the peak
                //
                // also check to ensure that we only have one peak above
                // threshold
                int peak_counter = 0;
                for (int i_tick = tick_window_left; i_tick < tick_window_right; i_tick++){

                    float value = wire_from_hit->Signal().at(i_tick);

                    h_wire_in_window->SetBinContent(i_tick - tick_window_left, value);

                    if (value > peak_finder_threshold){

                        // make sure we only look for the peak
                        if (wire_from_hit->Signal().at(i_tick-1) < value
                         && wire_from_hit->Signal().at(i_tick+1) < value){
                            std::cout << "found peak!" << std::endl;
                        
                            peak_counter++;
                        }

                    }

                }

                if (peak_counter != 1){
                    std::cout << "peak counter: " << peak_counter << ", skipping channel" << std::endl;
                    continue;
                }

                // get peak bin tick after t0 correction
                int maximum_tick = h_wire_in_window->GetMaximumBin()
                    + tick_window_left - t0_tick;

                // now the magic: 
                // loop over the drift bins and check to see if the 
                // t0-corrected pulse center is in each bin
                // if it is then sum the pulses

                std::cout << "maximum tick: " << maximum_tick << std::endl;

                for (int bin_it = 0; bin_it < number_drift_bins; bin_it++){

                    // set binning for histograms
                    h_wire_baseline_corrected->SetBins(h_wire_in_window->GetNbinsX(),
                            waveform_intime_start + (bin_it) * number_ticks_per_bin,
                            waveform_intime_start + (bin_it + 1) * number_ticks_per_bin);

                    if (maximum_tick >= (bin_it * number_ticks_per_bin)
                        && maximum_tick < ((bin_it + 1) * number_ticks_per_bin)) {

                        h_wire_in_window->GetXaxis()->SetLimits(bin_it * number_ticks_per_bin, (bin_it +1) * number_ticks_per_bin);

                    }


                    for (int i = 0; i < h_wire_in_window->GetNbinsX(); i++){

                        std::cout << "t0corr: " << i << " " << h_wire_in_window->GetBinContent(i) << std::endl;

                    }

                    h_wire_baseline_corrected = _waveform_func.applyGlobalBaselineCorrection(h_wire_in_window, h_wire_baseline_corrected); 
                    
                    for (int i = 0; i < h_wire_baseline_corrected->GetNbinsX(); i++){

                        std::cout << "baseline corrected: " << i << " " << h_wire_baseline_corrected->GetBinContent(i) << std::endl;

                    }

                }

            }

        }

    }

}

void diffmod::LArDiffusion::beginJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(diffmod::LArDiffusion)
