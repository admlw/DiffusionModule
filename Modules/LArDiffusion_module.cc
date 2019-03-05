////////////////////////////////////////////////////////////////////////
// Class:       LArDiffusion
// Plugin Type: analyzer (art v2_11_03)
// File:        LArDiffusion_module.cc
//
// Generated at Thu Nov 29 09:47:03 2018 by Adam Lister using cetskelgen
// from cetlib version v3_03_01.
// 
// authors: A. Lister (a.lister1@lancaster.ac.uk)
//          A. Mogan (amogan@vols.utk.edu)
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
#include "art/Framework/Services/Optional/TFileService.h"

// ROOT
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TStyle.h"

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

        art::ServiceHandle< art::TFileService > tfs;

        TTree* difftree;

        // variables
        int run;
        int sub_run;
        int event;
        int is_real_data;

        float hit_peak_time;
        float t0;
        float t0_tick;
        float track_x_correction;
        double pulse_height;
        double mean;
        double sigma;       
        double fit_chisq;   
        double waveform_x_correction;
        int bin_no;

        int number_ticks_per_bin;
        int tick_window_size;
        int tick_window_left;
        int tick_window_right;
        int waveform_drift_size;
        int track_startX;
        

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
        float drift_distance;
        float peak_finder_threshold;
        int number_dropped_ticks;

        // manipulation histograms
        // binnings are placeholders and will be modified later

        // histogram opened around the hit peak value 
        TH1D* h_wire_in_window = tfs->make<TH1D>("h_wire_in_window", "", 100, 0, 100);

        // after baseline correcting
        TH1D* h_wire_baseline_corrected = tfs->make<TH1D>("h_wire_baseline_corrected", "", 100, 0, 100);

        // Define troubleshooting histograms
        TH1D *h_trackLength;
        TH1D *h_cosTheta;
        TH1D *h_startX;
        TH1D *h_sigma;
        TH1D *h_pulseHeight;
        TH1D *h_nWvfmsInBin;

        TH2D *h_driftVsigma;
        TH2D *h_driftVPulseHeight;

        // output histograms
        std::vector<TH1D*> h_summed_wire_info_per_bin; 

        // other classes
        diffmod::WaveformFunctions _waveform_func;

};


diffmod::LArDiffusion::LArDiffusion(fhicl::ParameterSet const & p)
    :
        EDAnalyzer(p)  // ,
        // More initializers here.
{

    // defaults set to be MCC9 defaults
    track_label = p.get< std::string >("TrackLabel", "pandora");
    wire_label  = p.get< std::string >("WireLabel", "butcher");
    hit_label   = p.get< std::string >("HitLabel", "gaushit");

    track_hit_assn = p.get< std::string >("TrackHitAssn", "pandora");
    track_t0_assn  = p.get< std::string >("TrackT0Assn", "t0reco");
    hit_wire_assn  = p.get< std::string >("HitWireAssn", "gaushit");

    use_t0tagged_tracks   = p.get< bool >("UseT0TaggedTracks", true);
    drift_velocity        = p.get< float >("DriftVelocity");
    hit_GOF_cut           = p.get< double >("HitGOFCut", 1.1);
    waveform_size         = p.get< int >("WaveformSize", 6400);
    waveform_intime_start = p.get< int >("WaveformIntimeStart", 800);
    waveform_intime_end   = p.get< int >("WaveformIntimeEnd", 5400);
    number_drift_bins     = p.get< int >("NumberDriftBins", 25);
    number_dropped_ticks  = p.get< int >("NumberDroppedTicks", 2400);
    peak_finder_threshold = p.get< float >("PeakFinderThreshold", 3.0);

    waveform_drift_size  = waveform_intime_end - waveform_intime_start; // 4600
    number_ticks_per_bin = waveform_drift_size/number_drift_bins; // 184
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

        h_trackLength->Fill(thisTrack->Length() );
        h_cosTheta->Fill(thisTrack->Theta() );
        h_startX->Fill(thisTrack->Start().X() );

        std::vector< art::Ptr< anab::T0 > > t0_from_track;
        if (use_t0tagged_tracks) {
            //std::vector< art::Ptr< anab::T0 > > t0_from_track = t0_from_tracks.at(thisTrack.key());
            t0_from_track = t0_from_tracks.at(thisTrack.key());
        }

        std::vector< art::Ptr< recob::Hit > > hits_from_track = hits_from_tracks.at(thisTrack.key());

        // if a t0 exists (should, if using diffusion-filtered sample), then
        // go grab it and get the tick correction value
        if (t0_from_track.size() == 1 || !use_t0tagged_tracks){
            // correction time is the t0 of the track 
            // + the number of ticks we drop
            // + the tick at which we begin to be in time
            if (use_t0tagged_tracks) {
              art::Ptr< anab::T0 > thisT0 = t0_from_track.at(0);
              t0 = thisT0->Time();
              track_x_correction = t0 * drift_velocity;
            }
            else {
              t0 = 800.; // 800 microsecond default value for single muon samples
              track_x_correction = 0;
            }

            t0_tick = t0 * 2;

            /*
            std::cout << "t0: " << t0 << std::endl;
            std::cout << "drift velocity: " << drift_velocity << std::endl;
            std::cout << "track_x_correction: " << track_x_correction << std::endl;
            std::cout << "uncorr: track start x: " << thisTrack->Start().X() << " end x: " << thisTrack->End().X() << std::endl;
            std::cout << "corr+: track start x: " << thisTrack->Start().X()+track_x_correction << " end x: " << thisTrack->End().X()+track_x_correction << std::endl;
            std::cout << "corr-: track start x: " << thisTrack->Start().X()-track_x_correction << " end x: " << thisTrack->End().X()-track_x_correction << std::endl;
            */

            // loop hits
            for (size_t i_hit = 0; i_hit < hits_from_track.size(); i_hit++){

                art::Ptr< recob::Hit > thisHit = hits_from_track.at(i_hit);

                // if hit selection is not passed then ignore the hit
                if (!_waveform_func.passesHitSelection(thisHit, hit_GOF_cut)) continue;

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
                //if (tick_window_left > waveform_size){
                if (tick_window_right > waveform_size){
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

                    // TODO: Check this value. Where does it come from? Is it reasonable?
                    if (value > peak_finder_threshold){ // 3.0 by default

                        // define peak search region
                        if (tick_window_left == 0 || tick_window_right == (int)wire_from_hit->Signal().size())
                          continue;

                        // make sure we only look for the peak
                        if (wire_from_hit->Signal().at(i_tick-1) < value
                                && wire_from_hit->Signal().at(i_tick+1) < value){
                            //std::cout << "found peak!" << std::endl;

                            peak_counter++;
                        }

                    }

                }

                if (peak_counter != 1){
                    //std::cout << "peak counter: " << peak_counter << ", skipping channel" << std::endl;
                    continue;
                }

                // get peak bin tick after t0 correction
                int maximum_tick;
                if (use_t0tagged_tracks) 
                    maximum_tick = h_wire_in_window->GetMaximumBin()
                        + tick_window_left - t0_tick;
                else 
                    maximum_tick = h_wire_in_window->GetMaximumBin() + tick_window_left;

                // now the magic: 
                // loop over the drift bins and check to see if the 
                // t0-corrected pulse center is in each bin
                // if it is then sum the pulses

                for (int bin_it = 0; bin_it < number_drift_bins; bin_it++){

                    //std::cout << "maximum tick: " << maximum_tick << std::endl;
                    //std::cout << "Tick low: " << waveform_intime_start + bin_it*number_ticks_per_bin << std::endl;
                    //std::cout << "Tick high: " << waveform_intime_start + (bin_it+1)*number_ticks_per_bin << std::endl;
                    // set binning for histograms
                    h_wire_baseline_corrected->SetBins(h_wire_in_window->GetNbinsX(),
                            waveform_intime_start + (bin_it) * number_ticks_per_bin, // 800 + bin_it*184
                            waveform_intime_start + (bin_it + 1) * number_ticks_per_bin);

                    if (maximum_tick >= (waveform_intime_start + bin_it * number_ticks_per_bin)
                            && maximum_tick < (waveform_intime_start + (bin_it + 1) * number_ticks_per_bin)) {

                        bin_no = bin_it;
                        //std::cout << "bin_no: " << bin_no << std::endl;

                        h_wire_in_window->GetXaxis()->SetLimits(bin_it * number_ticks_per_bin, (bin_it +1) * number_ticks_per_bin);

                        // apply baseline correction
                        h_wire_baseline_corrected = 
                            _waveform_func.applyGlobalBaselineCorrection(
                                    h_wire_in_window, 
                                    h_wire_baseline_corrected); 

                        // calculate sigma
                        pulse_height = h_wire_baseline_corrected->GetMaximum();
                        mean         = _waveform_func.getSigma(h_wire_baseline_corrected).at(0);
                        sigma        = _waveform_func.getSigma(h_wire_baseline_corrected).at(1);
                        fit_chisq    = _waveform_func.getSigma(h_wire_baseline_corrected).at(2);

                        h_sigma->Fill(sigma);
                        //std::cout << "Baseline sigma: " << sigma << std::endl;
                        h_pulseHeight->Fill(pulse_height);
                        h_nWvfmsInBin->Fill(bin_it, 1);

                        h_driftVsigma->Fill(bin_no*10, sigma);
                        h_driftVPulseHeight->Fill(bin_it, pulse_height);

                        // now find the correction needed to minimise the rms of the sum of the 
                        // histograms
                        if (h_summed_wire_info_per_bin.at(bin_it)->Integral() == 0) 
                            waveform_x_correction = 0;
                        else
                            waveform_x_correction = 
                                _waveform_func.findXCorrection(
                                        _waveform_func, 
                                        h_summed_wire_info_per_bin.at(bin_it), 
                                        h_wire_baseline_corrected, 
                                        number_ticks_per_bin, mean);

                        

                        TH1D* h_waveform_x_correction = 
                            new TH1D("h_waveform_x_correction", 
                                    "", 
                                    number_ticks_per_bin, 
                                    h_wire_baseline_corrected->GetXaxis()->GetXmin(), 
                                    h_wire_baseline_corrected->GetXaxis()->GetXmax()); 

                        for (int ntick = 1; ntick <= h_wire_baseline_corrected->GetNbinsX(); ntick++)
                            h_waveform_x_correction->SetBinContent(
                                    ntick, 
                                    h_wire_baseline_corrected->GetBinContent(ntick+waveform_x_correction));

                        difftree->Fill();
                        //double sigma2 = _waveform_func.getSigma(h_waveform_x_correction).at(1);
                        //std::cout << "Corrected sigma: " << sigma2 << std::endl;

                        // finally add to output histograms
                        h_summed_wire_info_per_bin.at(bin_it)->Add(h_waveform_x_correction);
                        //std::cout << "Summed sigma: " << _waveform_func.getSigma(h_summed_wire_info_per_bin.at(bin_it)).at(1) << std::endl;
                    }
                }
            }
        }
    }
}

void diffmod::LArDiffusion::beginJob()
{

    difftree = tfs->make<TTree>("difftree", "diffusion tree");
    difftree->Branch("hit_peak_time", &hit_peak_time);
    difftree->Branch("hit_peak_time", &hit_peak_time);
    difftree->Branch("t0", &t0);
    difftree->Branch("t0_tick", &t0_tick);
    difftree->Branch("track_x_correction", &track_x_correction);
    difftree->Branch("pulse_height", &pulse_height);
    difftree->Branch("mean", &mean);        
    difftree->Branch("sigma", &sigma);       
    difftree->Branch("fit_chisq", &fit_chisq);  
    difftree->Branch("waveform_x_correction", &waveform_x_correction);
    difftree->Branch("bin_no", &bin_no);
    
    // Troubleshooting histograms
    //TH1I *h_nTrack = tfs->make<TH1I>("h_nTracks", ";No. Tracks/Event;", 100, 0, 100);
    h_trackLength = tfs->make<TH1D>("h_trackLength", ";Track Length (cm);", 25, 0, 256);
    h_cosTheta = tfs->make<TH1D>("h_cosTheta", ";Track cos(#theta);", 50, -1, 1); 
    h_startX = tfs->make<TH1D>("h_startX", ";Starting x-Position (cm);", 36, -50, 310);
    h_sigma = tfs->make<TH1D>("h_sigma", ";#sigma^{2};", 100, 0, 10);
    h_pulseHeight = tfs->make<TH1D>("h_pulseHeight", ";Pulse Height;", 50, 0, 100);
    h_nWvfmsInBin = tfs->make<TH1D>("h_nWvfmsInBin", ";Drift bin; No. Waveforms;", 25, 0, 25);
    
    h_driftVsigma = tfs->make<TH2D>("h_driftVsigma", ";Drift Bin; #sigma;", 25, 0, 256, 100, 0, 20);
    h_driftVPulseHeight = tfs->make<TH2D>("h_driftVPulseHeight", ";Drift Distance (cm); Pulse Height;", 25, 0, 25, 50, 0, 100);

    for (int i = 0; i < number_drift_bins; i++){

        TString histo_name = Form("histo_bin_%i", i);
        h_summed_wire_info_per_bin.push_back(tfs->make<TH1D>(histo_name, "", number_ticks_per_bin, waveform_intime_start + (i * number_ticks_per_bin), waveform_intime_start + ((i+1) * number_ticks_per_bin)));

    }

}

DEFINE_ART_MODULE(diffmod::LArDiffusion)
