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
#include "larana/TruncatedMean/Algorithm/TruncMean.h"

// art
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Services/Optional/TFileService.h"

// ROOT
#include "TH1.h"
#include "TH2D.h"
#include "TTree.h"
#include "TStyle.h"
#include "TFile.h"
#include "TString.h"

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

        // Variables for tree/histograms
        int run;
        int sub_run;
        int event;
        int is_real_data;

        double maximum_tick;
        double track_length;
        double cos_theta;
        double theta_xz;
        double theta_yz;
        double start_x;
        double hit_peak_time;
        double hit_peak_time_stddev;
        double hit_rms;
        double hit_charge;
        double hit_multiplicity;
        double hit_peak_time_postSel;
        double hit_peak_time_stddev_postSel;
        double hit_rms_postSel;
        double hit_charge_postSel;
        double hit_multiplicity_postSel;
        double t0;
        double t0_tick;
        double track_t_correction;
        double pulse_height;
        double mean;
        double sigma;       
        double fit_chisq;   
        double waveform_tick_correction;
        int bin_no;
        int num_waveforms;
        TVector3 track_start;
        
        // Truncated mean has to be a float
        float trunc_mean;


        // For calculations 
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
        bool make_sigma_map;
        double sigma_cut;
        double pulse_height_cut;
        float drift_velocity;
        double hit_GOF_cut;
        int hit_multiplicity_cut;
        int hit_view;
        int hit_min_channel;
        int waveform_size;
        int waveform_intime_start;
        int waveform_intime_end;
        int number_time_bins;
        float drift_distance;
        float peak_finder_threshold;
        int number_dropped_ticks;

        // manipulation histograms
        // binnings are placeholders and will be modified later

        // histogram opened around the hit peak value 
        TH1D* h_wire_in_window = tfs->make<TH1D>("h_wire_in_window", "", 100, 0, 100);

        // after baseline correcting
        TH1D* h_wire_baseline_corrected = tfs->make<TH1D>("h_wire_baseline_corrected", "", 100, 0, 100);

        //TH1D *h_single_waveform;
        TH1D *h_nWvfmsInBin;
        //TH1D *h_correctedTicks;

        // For dynamic sigma cut
        std::vector<double> sigmaMedians;
        std::vector<double> pulseHeightMedians;
        std::vector<double> sigmaMaxs;
        std::vector<double> pulseHeightMaxs;
      
        // For truncated mean calculation; needs to be float
        std::vector<float> sigmaDistsPerBin;

        // Output histograms
        TH1D *h_sigma_hist_medians;
        TH1D *h_sigma_hist_maxs;
        std::vector<TH1D*> h_summed_wire_info_per_bin; 
        std::vector<TH1D*> h_sigma_hists; 
        std::vector<TH1D*> h_pulse_height_hists; 
        TH2D *h_sigma_v_bin_precut;
        TH2D *h_sigma_v_bin_postcut;
        TH2D *h_pulse_height_v_bin_precut;
        TH2D *h_pulse_height_v_bin_postcut;
        TH2D *h_sigma_v_pulse_height_precut;
        TH2D *h_sigma_v_pulse_height_postcut;
        TH2D *h_theta_xz_v_bin;
        TH2D *h_theta_yz_v_bin;

        // other classes
        diffmod::WaveformFunctions _waveform_func;
        TruncMean _trunc_mean_func;

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
    make_sigma_map        = p.get< bool >("MakeSigmaMap", false);
    sigma_cut             = p.get< double >("SigmaCut", 1.0);
    pulse_height_cut      = p.get< double >("PulseHeightCut", 100.0);
    drift_velocity        = p.get< float >("DriftVelocity", 0.1098);
    hit_GOF_cut           = p.get< double >("HitGOFCut", 1.1);
    hit_multiplicity_cut  = p.get< int >("HitMultiplicityCut", 1);
    hit_view              = p.get< int >("HitView", 1);
    hit_min_channel       = p.get< unsigned int >("HitMinChannel", 6150);
    waveform_size         = p.get< int >("WaveformSize", 6400);
    waveform_intime_start = p.get< int >("WaveformIntimeStart", 800);
    waveform_intime_end   = p.get< int >("WaveformIntimeEnd", 5400);
    number_time_bins      = p.get< int >("NumberTimeBins", 25);
    number_dropped_ticks  = p.get< int >("NumberDroppedTicks", 2400);
    peak_finder_threshold = p.get< float >("PeakFinderThreshold", 3.0);
    waveform_drift_size  = waveform_intime_end - waveform_intime_start; // 4600
    number_ticks_per_bin = waveform_drift_size/number_time_bins; // 184
}

void diffmod::LArDiffusion::analyze(art::Event const & e) {
    run = e.run();
    sub_run = e.subRun();
    event = e.event();
    is_real_data = e.isRealData();

    std::cout << "[DIFFMOD] --- Processing event " 
        << run << "." << sub_run << "." << event << std::endl;

    // Tracks
    art::Handle< std::vector<recob::Track> > track_handle;
    e.getByLabel(track_label, track_handle);
    std::vector< art::Ptr<recob::Track> > track_ptr_vector;
    art::fill_ptr_vector(track_ptr_vector, track_handle);

    // Hits 
    art::Handle< std::vector<recob::Hit> > hit_handle;
    e.getByLabel(hit_label, hit_handle);

    // Associations
    art::FindManyP< recob::Hit > hits_from_tracks(track_handle, e, track_hit_assn);
    art::FindManyP< anab::T0 > t0_from_tracks(track_handle, e, track_t0_assn);
    art::FindManyP< recob::Wire > wire_from_hits(hit_handle, e, hit_wire_assn);

    // loop tracks, get associated hits
    for (size_t i_tr = 0; i_tr < track_ptr_vector.size(); i_tr++){

        art::Ptr< recob::Track > thisTrack = track_ptr_vector.at(i_tr);

        ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::GlobalCoordinateSystemTag > trkDir = thisTrack->StartDirection();
        theta_xz = std::abs(std::atan2(trkDir.X(), trkDir.Z()))* 180 / 3.14159;
        theta_yz = std::abs(std::atan2(trkDir.Y(), trkDir.Z()))* 180 / 3.14159;

        track_length = thisTrack->Length();
        cos_theta = thisTrack->Theta();
        track_start = thisTrack->Start<TVector3>();
        start_x = track_start.X();

        std::vector< art::Ptr< anab::T0 > > t0_from_track;
        if (use_t0tagged_tracks) {
            //std::vector< art::Ptr< anab::T0 > > t0_from_track = t0_from_tracks.at(thisTrack.key());
            t0_from_track = t0_from_tracks.at(thisTrack.key());
        }

        std::vector< art::Ptr< recob::Hit > > hits_from_track = hits_from_tracks.at(thisTrack.key());

        // if a t0 exists (should, if using diffusion-filtered sample), then
        // go grab it and get the tick correction value

        if (t0_from_track.size() == 1 && use_t0tagged_tracks) {
          art::Ptr< anab::T0 > thisT0 = t0_from_track.at(0);
          t0 = thisT0->Time();
          //track_x_correction = t0 * drift_velocity;
        }
        else {
          // TODO check units. Is it 800 ticks or microseconds? Want this in ticks
          t0 = 800.; // 800 microsecond (tick?) default value for single muon samples
        }

        t0_tick = t0 * 2;

        /*
         correction time is the t0 of the track 
         + the number of ticks we drop
         + the tick at which we begin to be in time
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
            if (!_waveform_func.passesHitSelection(thisHit, hit_GOF_cut, hit_multiplicity_cut, hit_view, hit_min_channel)) continue;

            // get wire information for hit
            art::Ptr< recob::Wire > wire_from_hit = wire_from_hits.at(thisHit.key()).at(0);            
            hit_peak_time = thisHit->PeakTime(); 
            hit_peak_time_stddev = thisHit->SigmaPeakTime(); 
            hit_rms = thisHit->RMS();
            hit_charge = thisHit->Integral();
            hit_multiplicity = thisHit->Multiplicity();

            // Hit information after quality cuts
            hit_peak_time_postSel = thisHit->PeakTime(); 
            hit_peak_time_stddev_postSel = thisHit->SigmaPeakTime(); 
            hit_rms_postSel = thisHit->RMS();
            hit_charge_postSel = thisHit->Integral();
            hit_multiplicity_postSel = thisHit->Multiplicity();

            tick_window_size = number_ticks_per_bin;
            tick_window_left  = hit_peak_time - tick_window_size/2;
            tick_window_right = hit_peak_time + tick_window_size/2;

            // make sure that the window stops at the edge of the waveform
            if (tick_window_left < 0){
                tick_window_left = 0;
                tick_window_size = tick_window_right; 
            }
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

                if (value > peak_finder_threshold) {

                    // define peak search region
                    if (tick_window_left == 0 || tick_window_right == (int)wire_from_hit->Signal().size())
                      continue;

                    // make sure we only look for the peak
                    if (wire_from_hit->Signal().at(i_tick-1) < value
                        && wire_from_hit->Signal().at(i_tick+1) < value){
                        peak_counter++;
                    }

                }

            }

            if (peak_counter != 1) continue;

            // get peak bin tick after t0 correction
            if (use_t0tagged_tracks) {
                maximum_tick = h_wire_in_window->GetMaximumBin()
                    + tick_window_left - t0_tick;
            }
            else {
                maximum_tick = h_wire_in_window->GetMaximumBin() + tick_window_left;
            }

            //h_correctedTicks->Fill(maximum_tick);

            // now the magic: 
            // loop over the drift bins and check to see if the 
            // t0-corrected pulse center is in each bin
            // if it is then sum the pulses

            for (int bin_it = 0; bin_it < number_time_bins; bin_it++){

                /*
                std::cout << "[DIFFMOD]: Bin no. " << bin_it << std::endl;
                std::cout << "[DIFFMOD]: Maximum tick: " << maximum_tick << std::endl;
                std::cout << "[DIFFMOD]: Tick low: " << waveform_intime_start + bin_it*number_ticks_per_bin << std::endl;
                std::cout << "[DIFFMOD]: Tick high: " << waveform_intime_start + (bin_it+1)*number_ticks_per_bin << std::endl;
                */

                // Set binning for histograms
                h_wire_baseline_corrected->SetBins(h_wire_in_window->GetNbinsX(),
                        waveform_intime_start + (bin_it) * number_ticks_per_bin, // 800 + bin_it*184
                        waveform_intime_start + (bin_it + 1) * number_ticks_per_bin);

                if (maximum_tick >= (waveform_intime_start + bin_it * number_ticks_per_bin)
                        && maximum_tick < (waveform_intime_start + (bin_it + 1) * number_ticks_per_bin)) {

                    bin_no = bin_it;

                    h_wire_in_window->GetXaxis()->SetLimits(bin_it * number_ticks_per_bin, (bin_it +1) * number_ticks_per_bin);

                    // apply baseline correction
                    h_wire_baseline_corrected = 
                        _waveform_func.applyGlobalBaselineCorrection(
                                h_wire_in_window, 
                                h_wire_baseline_corrected); 

                    // Save individual waveform for plotting
                    /*
                     * TODO: Why is this segfaulting?
                    if (run == 1 && sub_run == 785 && event == 35281) {
                        double lowVal = 0., highVal = 0.;
                        lowVal = h_wire_in_window->GetBinLowEdge(1);
                        std::cout << "Got low val " << lowVal << std::endl;
                        highVal = h_wire_in_window->GetBinLowEdge(h_wire_in_window->GetNbinsX() );
                        std::cout << "Got high val " << highVal << std::endl;
                        if (!h_single_waveform) {
                            std::cout << "Bad single waveform hist" << std::endl;
                            continue;
                        }
                        h_single_waveform->GetXaxis()->SetLimits(lowVal, highVal);
                        std::cout << "Set limits" << std::endl;
                        //h_single_waveform->GetXaxis()->SetLimits(4300, 4500);
                        h_single_waveform->GetYaxis()->SetRangeUser(-1, 8);
                        std::cout << "Set range" << std::endl;
                        for (int i_wv = 1; i_wv < 101; i_wv++) {
                            h_single_waveform->SetBinContent(i_wv, h_wire_in_window->GetBinContent(i_wv+50) );
                            //std::cout << "wvfm bin " << i_wv << " has " << h_single_waveform->GetBinContent(i_wv) << std::endl;
                            //std::cout << "Should be " << h_wire_in_window->GetBinContent(i_wv) << std::endl;
                        }
                    }   
                    */

                    // calculate sigma
                    pulse_height = h_wire_baseline_corrected->GetMaximum();
                    mean         = _waveform_func.getSigma(h_wire_baseline_corrected).at(0);
                    sigma        = _waveform_func.getSigma(h_wire_baseline_corrected).at(1);
                    fit_chisq    = _waveform_func.getSigma(h_wire_baseline_corrected).at(2);

                    if (make_sigma_map) {
                        h_sigma_hists.at(bin_no)->Fill(sigma);
                        h_sigma_v_bin_precut->Fill(bin_no, sigma);
                        h_pulse_height_hists.at(bin_no)->Fill(pulse_height);
                        h_pulse_height_v_bin_precut->Fill(bin_no, pulse_height);
                        h_sigma_v_pulse_height_precut->Fill(sigma, pulse_height);
                        if (theta_xz < 90)
                            h_theta_xz_v_bin->Fill(bin_no, theta_xz);
                        if (theta_xz > 90)
                            h_theta_xz_v_bin->Fill(bin_no, 180-theta_xz);
                        if (theta_yz < 90)
                            h_theta_yz_v_bin->Fill(bin_no, theta_yz);
                        if (theta_yz > 90)
                            h_theta_yz_v_bin->Fill(bin_no, 180-theta_yz);
                    }

                    else {
                        // For comparison
                        h_sigma_v_bin_precut->Fill(bin_no, sigma);
                        h_pulse_height_v_bin_precut->Fill(bin_no, pulse_height);
                        h_sigma_v_pulse_height_precut->Fill(sigma, pulse_height);

                        h_sigma_hist_medians->Fill(sigmaMedians.at(bin_it) );
                        h_sigma_hist_maxs->Fill(sigmaMaxs.at(bin_it) );

                        // Dynamic sigma cut: check if pulseHeight, sigma, 
                        // fall within some region around the median
                        /*
                        double sigma_lowerLimit = 
                        sigmaMedians.at(bin_it) - sigma_cut * h_sigma_hists.at(bin_no)->GetStdDev();
                        double sigma_higherLimit = 
                        sigmaMedians.at(bin_it) + sigma_cut * h_sigma_hists.at(bin_no)->GetStdDev();
                        double pulseHeight_lowerLimit = 
                        pulseHeightMedians.at(bin_it) - pulse_height_cut * h_pulse_height_hists.at(bin_no)->GetStdDev();
                        double pulseHeight_higherLimit = 
                        pulseHeightMedians.at(bin_it) + pulse_height_cut * h_pulse_height_hists.at(bin_no)->GetStdDev();
                        */

                        double sigma_lowerLimit = 
                        sigmaMaxs.at(bin_it) - sigma_cut * h_sigma_hists.at(bin_no)->GetStdDev();
                        double sigma_higherLimit = 
                        sigmaMaxs.at(bin_it) + sigma_cut * h_sigma_hists.at(bin_no)->GetStdDev();
                        double pulseHeight_lowerLimit = 
                        pulseHeightMaxs.at(bin_it) - pulse_height_cut * h_pulse_height_hists.at(bin_no)->GetStdDev();
                        double pulseHeight_higherLimit = 
                        pulseHeightMaxs.at(bin_it) + pulse_height_cut * h_pulse_height_hists.at(bin_no)->GetStdDev();
                        
                        if (sigma < sigma_lowerLimit 
                            || sigma > sigma_higherLimit 
                            || pulse_height < pulseHeight_lowerLimit 
                            || pulse_height > pulseHeight_higherLimit) {
                            //std::cout << "Bad sigma value" << std::endl;
                            continue;
                        }

                        h_sigma_hists.at(bin_no)->Fill(sigma);
                        h_sigma_v_bin_postcut->Fill(bin_no, sigma);
                        h_pulse_height_hists.at(bin_no)->Fill(pulse_height);
                        h_pulse_height_v_bin_postcut->Fill(bin_no, pulse_height);
                        h_sigma_v_pulse_height_postcut->Fill(pulse_height, sigma);
                        if (theta_xz < 90)
                            h_theta_xz_v_bin->Fill(bin_no, theta_xz);
                        if (theta_xz > 90)
                            h_theta_xz_v_bin->Fill(bin_no, 180-theta_xz);
                        if (theta_yz < 90)
                            h_theta_yz_v_bin->Fill(bin_no, theta_yz);
                        if (theta_yz > 90)
                            h_theta_yz_v_bin->Fill(bin_no, 180-theta_yz);

                        h_nWvfmsInBin->Fill(bin_it, 1);

                        // Now find the shift (in ticks) needed to minimise the rms 
                        // of the sum of the histograms
                        if (h_summed_wire_info_per_bin.at(bin_it)->Integral() == 0) 
                            waveform_tick_correction = 0;
                        else
                            waveform_tick_correction = 
                                _waveform_func.findXCorrection(
                                        _waveform_func, 
                                        h_summed_wire_info_per_bin.at(bin_it), 
                                        h_wire_baseline_corrected, 
                                        number_ticks_per_bin, mean);

                        
                        TH1D* h_waveform_tick_correction = 
                            new TH1D("h_waveform_tick_correction", 
                                    "", 
                                    number_ticks_per_bin, 
                                    h_wire_baseline_corrected->GetXaxis()->GetXmin(), 
                                    h_wire_baseline_corrected->GetXaxis()->GetXmax()); 

                        for (int ntick = 1; ntick <= h_wire_baseline_corrected->GetNbinsX(); ntick++)
                            h_waveform_tick_correction->SetBinContent(
                                    ntick, 
                                    h_wire_baseline_corrected->GetBinContent(ntick+waveform_tick_correction));

                        difftree->Fill();

                        // finally add to output histograms
                        h_summed_wire_info_per_bin.at(bin_it)->Add(h_waveform_tick_correction);
                        //std::cout << "[DIFFMOD]: Summed sigma: " << _waveform_func.getSigma(h_summed_wire_info_per_bin.at(bin_it)).at(1) << std::endl;
                    }
                }
            }
        }
    }
}

void diffmod::LArDiffusion::beginJob()
{

        h_sigma_v_bin_precut = tfs->make<TH2D>("h_sigma_v_bin_precut", ";Bin no. ; #sigma_{t}^{2} (#mus^{2});", number_time_bins, 0, number_time_bins, 100, 0, 10);
        h_sigma_v_bin_postcut = tfs->make<TH2D>("h_sigma_v_bin_postcut", ";Bin no. ; #sigma_{t}^{2} (#mus^{2});", number_time_bins, 0, number_time_bins, 100, 0, 10);
        h_pulse_height_v_bin_precut = tfs->make<TH2D>("h_pulse_height_v_bin_precut", ";Bin no. ; Pulse Height (Arb. Units);", number_time_bins, 0, number_time_bins, 100, 0, 20);
        h_pulse_height_v_bin_postcut = tfs->make<TH2D>("h_pulse_height_v_bin_postcut", ";Bin no. ; Pulse Height (Arb. Units);", number_time_bins, 0, number_time_bins, 100, 0, 20);
        h_sigma_v_pulse_height_precut = tfs->make<TH2D>("h_sigma_v_pulse_height_precut", ";#sigma_{t}^{2} (#mus^{2}); Pulse Height (Arb. Units);", 100, 0, 10, 100, 0, 20);
        h_sigma_v_pulse_height_postcut = tfs->make<TH2D>("h_sigma_v_pulse_height_postcut", ";#sigma_{t}^{2} (#mus^{2}); Pulse Height (Arb. Units);", 100, 0, 10, 100, 0, 20);
        h_theta_xz_v_bin = tfs->make<TH2D>("h_theta_xz_v_bin", ";Bin no. ; #theta_{xz} (Deg.);", number_time_bins, 0, number_time_bins, 100, 0, 20);
        h_theta_yz_v_bin = tfs->make<TH2D>("h_theta_yz_v_bin", ";Bin no. ; #theta_{yz} (Deg.);", number_time_bins, 0, number_time_bins, 100, 0, 20);
        h_sigma_hist_medians = tfs->make<TH1D>("h_sigma_hist_medians", ";Median #sigma per bin;", number_time_bins, 0, number_time_bins);
        h_sigma_hist_maxs = tfs->make<TH1D>("h_sigma_hist_maxs", ";Max #sigma per bin;", number_time_bins, 0, number_time_bins);
        
    if (!make_sigma_map) {
        difftree = tfs->make<TTree>("difftree", "diffusion tree");
        difftree->Branch("maximum_tick", &maximum_tick);
        difftree->Branch("track_length", &track_length);
        difftree->Branch("cos_theta", &cos_theta);
        difftree->Branch("theta_xz", &theta_xz);
        difftree->Branch("theta_yz", &theta_yz);
        difftree->Branch("start_x", &start_x);
        difftree->Branch("hit_peak_time", &hit_peak_time);
        difftree->Branch("hit_peak_time_stddev", &hit_peak_time_stddev);
        difftree->Branch("hit_rms", &hit_rms);
        difftree->Branch("hit_charge", &hit_charge);
        difftree->Branch("hit_multiplicity", &hit_multiplicity);
        difftree->Branch("hit_peak_time_postSel", &hit_peak_time_postSel);
        difftree->Branch("hit_peak_time_stddev_postSel", &hit_peak_time_stddev_postSel);
        difftree->Branch("hit_rms_postSel", &hit_rms_postSel);
        difftree->Branch("hit_charge_postSel", &hit_charge_postSel);
        difftree->Branch("hit_multiplicity_postSel", &hit_multiplicity_postSel);
        difftree->Branch("t0", &t0);
        difftree->Branch("t0_tick", &t0_tick);
        difftree->Branch("pulse_height", &pulse_height);
        difftree->Branch("mean", &mean);        
        difftree->Branch("sigma", &sigma);       
        difftree->Branch("fit_chisq", &fit_chisq);  
        difftree->Branch("waveform_tick_correction", &waveform_tick_correction);
        difftree->Branch("bin_no", &bin_no);
        difftree->Branch("num_waveforms", &num_waveforms);

        //h_single_waveform = tfs->make<TH1D>("h_single_waveform", ";Time (ticks); Arb. Units;", 100, 0, 100);
        h_nWvfmsInBin = tfs->make<TH1D>("h_nWvfmsInBin", ";Drift bin; No. Waveforms;", 25, 0, 25);
        //h_correctedTicks = tfs->make<TH1D>("h_correctedTicks", ";Corrected tick value;", 1150, 0, 4600);
        
        // Troubleshooting histograms
        /*
        //TH1I *h_nTrack = tfs->make<TH1I>("h_nTracks", ";No. Tracks/Event;", 100, 0, 100);
        h_trackLength = tfs->make<TH1D>("h_trackLength", ";Track Length (cm);", 25, 0, 256);
        h_cosTheta = tfs->make<TH1D>("h_cosTheta", ";Track cos(#theta);", 50, -1, 1); 
        h_startX = tfs->make<TH1D>("h_startX", ";Starting x-Position (cm);", 36, -50, 310);
        h_sigma = tfs->make<TH1D>("h_sigma", ";#sigma^{2};", 100, 0, 10);
        h_pulseHeight = tfs->make<TH1D>("h_pulseHeight", ";Pulse Height;", 50, 0, 100);
        
        h_driftVsigma = tfs->make<TH2D>("h_driftVsigma", ";Drift Bin; #sigma;", 25, 0, 256, 100, 0, 20);
        h_driftVPulseHeight = tfs->make<TH2D>("h_driftVPulseHeight", ";Drift Distance (cm); Pulse Height;", 25, 0, 25, 50, 0, 100);
        */

        for (int i = 0; i < number_time_bins; i++){
            TString histo_name = Form("summed_waveform_bin_%i", i);
            h_summed_wire_info_per_bin.push_back(tfs->make<TH1D>(histo_name, "", number_ticks_per_bin, waveform_intime_start + (i * number_ticks_per_bin), waveform_intime_start + ((i+1) * number_ticks_per_bin)));

        }
        // Import sigmaMap, assuming it already exists
        if(!make_sigma_map) {

            std::cout << "[DIFFMOD]: Running without producing sigma map. Checking that it exists..." << std::endl;
            std::cout << "[DIFFMOD]: Getting sigma map..." << std::endl;
            TString sigma_map_dir = "";
            TFile sigmaMap(sigma_map_dir+"sigma_map.root", "READ");

            if (sigmaMap.IsOpen() == false){
                std::cout << "[DIFFMOD]: No sigma map! Run module using run_sigma_map.fcl first, " << 
                             "or check that you're in the right directory.\n" << std::endl;
            }

            std::cout << "[DIFFMOD]: Got sigma map" << std::endl;
            std::cout << "[DIFFMOD]: Geting sigma and pulse heights hists..." << std::endl;


            for (int i = 0; i < number_time_bins; i++){

                TString sigmaMapHistoName = Form("h_sigma_%i", i); 
                h_sigma_hists.push_back((TH1D*)sigmaMap.Get("DiffusionModule/"+sigmaMapHistoName) );

                TString pulseHeightHistoName = Form("h_pulse_height_%i", i); 
                h_pulse_height_hists.push_back((TH1D*)sigmaMap.Get("DiffusionModule/"+pulseHeightHistoName) );

                // Calculate medians in each bin
                sigmaMedians.push_back(_waveform_func.getMedian(h_sigma_hists.at(i) ) );
                pulseHeightMedians.push_back(_waveform_func.getMedian(h_pulse_height_hists.at(i) ) );
                
                // Calculate maximum in each bin
                int sigmaMaxBin = h_sigma_hists.at(i)->GetMaximumBin();
                int pulseHeightMaxBin = h_pulse_height_hists.at(i)->GetMaximumBin();
                sigmaMaxs.push_back(h_sigma_hists.at(i)->GetXaxis()->GetBinCenter(sigmaMaxBin) );
                pulseHeightMaxs.push_back(h_pulse_height_hists.at(i)->GetXaxis()->GetBinCenter(pulseHeightMaxBin) );

                // Take sigma hist and calculate truncated mean 
                trunc_mean = 0.;
                for (int j = 1; j < h_sigma_hists.at(i)->GetNbinsX()+1; j++) {
                    if (h_sigma_hists.at(i)->GetBinContent(j) > 0) {
                        /*
                        std::cout << "Filling sigma dist bin " << i << " with " 
                                  << h_sigma_hists.at(i)->GetXaxis()->GetBinCenter(j) 
                                  << " " << h_sigma_hists.at(i)->GetBinContent(j) << " times " << std::endl;
                        */

                        for (int k = 0; k < h_sigma_hists.at(i)->GetBinContent(j); k++ ) {
                                sigmaDistsPerBin.push_back(h_sigma_hists.at(i)->GetXaxis()->GetBinCenter(j) );
                        }
                    }
                }
                std::sort(sigmaDistsPerBin.begin(), sigmaDistsPerBin.end() );
                trunc_mean = _trunc_mean_func.CalcIterativeTruncMean(sigmaDistsPerBin,20,100,0,100,0.02,sigmaMedians.at(i) );

                /*
                std::cout << "sigma max bin = " << sigmaMaxBin << std::endl;
                std::cout << "Median = " << sigmaMedians.at(i) << std::endl;
                std::cout << "Max = " << sigmaMaxs.at(i) << std::endl;
                std::cout << "TruncMean = " << trunc_mean << std::endl;
                std::cout << "-------------------" << std::endl;
                */
                std::cout << sigmaMedians.at(i) << "\t" << sigmaMaxs.at(i) << "\t" << trunc_mean << std::endl;

            }
        }
    }
    else {
        for (int n = 0; n < number_time_bins; n++) {
            TString sigmaHistName = Form("h_sigma_%i", n);
            h_sigma_hists.push_back(tfs->make<TH1D>(sigmaHistName, ";#sigma_{t};", 250, 0, 10) );
            TString pulseHeightHistName = Form("h_pulse_height_%i", n);
            h_pulse_height_hists.push_back(tfs->make<TH1D>(pulseHeightHistName, ";Pulse Height;", 250, 0, 20) );
        }
      
    }

}

DEFINE_ART_MODULE(diffmod::LArDiffusion)
