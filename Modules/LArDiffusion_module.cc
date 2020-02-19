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
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larana/TruncatedMean/Algorithm/TruncMean.h"
#include "larcore/Geometry/Geometry.h"

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
#include "TVector3.h"

// local
#include "ubana/DiffusionModule/Algorithms/WaveformFunctions.h"
#include "ubana/UBXSec/Algorithms/FiducialVolume.h"        


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

    // prints a histogram for troublshooting
    void printHistogram(TH1D* h);

  private:

    art::ServiceHandle< art::TFileService > tfs;
    art::ServiceHandle< geo::Geometry > geo;

    TTree* difftree;

    // Variables for tree/histograms
    int run = -9999;
    int sub_run = -9999;
    int event = -9999;
    int is_real_data;
    double maximum_tick             = -9999;
    double track_length             = -9999;
    double track_direction_rms      = -9999;
    double cos_theta                = -9999;
    double theta_xz                 = -9999;
    double theta_yz                 = -9999;
    double start_x                  = -9999;
    double start_x_t0corr           = -9999;
    double start_y                  = -9999;
    double start_z                  = -9999;
    double end_x                    = -9999;
    double end_x_t0corr             = -9999;
    double end_y                    = -9999;
    double end_z                    = -9999;
    double hit_peak_time            = -9999;
    double hit_peak_time_t0corr     = -9999;
    double hit_peak_time_stddev     = -9999;
    double hit_rms                  = -9999;
    double hit_charge               = -9999;
    int hit_multiplicity            = -9999;
    int hit_view                    = -9999;
    double sp_x                     = -9999;
    double sp_y                     = -9999;
    double sp_z                     = -9999;
    double sp_x_t0                  = -9999;
    double t0                       = -9999;
    double t0_tick_shift            = -9999;
    double t0_x_shift               = -9999;
    double track_t_correction       = -9999;
    double pulse_height             = -9999;
    double fit_mean                 = -9999;
    double fit_sigma                = -9999;
    double fit_chisq                = -9999;
    double waveform_tick_correction = -9999;
    int bin_no                      = -9999;
    int num_waveforms               = -9999;

    TVector3 track_start;
    TVector3 track_end;
    // Truncated mean has to be a float
    float trunc_mean;

    // For calculations 
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
    std::string hit_sp_assn;
    std::string track_t0_assn;
    std::string sigma_map_file_path;
    std::string sigma_map_directory_file;
    const char* uboonedata_env;
    bool use_t0tagged_tracks;
    bool make_sigma_map;
    double track_dir_rms_cut;
    double sigma_cut;
    double pulse_height_cut;
    float drift_velocity;
    double hit_GOF_cut;
    int hit_multiplicity_cut;
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

    TH1D *h_nWvfmsInBin;

    //TH1D *h_single_waveform;

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
    ubana::FiducialVolume _fiducialVolume;

};


diffmod::LArDiffusion::LArDiffusion(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{



  // defaults set to be MCC9 defaults
  track_label           = p.get< std::string >("TrackLabel"          , "pandora");
  wire_label            = p.get< std::string >("WireLabel"           , "butcher");
  hit_label             = p.get< std::string >("HitLabel"            , "gaushit");
  track_hit_assn        = p.get< std::string >("TrackHitAssn"        , "pandora");
  hit_sp_assn           = p.get< std::string >("HitSpacePointAssn"   , "pandora");
  track_t0_assn         = p.get< std::string >("TrackT0Assn"         , "t0reco");

  sigma_map_file_path      = p.get< std::string >("SigmaMapFilePath"     , "");
  sigma_map_directory_file = p.get< std::string >("SigmaMapDirectoryFile", "DiffusionModule");

  drift_velocity        = p.get< float        >("DriftVelocity"       , 0.1098);
  use_t0tagged_tracks   = p.get< bool         >("UseT0TaggedTracks"   , true);
  make_sigma_map        = p.get< bool         >("MakeSigmaMap"        , false);
  track_dir_rms_cut     = p.get< double       >("TrackDirRMSCut"      , 1e10);
  sigma_cut             = p.get< double       >("SigmaCut"            , 1.0);
  pulse_height_cut      = p.get< double       >("PulseHeightCut"      , 100.0);
  hit_GOF_cut           = p.get< double       >("HitGOFCut"           , 1.1);
  peak_finder_threshold = p.get< float        >("PeakFinderThreshold" , 3.0);
  hit_min_channel       = p.get< unsigned int >("HitMinChannel"       , 6150);
  hit_multiplicity_cut  = p.get< int          >("HitMultiplicityCut"  , 1);
  waveform_size         = p.get< int          >("WaveformSize"        , 6400);
  waveform_intime_start = p.get< int          >("WaveformIntimeStart" , 800);
  waveform_intime_end   = p.get< int          >("WaveformIntimeEnd"   , 5400);
  number_time_bins      = p.get< int          >("NumberTimeBins"      , 25);
  number_dropped_ticks  = p.get< int          >("NumberDroppedTicks"  , 2400);
  waveform_drift_size   = waveform_intime_end - waveform_intime_start; // 4600
  number_ticks_per_bin  = waveform_drift_size/number_time_bins; // 184

  MF_LOG_VERBATIM("LArDiffusion") 
    << "PRINTING FHICL CONFIGURATION"
    << "\n-- track_label           : " << track_label
    << "\n-- wire_label            : " << wire_label
    << "\n-- hit_label             : " << hit_label
    << "\n-- track_hit_assn        : " << track_hit_assn
    << "\n-- hit_sp_assn           : " << hit_sp_assn
    << "\n-- track_t0_assn         : " << track_t0_assn
    << "\n-- drift_velocity        : " << drift_velocity
    << "\n-- use_t0tagged_tracks   : " << use_t0tagged_tracks
    << "\n-- make_sigma_map        : " << make_sigma_map
    << "\n-- track_dir_rms_cut     : " << track_dir_rms_cut
    << "\n-- sigma_cut             : " << sigma_cut
    << "\n-- pulse_height_cut      : " << pulse_height_cut
    << "\n-- hit_GOF_cut           : " << hit_GOF_cut
    << "\n-- peak_finder_threshold : " << peak_finder_threshold
    << "\n-- hit_min_channel       : " << hit_min_channel
    << "\n-- hit_multiplicity_cut  : " << hit_multiplicity_cut
    << "\n-- waveform_size         : " << waveform_size
    << "\n-- waveform_intime_start : " << waveform_intime_start
    << "\n-- waveform_intime_end   : " << waveform_intime_end
    << "\n-- number_time_bins      : " << number_time_bins
    << "\n-- number_dropped_ticks  : " << number_dropped_ticks
    << "\n-- waveform_drift_size   : " << waveform_drift_size
    << "\n-- number_ticks_per_bin  : " << number_ticks_per_bin;

  // define fiducial volume for analysis
  fhicl::ParameterSet const p_fv     = p.get<fhicl::ParameterSet>("FiducialVolume");
  _fiducialVolume.Configure(p_fv,
      geo->DetHalfHeight(),
      2.*geo->DetHalfWidth(),
      geo->DetLength());

  _fiducialVolume.PrintConfig();

}

void diffmod::LArDiffusion::analyze(art::Event const & e) {
  run          = e.run();
  sub_run      = e.subRun();
  event        = e.event();
  is_real_data = e.isRealData();

  std::cout << "[DIFFMOD] --- Processing event " 
    << run << "." << sub_run << "." << event << std::endl;

  // Tracks
  art::Handle< std::vector<recob::Track> > track_handle;
  e.getByLabel(track_label, track_handle);
  std::vector< art::Ptr<recob::Track> > track_ptr_vector;
  art::fill_ptr_vector(track_ptr_vector, track_handle);

  // Wires
  art::Handle< std::vector<recob::Wire> > wire_handle;
  e.getByLabel(wire_label, wire_handle);
  std::vector< art::Ptr<recob::Wire> > wire_ptr_vector;
  art::fill_ptr_vector(wire_ptr_vector, wire_handle);

  // Hits 
  art::Handle< std::vector<recob::Hit> > hit_handle;
  e.getByLabel(hit_label, hit_handle);

  // Associations
  art::FindManyP< recob::Hit >        hits_from_tracks(track_handle, e, track_hit_assn);
  art::FindManyP< recob::SpacePoint > sp_from_hits    (hit_handle  , e, hit_sp_assn);
  art::FindManyP< anab::T0 >          t0_from_tracks  (track_handle, e, track_t0_assn);

  // loop tracks, get associated hits
  for (size_t i_tr = 0; i_tr < track_ptr_vector.size(); i_tr++){

    art::Ptr< recob::Track > thisTrack = track_ptr_vector.at(i_tr);

    bool printHitMsg = true;
    std::vector< art::Ptr< anab::T0 > > t0_from_track;
    if (use_t0tagged_tracks) {
      t0_from_track = t0_from_tracks.at(thisTrack.key());

      if (t0_from_track.size() != 1) {
        continue;
      }
      else if (t0_from_track.size() == 1) {
        art::Ptr< anab::T0 > thisT0 = t0_from_track.at(0);
        t0 = thisT0->Time(); // us
      }

    }

    else {
      // TODO: Make sure this matches the t0 in the generator fcl file if using single muons 
      t0 = 800.; // ticks
    }

    // multiply t0 by drift velocity to get the x shift value
    t0_x_shift = t0 * drift_velocity; 
    t0_tick_shift = t0 * 2; // convert to ticks

    ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::GlobalCoordinateSystemTag > trkDir = thisTrack->StartDirection();
    theta_xz = std::abs(std::atan2(trkDir.X(), trkDir.Z()))* 180 / 3.14159;
    theta_yz = std::abs(std::atan2(trkDir.Y(), trkDir.Z()))* 180 / 3.14159;

    track_length        = thisTrack->Length();
    cos_theta           = thisTrack->Theta();
    track_start         = thisTrack->Start<TVector3>();
    track_end           = thisTrack->End<TVector3>();
    start_x             = track_start.X();
    start_x_t0corr      = track_start.X() - t0_x_shift;
    start_y             = track_start.Y();
    start_z             = track_start.Z();
    end_x               = track_end.X();
    end_x_t0corr        = track_end.X() - t0_x_shift;
    end_y               = track_end.Y();
    end_z               = track_end.Z();

    std::vector<double> dotProds;
    dotProds.resize(0);
    for (size_t i_pt = 1; i_pt < thisTrack->NumberTrajectoryPoints(); i_pt++){
      recob::Track::Vector_t firstDir = thisTrack->DirectionAtPoint(0);
      recob::Track::Vector_t thisDir = thisTrack->DirectionAtPoint(i_pt);

      double dotProd = firstDir.Dot(thisDir);
      std::cout << "dot product is... " << dotProd << std::endl;
      dotProds.push_back(dotProd);
    }
    track_direction_rms = TMath::RMS(dotProds.size(), &dotProds[0]);
    std::cout << "rms is ... " << track_direction_rms << std::endl;
    if (track_direction_rms > track_dir_rms_cut) {
      std::cout << "[DIFFMOD]: Track RMS of " << track_direction_rms << 
        " falls above cut value of " << track_dir_rms_cut << 
        ". Skipping it" << std::endl;
    }

    std::vector< art::Ptr< recob::Hit > > hits_from_track = hits_from_tracks.at(thisTrack.key());

    // loop hits
    for (size_t i_hit = 0; i_hit < hits_from_track.size(); i_hit++){

      art::Ptr< recob::Hit > thisHit = hits_from_track.at(i_hit);

      std::vector< art::Ptr< recob::SpacePoint > > sps_from_hit = sp_from_hits.at(thisHit.key());
      if (sps_from_hit.size() != 1){
        MF_LOG_VERBATIM("LArDiffusion")
          << sps_from_hit.size() 
          << " spacepoints from the hit, just taking first one";
      }
      if (sps_from_hit.size() == 0){
        MF_LOG_VERBATIM("LArDiffusion")
          << "zero spacepoints associated with hit, skip this hit";
          continue;
      }

      art::Ptr< recob::SpacePoint > thisSpacePoint = sps_from_hit.at(0);

      const double* spXYZ = thisSpacePoint->XYZ();

      sp_x = spXYZ[0];
      sp_x_t0 = sp_x - t0_x_shift;
      sp_y = spXYZ[1];
      sp_z = spXYZ[2];

      bool  isInFV = _fiducialVolume.InFV(sp_x_t0,
                                          sp_y,
                                          sp_z);

      // if hit selection is not passed then ignore the hit
      if (!isInFV) continue;
      if (!_waveform_func.passesHitSelection(thisHit, 
                                             hit_GOF_cut, 
                                             hit_multiplicity_cut)) continue;

      // get wire information for hit
      art::Ptr< recob::Wire > wire_from_hit;

      for (size_t i_w = 0; i_w < wire_ptr_vector.size(); i_w++) {

        if ( wire_ptr_vector.at(i_w)->Channel() == thisHit->Channel())
          wire_from_hit = wire_ptr_vector.at(i_w);

      }

      hit_peak_time        = thisHit->PeakTime();
      hit_peak_time_t0corr = thisHit->PeakTime() - t0_tick_shift;

      if (hit_peak_time_t0corr>3000 && hit_peak_time_t0corr<4500 && printHitMsg) {
        std::cout << "[BADEVENT] Weird hit_rms value in run " << run << " subrun " << 
          sub_run << " event " << event << std::endl;
        std::cout << "[BADEVENT] Channel " << thisHit->Channel() << std::endl;
        printHitMsg = false;
      }

      hit_peak_time_stddev = thisHit->SigmaPeakTime();
      hit_rms              = thisHit->RMS();
      hit_charge           = thisHit->Integral();
      hit_view             = thisHit->View();
      hit_multiplicity     = thisHit->Multiplicity();

      tick_window_size  = number_ticks_per_bin;
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
          + tick_window_left - t0_tick_shift;
      }
      else {
        maximum_tick = h_wire_in_window->GetMaximumBin() + tick_window_left;
      }

      /*
         MF_LOG_VERBATIM("LArDiffusion")
         << "PRINTING INFORMATION FOR HIT " << i_hit
         << "\n-- hit_peak_time        : "       << hit_peak_time
         << "\n-- t0_tick_shift        : "       << t0_tick_shift
         << "\n-- hit_peak_time_t0corr : "       << hit_peak_time_t0corr
         << "\n-- tick_window_size     : "       << tick_window_size
         << "\n-- tick_window_left     : "       << tick_window_left
         << "\n-- tick_window_right    : "       << tick_window_right
         << "\n-- maximum_tick         : "       << maximum_tick; 
         */

      // now the magic: 
      // loop over the drift bins and check to see if the 
      // t0-corrected pulse center is in each bin
      // if it is then sum the pulses

      for (int bin_it = 0; bin_it < number_time_bins; bin_it++){

        // get bin edges
        double binEdgeLeft  = waveform_intime_start + ((bin_it)   * number_ticks_per_bin);
        double binEdgeRight = waveform_intime_start + ((bin_it+1) * number_ticks_per_bin);

        // Set binning for histograms
        h_wire_baseline_corrected->SetBins(
            h_wire_in_window->GetNbinsX(),
            binEdgeLeft, 
            binEdgeRight);

        // if maximum tick of the histogram is in this deift bin, grab it! 
        if (maximum_tick >= binEdgeLeft && maximum_tick < binEdgeRight) {

          bin_no = bin_it;

          /*
             MF_LOG_VERBATIM("LArDiffusion")
             << "-- Falls into bin "    << bin_no
             << "\n-- histogram set to"
             << "\n---- binEdgeLeft: "  << binEdgeLeft
             << "\n---- binEdgeRight: " << binEdgeRight;
             */

          //h_wire_in_window->GetXaxis()->SetLimits(bin_it * number_ticks_per_bin, (bin_it +1) * number_ticks_per_bin);

          // apply baseline correction
          h_wire_baseline_corrected = 
            _waveform_func.applyGlobalBaselineCorrection(
                h_wire_in_window, 
                h_wire_baseline_corrected); 

          // Save individual waveform for plotting
          //  TODO: Why is this segfaulting?
          //if (run == 1 && sub_run == 785 && event == 35281) {
          //  double lowVal = 0., highVal = 0.;
          //  lowVal = h_wire_in_window->GetBinLowEdge(1);
          //  std::cout << "Got low val " << lowVal << std::endl;
          //  highVal = h_wire_in_window->GetBinLowEdge(h_wire_in_window->GetNbinsX() );
          //  std::cout << "Got high val " << highVal << std::endl;
          //  if (!h_single_waveform) {
          //    std::cout << "Bad single waveform hist" << std::endl;
          //    continue;
          //  }
          //  h_single_waveform->GetXaxis()->SetLimits(lowVal, highVal);
          //  std::cout << "Set limits" << std::endl;
          //  //h_single_waveform->GetXaxis()->SetLimits(4300, 4500);
          //  h_single_waveform->GetYaxis()->SetRangeUser(-1, 8);
          //  std::cout << "Set range" << std::endl;
          //  for (int i_wv = 1; i_wv < 101; i_wv++) {
          //    h_single_waveform->SetBinContent(i_wv, h_wire_in_window->GetBinContent(i_wv+50) );
          //    //std::cout << "wvfm bin " << i_wv << " has " << h_single_waveform->GetBinContent(i_wv) << std::endl;
          //    //std::cout << "Should be " << h_wire_in_window->GetBinContent(i_wv) << std::endl;
          //  }
          //}   


          // calculate sigma
          pulse_height = h_wire_baseline_corrected->GetMaximum();
          fit_mean     = _waveform_func.getSigma(h_wire_baseline_corrected).at(0);
          fit_sigma    = _waveform_func.getSigma(h_wire_baseline_corrected).at(1);
          fit_chisq    = _waveform_func.getSigma(h_wire_baseline_corrected).at(2);

          difftree->Fill();

          if (make_sigma_map) {
            h_sigma_hists.at(bin_no)->Fill(fit_sigma);
            h_sigma_v_bin_precut->Fill(bin_no, fit_sigma);
            h_pulse_height_hists.at(bin_no)->Fill(pulse_height);
            h_pulse_height_v_bin_precut->Fill(bin_no, pulse_height);
            h_sigma_v_pulse_height_precut->Fill(fit_sigma, pulse_height);
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
            h_sigma_v_bin_precut->Fill(bin_no, fit_sigma);
            h_pulse_height_v_bin_precut->Fill(bin_no, pulse_height);
            h_sigma_v_pulse_height_precut->Fill(fit_sigma, pulse_height);

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

            // Maximum 
            double sigma_lowerLimit = 
              sigmaMaxs.at(bin_it) - sigma_cut * h_sigma_hists.at(bin_no)->GetStdDev();
            double sigma_higherLimit = 
              sigmaMaxs.at(bin_it) + sigma_cut * h_sigma_hists.at(bin_no)->GetStdDev();
            double pulseHeight_lowerLimit = 
              pulseHeightMaxs.at(bin_it) - pulse_height_cut * h_pulse_height_hists.at(bin_no)->GetStdDev();
            double pulseHeight_higherLimit = 
              pulseHeightMaxs.at(bin_it) + pulse_height_cut * h_pulse_height_hists.at(bin_no)->GetStdDev();


            if (fit_sigma < sigma_lowerLimit 
                || fit_sigma > sigma_higherLimit 
                || pulse_height < pulseHeight_lowerLimit 
                || pulse_height > pulseHeight_higherLimit) {
              //std::cout << "Bad sigma value" << std::endl;
              continue;
            }

            h_sigma_hists.at(bin_no)->Fill(fit_sigma);
            h_sigma_v_bin_postcut->Fill(bin_no, fit_sigma);
            h_pulse_height_hists.at(bin_no)->Fill(pulse_height);
            h_pulse_height_v_bin_postcut->Fill(bin_no, pulse_height);
            h_sigma_v_pulse_height_postcut->Fill(pulse_height, fit_sigma);
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
            if (h_summed_wire_info_per_bin.at(bin_it)->Integral() == 0){
              waveform_tick_correction = 0;
            }
            else {
              waveform_tick_correction = 
                _waveform_func.findXCorrection(
                    h_summed_wire_info_per_bin.at(bin_it), 
                    h_wire_baseline_corrected, 
                    number_ticks_per_bin, 
                    fit_mean);
            }

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

            // finally add to output histograms
            h_summed_wire_info_per_bin.at(bin_it)->Add(h_waveform_tick_correction);
            //std::cout << "[DIFFMOD]: Summed sigma: " << _waveform_func.getSigma(h_summed_wire_info_per_bin.at(bin_it)).at(1) << std::endl;

          } // !make_sigma_map
        } // if maximum tick in bin
      } // loop bins
    } // loop hits
  } // loop tracks
} // LArDiffusion::analyze

void diffmod::LArDiffusion::beginJob()
{

  h_sigma_v_bin_precut = tfs->make<TH2D>(
      "h_sigma_v_bin_precut", 
      ";Bin no. ; #sigma_{t}^{2} (#mus^{2});", 
      number_time_bins, 0, number_time_bins, 
      100, 0, 10);

  h_sigma_v_bin_postcut = tfs->make<TH2D>(
      "h_sigma_v_bin_postcut", 
      ";Bin no. ; #sigma_{t}^{2} (#mus^{2});", 
      number_time_bins, 0, number_time_bins, 
      100, 0, 10);

  h_pulse_height_v_bin_precut = tfs->make<TH2D>(
      "h_pulse_height_v_bin_precut", 
      ";Bin no. ; Pulse Height (Arb. Units);", 
      number_time_bins, 0, number_time_bins, 
      100, 0, 20);

  h_pulse_height_v_bin_postcut = tfs->make<TH2D>(
      "h_pulse_height_v_bin_postcut", 
      ";Bin no. ; Pulse Height (Arb. Units);", 
      number_time_bins, 0, number_time_bins, 
      100, 0, 20);

  h_sigma_v_pulse_height_precut = tfs->make<TH2D>(
      "h_sigma_v_pulse_height_precut", 
      ";#sigma_{t}^{2} (#mus^{2}); Pulse Height (Arb. Units);", 
      100, 0, 10, 
      100, 0, 20);

  h_sigma_v_pulse_height_postcut = tfs->make<TH2D>(
      "h_sigma_v_pulse_height_postcut", 
      ";#sigma_{t}^{2} (#mus^{2}); Pulse Height (Arb. Units);", 
      100, 0, 10, 
      100, 0, 20);

  h_theta_xz_v_bin = tfs->make<TH2D>(
      "h_theta_xz_v_bin", 
      ";Bin no. ; #theta_{xz} (Deg.);", 
      number_time_bins, 0, number_time_bins, 
      100, 0, 20);

  h_theta_yz_v_bin = tfs->make<TH2D>(
      "h_theta_yz_v_bin", 
      ";Bin no. ; #theta_{yz} (Deg.);", 
      number_time_bins, 0, number_time_bins, 
      250, 0, 50);

  h_sigma_hist_medians = tfs->make<TH1D>(
      "h_sigma_hist_medians", 
      ";Median #sigma per bin;", 
      number_time_bins, 0, number_time_bins);

  h_sigma_hist_maxs = tfs->make<TH1D>(
      "h_sigma_hist_maxs", 
      ";Max #sigma per bin;", 
      number_time_bins, 0, number_time_bins);

  difftree = tfs->make<TTree>("difftree"          , "diffusion tree");
  difftree->Branch("run"                          , &run);
  difftree->Branch("sub_run"                      , &sub_run);
  difftree->Branch("event"                        , &event);
  difftree->Branch("maximum_tick"                 , &maximum_tick);
  difftree->Branch("track_length"                 , &track_length);
  difftree->Branch("track_direction_rms"          , &track_direction_rms);
  difftree->Branch("cos_theta"                    , &cos_theta);
  difftree->Branch("theta_xz"                     , &theta_xz);
  difftree->Branch("theta_yz"                     , &theta_yz);
  difftree->Branch("start_x"                      , &start_x);
  difftree->Branch("start_x_t0corr"               , &start_x_t0corr);
  difftree->Branch("start_y"                      , &start_y);
  difftree->Branch("start_z"                      , &start_z);
  difftree->Branch("end_x"                        , &end_x);
  difftree->Branch("end_x_t0corr"                 , &end_x_t0corr);
  difftree->Branch("end_y"                        , &end_y);
  difftree->Branch("end_z"                        , &end_z);
  difftree->Branch("hit_peak_time"                , &hit_peak_time);
  difftree->Branch("hit_peak_time_t0corr"         , &hit_peak_time_t0corr);
  difftree->Branch("hit_peak_time_stddev"         , &hit_peak_time_stddev);
  difftree->Branch("hit_rms"                      , &hit_rms);
  difftree->Branch("hit_charge"                   , &hit_charge);
  difftree->Branch("hit_multiplicity"             , &hit_multiplicity);
  difftree->Branch("hit_view"                     , &hit_view);
  difftree->Branch("sp_x"                         , &sp_x);
  difftree->Branch("sp_x_t0"                      , &sp_x_t0);
  difftree->Branch("sp_y"                         , &sp_y);
  difftree->Branch("sp_z"                         , &sp_z);
  difftree->Branch("t0"                           , &t0);
  difftree->Branch("t0_x_shift"                   , &t0_x_shift);
  difftree->Branch("pulse_height"                 , &pulse_height);
  difftree->Branch("fit_mean"                     , &fit_mean);
  difftree->Branch("fit_sigma"                    , &fit_sigma);
  difftree->Branch("fit_chisq"                    , &fit_chisq);
  difftree->Branch("waveform_tick_correction"     , &waveform_tick_correction);
  difftree->Branch("bin_no"                       , &bin_no);
  difftree->Branch("num_waveforms"                , &num_waveforms);

  if (!make_sigma_map) {
    //h_single_waveform = tfs->make<TH1D>("h_single_waveform", ";Time (ticks); Arb. Units;", 100, 0, 100);

    h_nWvfmsInBin = tfs->make<TH1D>("h_nWvfmsInBin", ";Drift bin; No. Waveforms;", 25, 0, 25);

    // Import sigmaMap, assuming it already exists
    char fullPath[200];
    uboonedata_env = getenv("UBOONEDATA_DIR");
    if (uboonedata_env!=NULL) {
      std::cout << "[DIFFMOD]: Got uboonedata env path " << uboonedata_env << std::endl;
    }
    strcpy(fullPath, uboonedata_env);
    strcat(fullPath, sigma_map_file_path.c_str() );

    std::cout << "[DIFFMOD]: Running without producing sigma map. Checking that it exists..." << std::endl;
    std::cout << "[DIFFMOD]: Getting sigma map from file path " << fullPath << std::endl;
    TFile sigmaMap(fullPath, "READ");
    //TFile sigmaMap(sigma_map_file_path.c_str(), "READ");

    if (sigmaMap.IsOpen() == false){
      std::cout << "[DIFFMOD]: No sigma map! Run module using run_diffusion_module_sigmamap.fcl first, " << 
        "or check that you're in the right directory.\n" << std::endl;
    }

    std::cout << "[DIFFMOD]: Got sigma map" << std::endl;
    std::cout << "[DIFFMOD]: Getting sigma and pulse heights hists..." << std::endl;

    for (int i = 0; i < number_time_bins; i++){

      // declare summed waveforms
      TString histo_name = Form("summed_waveform_bin_%i", i);
      h_summed_wire_info_per_bin.push_back(tfs->make<TH1D>(histo_name, "", number_ticks_per_bin, waveform_intime_start + (i * number_ticks_per_bin), waveform_intime_start + ((i+1) * number_ticks_per_bin)));

      sigmaDistsPerBin.resize(0);

      TString sigmaMapFilePath(sigma_map_directory_file);
      TString sigmaMapHistoName = Form("/h_sigma_%i", i); 
      TString t = sigmaMapFilePath+sigmaMapHistoName;
      h_sigma_hists.push_back((TH1D*)sigmaMap.Get(t.Data()));

      TString pulseHeightHistoName = Form("/h_pulse_height_%i", i); 
      TString t2 = sigmaMapFilePath+pulseHeightHistoName;
      h_pulse_height_hists.push_back((TH1D*)sigmaMap.Get(t2.Data()));

      // Three options for picking out waveforms, either
      // 1) get waveforms around the median
      // 2) get the peak bin value
      // 3) calculate the truncated mean

      // OPTION 1) Calculate medians in each bin
      sigmaMedians.push_back(_waveform_func.getMedian(h_sigma_hists.at(i) ) );
      pulseHeightMedians.push_back(_waveform_func.getMedian(h_pulse_height_hists.at(i) ) );

      // OPTION 2) Calculate maximum in each bin
      int sigmaMaxBin = h_sigma_hists.at(i)->GetMaximumBin();
      int pulseHeightMaxBin = h_pulse_height_hists.at(i)->GetMaximumBin();
      sigmaMaxs.push_back(h_sigma_hists.at(i)->GetXaxis()->GetBinCenter(sigmaMaxBin) );
      pulseHeightMaxs.push_back(h_pulse_height_hists.at(i)->GetXaxis()->GetBinCenter(pulseHeightMaxBin) );

      // OPTION 3) Take sigma hist and calculate truncated mean 
      for (int j = 1; j < h_sigma_hists.at(i)->GetNbinsX()+1; j++) {
        if (h_sigma_hists.at(i)->GetBinContent(j) > 0) {
          for (int k = 0; k < h_sigma_hists.at(i)->GetBinContent(j); k++ ) {
            sigmaDistsPerBin.push_back(h_sigma_hists.at(i)->GetXaxis()->GetBinCenter(j) );
          }
        }
      }
      std::sort(sigmaDistsPerBin.begin(), sigmaDistsPerBin.end() );
      if (sigmaDistsPerBin.size() > 0){
        trunc_mean = _trunc_mean_func.CalcIterativeTruncMean(
            sigmaDistsPerBin,    // v
            0,                   // nmin
            1000,                // nmax
            0,                   // currentiteration
            0,                   // lmin
            0.02,                // convergence limit
            1);                  // nsigma
      }
      std::cout << sigmaMedians.at(i) << "\t" << sigmaMaxs.at(i) << "\t" << trunc_mean << std::endl;
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

void diffmod::LArDiffusion::printHistogram(TH1D* h){
  for (int i = 0; i < h->GetNbinsX(); i++){
    std::cout << i << " " << h->GetBinLowEdge(i) << " " << h->GetBinContent(i) << std::endl;
  }
}

DEFINE_ART_MODULE(diffmod::LArDiffusion)
