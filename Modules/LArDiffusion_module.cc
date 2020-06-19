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
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
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
#include "ubana/DiffusionModule/Algorithms/Utilities.h"
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

    // clears vectors
    void clearVectors();

  private:

    art::ServiceHandle< art::TFileService > tfs;
    art::ServiceHandle< geo::Geometry > geo;

    TTree* difftree;

    // Variables for tree/histograms
    int run                                                    = -9999;
    int sub_run                                                = -9999;
    int event                                                  = -9999;
    std::vector< double >*                track_length         = nullptr;
    std::vector< double >*                track_avg_trans_dist = nullptr;
    std::vector< double >*                track_cos_theta      = nullptr;
    std::vector< double >*                track_theta_xz       = nullptr;
    std::vector< double >*                track_theta_yz       = nullptr;
    std::vector< double >*                track_start_x        = nullptr;
    std::vector< double >*                track_start_x_t0corr = nullptr;
    std::vector< double >*                track_start_y        = nullptr;
    std::vector< double >*                track_start_z        = nullptr;
    std::vector< double >*                track_end_x          = nullptr;
    std::vector< double >*                track_end_x_t0corr   = nullptr;
    std::vector< double >*                track_end_y          = nullptr;
    std::vector< double >*                track_end_z          = nullptr;
    std::vector< double >*                track_t0             = nullptr;
    std::vector< double >*                track_t0_tick_shift  = nullptr;
    std::vector< double >*                track_t0_x_shift     = nullptr;
    std::vector< std::vector< double > >* hit_peak_time        = nullptr;
    std::vector< std::vector< double > >* hit_peak_time_t0corr = nullptr;
    std::vector< std::vector< double > >* hit_peak_time_stddev = nullptr;
    std::vector< std::vector< double > >* hit_rms              = nullptr;
    std::vector< std::vector< double > >* hit_charge           = nullptr;
    std::vector< std::vector< int    > >* hit_multiplicity     = nullptr;
    std::vector< std::vector< int    > >* hit_view             = nullptr;
    std::vector< std::vector< int    > >* hit_channel          = nullptr;
    std::vector< std::vector< int    > >* hit_maximum_tick     = nullptr;
    std::vector< std::vector< double > >* sp_x                 = nullptr;
    std::vector< std::vector< double > >* sp_y                 = nullptr;
    std::vector< std::vector< double > >* sp_z                 = nullptr;
    std::vector< std::vector< double > >* sp_x_t0              = nullptr;
    std::vector< std::vector< double > >* wvfm_pulse_height    = nullptr;
    std::vector< std::vector< double > >* wvfm_fit_mean        = nullptr;
    std::vector< std::vector< double > >* wvfm_fit_sigma       = nullptr;
    std::vector< std::vector< double > >* wvfm_fit_chisq       = nullptr;
    std::vector< std::vector< double > >* wvfm_tick_correction = nullptr;
    std::vector< std::vector< int    > >* wvfm_bin_no          = nullptr;

    TVector3 track_start;
    TVector3 track_start_t0corr;
    TVector3 track_end;
    TVector3 track_end_t0corr;
    float trunc_mean; // Truncated mean has to be a float

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
    bool        use_t0tagged_tracks;
    bool        make_sigma_map;
    bool        fill_tree;
    double      sigma_cut;
    double      wvfm_pulse_height_cut;
    float       drift_velocity;
    double      hit_GOF_cut;
    int         hit_multiplicity_cut;
    int         waveform_size;
    int         waveform_intime_start;
    int         waveform_intime_end;
    int         number_time_bins;
    float       drift_distance;
    float       peak_finder_threshold;
    float       track_rms_cut;
    float       track_trans_dist_cut;
		float       track_theta_xz_cut_low;
		float       track_theta_xz_cut_high;
    int         number_dropped_ticks;

    std::vector<art::TFileDirectory> theseTDs;

    // manipulation histograms
    // binnings are placeholders and will be modified later

    // histogram opened around the hit peak value 
    TH1D* h_wire_in_window = tfs->make<TH1D>("h_wire_in_window", "", 100, 0, 100);

    // after baseline correcting
    TH1D* h_wire_baseline_corrected = tfs->make<TH1D>("h_wire_baseline_corrected", "", 100, 0, 100);

    // For dynamic sigma cut
    std::vector<double> sigmaMedians;
    std::vector<double> pulseHeightMedians;
    //std::vector<double> sigmaMaxs;
    //std::vector<double> pulseHeightMaxs;

    // For truncated mean calculation; needs to be float
    std::vector<float> sigmaDistsPerBin;

    // Hit vectors - separated by plane
    std::vector<recob::Hit> hitsPlane0;
    std::vector<recob::Hit> hitsPlane1;
    std::vector<recob::Hit> hitsPlane2;

    // Output histograms - separated per plane
    std::vector<TH1D*>              h_nWvfmsInBin;
    std::vector<std::vector<TH1D*>> h_summed_wire_info_per_bin; 
    std::vector<std::vector<TH1D*>> h_sigma_hists; 
    std::vector<std::vector<TH1D*>> h_wvfm_pulse_height_hists;
    std::vector<TH1D*>              h_sigma_hist_medians;
    std::vector<TH1D*>              h_sigma_hist_maxs;
    std::vector<TH2D*>              h_sigma_v_bin_precut;
    std::vector<TH2D*>              h_sigma_v_bin_postcut;
    std::vector<TH2D*>              h_wvfm_pulse_height_v_bin_precut;
    std::vector<TH2D*>              h_wvfm_pulse_height_v_bin_postcut;
    std::vector<TH2D*>              h_sigma_v_wvfm_pulse_height_precut;
    std::vector<TH2D*>              h_sigma_v_wvfm_pulse_height_postcut;
    std::vector<TH2D*>              h_track_theta_xz_v_bin;
    std::vector<TH2D*>              h_track_theta_yz_v_bin;

    // other classes
    diffmod::WaveformFunctions _waveform_func;
    diffmod::Utilities         _util;
    TruncMean                  _trunc_mean_func;
    ubana::FiducialVolume      _fiducial_vol;

};


diffmod::LArDiffusion::LArDiffusion(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{

  // defaults set to be MCC9 defaults
  track_label              = p.get< std::string  > ("TrackLabel"           , "pandora");
  wire_label               = p.get< std::string  > ("WireLabel"            , "butcher");
  hit_label                = p.get< std::string  > ("HitLabel"             , "gaushit");
  track_hit_assn           = p.get< std::string  > ("TrackHitAssn"         , "pandora");
  hit_sp_assn              = p.get< std::string  > ("HitSpacePointAssn"    , "pandora");
  track_t0_assn            = p.get< std::string  > ("TrackT0Assn"          , "t0reco");
  sigma_map_file_path      = p.get< std::string  > ("SigmaMapFilePath"     , "");
  sigma_map_directory_file = p.get< std::string  > ("SigmaMapFileDirectory", "DiffusionModule");
  use_t0tagged_tracks      = p.get< bool         > ("UseT0TaggedTracks"    , true);
  make_sigma_map           = p.get< bool         > ("MakeSigmaMap"         , false);
  fill_tree                = p.get< bool         > ("FillTree"             , false);
  sigma_cut                = p.get< double       > ("SigmaCut"             , 1.0);
  wvfm_pulse_height_cut    = p.get< double       > ("PulseHeightCut"       , 100.0);
  hit_GOF_cut              = p.get< double       > ("HitGOFCut"            , 1.1);
  peak_finder_threshold    = p.get< float        > ("PeakFinderThreshold"  , 3.0);
  track_rms_cut            = p.get< float        > ("TrackRMSCut"          , 1000);
  track_trans_dist_cut     = p.get< float        > ("TrackTransDistCut"    , 1000);
	track_theta_xz_cut_low   = p.get< float        > ("TrackThetaXZCutLow"   , 0);       
	track_theta_xz_cut_high  = p.get< float        > ("TrackThetaXZCutHigh"  , 0);       
  hit_multiplicity_cut     = p.get< int          > ("HitMultiplicityCut"   , 1);
  number_time_bins         = p.get< int          > ("NumberTimeBins"       , 25);

  fhicl::ParameterSet const p_const = p.get<fhicl::ParameterSet>("Constants");
  drift_velocity           = p_const.get< float        > ("DriftVelocity");
  waveform_size            = p_const.get< int          > ("WaveformSize"         , 6400);
  waveform_intime_start    = p_const.get< int          > ("WaveformIntimeStart"  , 800);
  waveform_intime_end      = p_const.get< int          > ("WaveformIntimeEnd"    , 5400);
  number_dropped_ticks     = p_const.get< int          > ("NumberDroppedTicks"   , 2400);
  waveform_drift_size      = waveform_intime_end - waveform_intime_start; // 4600
  number_ticks_per_bin     = waveform_drift_size/number_time_bins;        // 184

  MF_LOG_VERBATIM("LArDiffusion") 
    << "PRINTING FHICL CONFIGURATION"
    << "\n-- track_label            : " << track_label
    << "\n-- wire_label             : " << wire_label
    << "\n-- hit_label              : " << hit_label
    << "\n-- track_hit_assn         : " << track_hit_assn
    << "\n-- hit_sp_assn            : " << hit_sp_assn
    << "\n-- track_t0_assn          : " << track_t0_assn
    << "\n-- drift_velocity         : " << drift_velocity
    << "\n-- use_t0tagged_tracks    : " << use_t0tagged_tracks
    << "\n-- make_sigma_map         : " << make_sigma_map
    << "\n-- fill_tree              : " << fill_tree
    << "\n-- sigma_cut              : " << sigma_cut
    << "\n-- wvfm_pulse_height_cut  : " << wvfm_pulse_height_cut
    << "\n-- hit_GOF_cut            : " << hit_GOF_cut
    << "\n-- peak_finder_threshold  : " << peak_finder_threshold
    << "\n-- track rms cut          : " << track_rms_cut
    << "\n-- track trans dist cut   : " << track_trans_dist_cut
		<< "\n-- track_theta_xz_cut_low : " << track_theta_xz_cut_low
		<< "\n-- track_theta_xz_cut_high: " << track_theta_xz_cut_high
    << "\n-- hit_multiplicity_cut   : " << hit_multiplicity_cut
    << "\n-- waveform_size          : " << waveform_size
    << "\n-- waveform_intime_start  : " << waveform_intime_start
    << "\n-- waveform_intime_end    : " << waveform_intime_end
    << "\n-- number_time_bins       : " << number_time_bins
    << "\n-- number_dropped_ticks   : " << number_dropped_ticks
    << "\n-- waveform_drift_size    : " << waveform_drift_size
    << "\n-- number_ticks_per_bin   : " << number_ticks_per_bin;

  // define fiducial volume for analysis
  fhicl::ParameterSet const p_fv = p.get<fhicl::ParameterSet>("FiducialVolume");

  _fiducial_vol.Configure(p_fv,
      geo->DetHalfHeight(),
      2.*geo->DetHalfWidth(),
      geo->DetLength());

  _fiducial_vol.PrintConfig();

}

void diffmod::LArDiffusion::analyze(art::Event const & e) {
  run     = e.run();
  sub_run = e.subRun();
  event   = e.event();

  this->clearVectors();

  MF_LOG_VERBATIM("LArDiffusion::analyze")
    << " --- Processing event " 
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

    // if using a cosmic sample, ensure only a single t0 is
    // associated with the given track, else there's confusion in what
    // the t0 is
    std::vector< art::Ptr< anab::T0 > > t0_from_track;
    if (use_t0tagged_tracks) {
      t0_from_track = t0_from_tracks.at(thisTrack.key());

      if (t0_from_track.size() != 1) {
        continue;
      }
      else if (t0_from_track.size() == 1) {
        art::Ptr< anab::T0 > thisT0 = t0_from_track.at(0);
        track_t0->push_back(thisT0->Time()); // time in us
      }

    }
    else {
      // TODO: Make sure this matches the t0 in the generator fcl file if using single muons 
      //track_t0->push_back(waveform_intime_start); // time in us
      track_t0->push_back(0.); // time in us
    }

    track_t0_x_shift    ->push_back(track_t0->back() * drift_velocity); // convert to x offset
    track_t0_tick_shift ->push_back(track_t0->back() * 2);              // convert to ticks

    track_theta_xz->push_back(_util.getThetaXZ(thisTrack));
    track_theta_yz->push_back(_util.getThetaYZ(thisTrack));
    track_length        ->push_back(thisTrack->Length());
    track_cos_theta     ->push_back(thisTrack->Theta());
    track_start         = thisTrack->Start<TVector3>();
    track_end           = thisTrack->End<TVector3>();
    track_start_x       ->push_back(track_start.X());
    track_start_x_t0corr->push_back(track_start.X() - track_t0_x_shift->back());
    track_start_y       ->push_back(track_start.Y());
    track_start_z       ->push_back(track_start.Z());
    track_end_x         ->push_back(track_end.X());
    track_end_x_t0corr  ->push_back(track_end.X() - track_t0_x_shift->back());
    track_end_y         ->push_back(track_end.Y());
    track_end_z         ->push_back(track_end.Z());

    // get information about how stright the track is
    // do this by taking the dot product of each trajecory point direction
    // and the first trajectory point direction, and then take the RMS
    std::vector<double> dotProds;
    dotProds.resize(0);
    std::vector<double> transDists;
    transDists.resize(0);

    // Get direction vector between first and last point
    TVector3 track_dirv = track_end - track_start;
    for (size_t i_pt = 0; i_pt < thisTrack->CountValidPoints(); i_pt++){

      if (!thisTrack->HasValidPoint(i_pt)) continue;

      recob::Track::Vector_t firstDir = thisTrack->DirectionAtPoint(0);
      recob::Track::Vector_t thisDir  = thisTrack->DirectionAtPoint(i_pt);

      double dotProd = firstDir.Dot(thisDir);
      dotProds.push_back(dotProd);
      
      // Find track point at i_th position, and its components
      recob::Track::TrajectoryPoint_t thisTrackTrajPoint = thisTrack->TrajectoryPoint(i_pt);
      TVector3 trajPoint(thisTrackTrajPoint.position.X(),
                         thisTrackTrajPoint.position.Y(),
                         thisTrackTrajPoint.position.Z()
      );

      // Find minimum distance between trajectory point and a point on 
      // the straight line connecting the start and end points. For derivation,
      // see http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
      double thisTransDist = ((trajPoint - track_start).Cross(trajPoint - track_end)).Mag() / 
                             ((track_end - track_start).Mag() );
      transDists.push_back(thisTransDist);

    }

    track_avg_trans_dist->push_back(TMath::Mean(transDists.size(), &transDists[0])); 

    // ensure track straightness
    if (track_avg_trans_dist->back() > track_trans_dist_cut) 
      continue;
    // Why is this still here???
		//if (track_theta_xz->back() < track_theta_xz_cut_low || track_theta_xz->back() > track_theta_xz_cut_high)
	  //	continue;

    std::vector< art::Ptr< recob::Hit > > hits_from_track = hits_from_tracks.at(thisTrack.key());

    std::vector<double> thit_peak_time        = {};
    std::vector<double> thit_peak_time_t0corr = {};
    std::vector<double> thit_peak_time_stddev = {};
    std::vector<double> thit_rms              = {};
    std::vector<double> thit_charge           = {};
    std::vector<int>    thit_multiplicity     = {};
    std::vector<int>    thit_view             = {};
    std::vector<int>    thit_channel          = {};
    std::vector<int>    thit_maximum_tick     = {};
    std::vector<double> tsp_x                 = {};
    std::vector<double> tsp_y                 = {};
    std::vector<double> tsp_z                 = {};
    std::vector<double> tsp_x_t0              = {};
    std::vector<double> twvfm_pulse_height    = {};
    std::vector<double> twvfm_fit_mean        = {};
    std::vector<double> twvfm_fit_sigma       = {};
    std::vector<double> twvfm_fit_chisq       = {};
    std::vector<double> twvfm_tick_correction = {};
    std::vector<int>    twvfm_bin_no          = {};

    // loop hits
    for (size_t i_hit = 0; i_hit < hits_from_track.size(); i_hit++){

      art::Ptr< recob::Hit > thisHit = hits_from_track.at(i_hit);
      std::vector< art::Ptr< recob::SpacePoint > > sps_from_hit = sp_from_hits.at(thisHit.key());

      // if there's more than one hit associated with the spacepoint, 
      // assume they're localised to a small region and take only the 
      // first one
      if (sps_from_hit.size() == 0){
        MF_LOG_VERBATIM("LArDiffusion")
          << "zero spacepoints associated with hit, skip this hit";
        continue;
      }
      if (sps_from_hit.size() != 1){
        MF_LOG_VERBATIM("LArDiffusion")
          << sps_from_hit.size() 
          << " spacepoints from the hit, just taking first one";
      }
      art::Ptr< recob::SpacePoint > thisSpacePoint = sps_from_hit.at(0);

      const double* spXYZ = thisSpacePoint->XYZ();

      bool  isInFV = _fiducial_vol.InFV(spXYZ[0] - track_t0_x_shift->back(),
                                        spXYZ[1],
                                        spXYZ[2]);

      if (!use_t0tagged_tracks)
        isInFV = _fiducial_vol.InFV(spXYZ[0],
                                    spXYZ[1],
                                    spXYZ[2]);

      // if hit selection is not passed then ignore the hit
      if (!isInFV) continue;
      if (!_waveform_func.passesHitSelection(thisHit, 
            hit_GOF_cut, 
            hit_multiplicity_cut)) continue;

      tsp_x   .push_back(spXYZ[0]);
      tsp_x_t0.push_back(tsp_x.back() - track_t0_x_shift->back());
      tsp_y   .push_back(spXYZ[1]);
      tsp_z   .push_back(spXYZ[2]);

      // looping over planes 
      //for (int thisHit->View() = 0; thisHit->View() < 3; thisHit->View()++){

        // get wire information for hit
        art::Ptr< recob::Wire > wire_from_hit;

        for (size_t i_w = 0; i_w < wire_ptr_vector.size(); i_w++) {

          if ( wire_ptr_vector.at(i_w)->Channel() == thisHit->Channel())
            wire_from_hit = wire_ptr_vector.at(i_w);

        }

        thit_peak_time       .push_back(thisHit->PeakTime());
        thit_peak_time_t0corr.push_back(thisHit->PeakTime() - track_t0_tick_shift->back());

        if (thit_peak_time_t0corr.back() > 3000 && thit_peak_time_t0corr.back() < 4500 && printHitMsg) {
          MF_LOG_VERBATIM("LArDiffusion::anlyze") 
            << "Weird hit_rms value in run " << run 
            << " subrun " << sub_run 
            << " event " << event
            << "\nChannel " << thisHit->Channel();
          printHitMsg = false;
        }

        thit_peak_time_stddev.push_back(thisHit->SigmaPeakTime());
        thit_rms             .push_back(thisHit->RMS());
        thit_charge          .push_back(thisHit->Integral());
        thit_view            .push_back(thisHit->View());
        thit_channel         .push_back(thisHit->Channel());
        thit_multiplicity    .push_back(thisHit->Multiplicity());

        tick_window_size  = number_ticks_per_bin;
        tick_window_left  = thit_peak_time.back() - tick_window_size/2;
        tick_window_right = thit_peak_time.back() + tick_window_size/2;

        // make sure that the window stops at the edge of the waveform
        if (tick_window_left < 0){
          tick_window_left = 0;
          tick_window_size = tick_window_right; 
        }
        if (tick_window_right > waveform_size){
          tick_window_right = waveform_size;
          tick_window_size  = tick_window_right - tick_window_left;
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
          thit_maximum_tick.push_back(h_wire_in_window->GetMaximumBin() + tick_window_left - track_t0_tick_shift->back());
        }
        else {
          thit_maximum_tick.push_back(h_wire_in_window->GetMaximumBin() + tick_window_left);
        }

        MF_LOG_DEBUG("LArDiffusion")
          << "PRINTING INFORMATION FOR HIT " << i_hit
          << "\n-- hit_peak_time        : "  << thit_peak_time.back()
          << "\n-- t0_tick_shift        : "  << track_t0_tick_shift->back()
          << "\n-- hit_peak_time_t0corr : "  << thit_peak_time_t0corr.back()
          << "\n-- tick_window_size     : "  << tick_window_size
          << "\n-- tick_window_left     : "  << tick_window_left
          << "\n-- tick_window_right    : "  << tick_window_right
          << "\n-- hit_maximum_tick     : "  << thit_maximum_tick.back()
          << "\n-- hit_view             : "  << thit_view.back()
          << "\n-- hit_channel          : "  << thit_channel.back();

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
          if (thit_maximum_tick.back() >= binEdgeLeft && thit_maximum_tick.back() < binEdgeRight) {

            twvfm_bin_no.push_back(bin_it);

            MF_LOG_DEBUG("LArDiffusion")
              << "-- Falls into bin "    << twvfm_bin_no.back()
              << "\n-- histogram set to"
              << "\n---- binEdgeLeft: "  << binEdgeLeft
              << "\n---- binEdgeRight: " << binEdgeRight;

            //h_wire_in_window->GetXaxis()->SetLimits(bin_it * number_ticks_per_bin, (bin_it +1) * number_ticks_per_bin);

            // apply baseline correction
            h_wire_baseline_corrected = 
              _waveform_func.applyGlobalBaselineCorrection(
                  h_wire_in_window, 
                  h_wire_baseline_corrected); 

            // calculate sigma
            twvfm_pulse_height.push_back(h_wire_baseline_corrected->GetMaximum());
            twvfm_fit_mean    .push_back(_waveform_func.getSigma(h_wire_baseline_corrected).at(0));
            twvfm_fit_sigma   .push_back(_waveform_func.getSigma(h_wire_baseline_corrected).at(1));
            twvfm_fit_chisq   .push_back(_waveform_func.getSigma(h_wire_baseline_corrected).at(2));


            if (make_sigma_map) {
              h_sigma_hists                     .at(thisHit->View()).at(twvfm_bin_no.back())->Fill(twvfm_fit_sigma.back());
              h_wvfm_pulse_height_hists         .at(thisHit->View()).at(twvfm_bin_no.back())->Fill(twvfm_pulse_height.back());
              h_sigma_v_bin_precut              .at(thisHit->View())->Fill(twvfm_bin_no.back()   , twvfm_fit_sigma.back());
              h_wvfm_pulse_height_v_bin_precut  .at(thisHit->View())->Fill(twvfm_bin_no.back()   , twvfm_pulse_height.back());
              h_sigma_v_wvfm_pulse_height_precut.at(thisHit->View())->Fill(twvfm_fit_sigma.back(), twvfm_pulse_height.back());

              h_track_theta_xz_v_bin.at(thisHit->View())->Fill(twvfm_bin_no.back(), track_theta_xz->back());
              h_track_theta_yz_v_bin.at(thisHit->View())->Fill(twvfm_bin_no.back(), track_theta_yz->back());
              
              //if (track_theta_xz->back() < 90)
              //  h_track_theta_xz_v_bin.at(thisHit->View())->Fill(twvfm_bin_no.back(), track_theta_xz->back());
              //if (track_theta_xz->back() > 90)
              //  h_track_theta_xz_v_bin.at(thisHit->View())->Fill(twvfm_bin_no.back(), 180-track_theta_xz->back());
              //if (track_theta_yz->back() < 90)
              //  h_track_theta_yz_v_bin.at(thisHit->View())->Fill(twvfm_bin_no.back(), track_theta_yz->back());
              //if (track_theta_yz->back() > 90)
              //  h_track_theta_yz_v_bin.at(thisHit->View())->Fill(twvfm_bin_no.back(), 180-track_theta_yz->back());
            }

            else {
              // For comparison
              h_sigma_v_bin_precut              .at(thisHit->View())->Fill(twvfm_bin_no.back()   , twvfm_fit_sigma.back());
              h_wvfm_pulse_height_v_bin_precut  .at(thisHit->View())->Fill(twvfm_bin_no.back()   , twvfm_pulse_height.back());
              h_sigma_v_wvfm_pulse_height_precut.at(thisHit->View())->Fill(twvfm_fit_sigma.back(), twvfm_pulse_height.back());
              h_sigma_hist_medians              .at(thisHit->View())->Fill(sigmaMedians.at(bin_it));
              //h_sigma_hist_maxs                 .at(thisHit->View())->Fill(sigmaMaxs   .at(bin_it));

              // Dynamic sigma cut: check if pulseHeight, sigma, 
              // fall within some region around the median
              /*
              double sigma_lowerLimit = 
              sigmaMedians.at(bin_it) - sigma_cut * h_sigma_hists.at(twvfm_bin_no)->GetStdDev();
              double sigma_higherLimit = 
              sigmaMedians.at(bin_it) + sigma_cut * h_sigma_hists.at(twvfm_bin_no)->GetStdDev();
              double pulseHeight_lowerLimit = 
              pulseHeightMedians.at(bin_it) - wvfm_pulse_height_cut * h_wvfm_pulse_height_hists.at(wvfm_bin_no)->GetStdDev();
              double pulseHeight_higherLimit = 
              pulseHeightMedians.at(bin_it) + wvfm_pulse_height_cut * h_wvfm_pulse_height_hists.at(wvfm_bin_no)->GetStdDev();
              */

              double sigma_lowerLimit = 
                sigmaMedians.at(bin_it) - sigma_cut * h_sigma_hists.at(thisHit->View()).at(twvfm_bin_no.back())->GetStdDev();
              double sigma_higherLimit = 
                sigmaMedians.at(bin_it) + sigma_cut * h_sigma_hists.at(thisHit->View()).at(twvfm_bin_no.back())->GetStdDev();
              double pulseHeight_lowerLimit = 
                pulseHeightMedians.at(bin_it) - wvfm_pulse_height_cut * h_wvfm_pulse_height_hists.at(thisHit->View()).at(twvfm_bin_no.back())->GetStdDev();
              double pulseHeight_higherLimit = 
                pulseHeightMedians.at(bin_it) + wvfm_pulse_height_cut * h_wvfm_pulse_height_hists.at(thisHit->View()).at(twvfm_bin_no.back())->GetStdDev();

              // Maximum 
              /*
              double sigma_lowerLimit = 
                sigmaMaxs.at(bin_it) - sigma_cut * h_sigma_hists.at(thisHit->View()).at(twvfm_bin_no.back())->GetStdDev();
              double sigma_higherLimit = 
                sigmaMaxs.at(bin_it) + sigma_cut * h_sigma_hists.at(thisHit->View()).at(twvfm_bin_no.back())->GetStdDev();
              double pulseHeight_lowerLimit = 
                pulseHeightMaxs.at(bin_it) - wvfm_pulse_height_cut * h_wvfm_pulse_height_hists.at(thisHit->View()).at(twvfm_bin_no.back())->GetStdDev();
              double pulseHeight_higherLimit = 
                pulseHeightMaxs.at(bin_it) + wvfm_pulse_height_cut * h_wvfm_pulse_height_hists.at(thisHit->View()).at(twvfm_bin_no.back())->GetStdDev();
                */

              if (twvfm_fit_sigma.back() < sigma_lowerLimit 
                  || twvfm_fit_sigma.back() > sigma_higherLimit 
                  || twvfm_pulse_height.back() < pulseHeight_lowerLimit 
                  || twvfm_pulse_height.back() > pulseHeight_higherLimit) {
                continue;
              }

              h_sigma_hists                      .at(thisHit->View()).at(twvfm_bin_no.back())->Fill(twvfm_fit_sigma.back());
              h_wvfm_pulse_height_hists          .at(thisHit->View()).at(twvfm_bin_no.back())->Fill(twvfm_pulse_height.back());
              h_sigma_v_bin_postcut              .at(thisHit->View())->Fill(twvfm_bin_no.back()      , twvfm_fit_sigma.back());
              h_wvfm_pulse_height_v_bin_postcut  .at(thisHit->View())->Fill(twvfm_bin_no.back()      , twvfm_pulse_height.back());
              h_sigma_v_wvfm_pulse_height_postcut.at(thisHit->View())->Fill(twvfm_pulse_height.back(), twvfm_fit_sigma.back());

              h_track_theta_xz_v_bin.at(thisHit->View())->Fill(twvfm_bin_no.back(), track_theta_xz->back());
              h_track_theta_yz_v_bin.at(thisHit->View())->Fill(twvfm_bin_no.back(), track_theta_yz->back());

              //h_nWvfmsInBin.at(thisHit->View())->Fill(twvfm_bin_no.back(), 1);
              h_nWvfmsInBin.at(thisHit->View())->Fill(twvfm_bin_no.back() );

              //if (track_theta_xz->back() < 90)
              //  h_track_theta_xz_v_bin.at(thisHit->View())->Fill(twvfm_bin_no.back(), track_theta_xz->back());
              //if (track_theta_xz->back() > 90)
              //  h_track_theta_xz_v_bin.at(thisHit->View())->Fill(twvfm_bin_no.back(), 180-track_theta_xz->back());
              //if (track_theta_yz->back() < 90)
              //  h_track_theta_yz_v_bin.at(thisHit->View())->Fill(twvfm_bin_no.back(), track_theta_yz->back());
              //if (track_theta_yz->back() > 90)
              //  h_track_theta_yz_v_bin.at(thisHit->View())->Fill(twvfm_bin_no.back(), 180-track_theta_yz->back());

              // Now find the shift (in ticks) needed to minimise the rms 
              // of the sum of the histograms
              if (h_summed_wire_info_per_bin.at(thisHit->View()).at(bin_it)->Integral() == 0){
                twvfm_tick_correction.push_back(0);
              }
              else {
                twvfm_tick_correction.push_back(  
                  _waveform_func.findXCorrection(h_summed_wire_info_per_bin.at(thisHit->View()).at(bin_it), 
                      h_wire_baseline_corrected, 
                      number_ticks_per_bin, 
                      twvfm_fit_mean.back()));
              }

              TH1D* h_wvfm_tick_correction = 
                new TH1D("h_wvfm_tick_correction", 
                    "", 
                    number_ticks_per_bin, 
                    h_wire_baseline_corrected->GetXaxis()->GetXmin(), 
                    h_wire_baseline_corrected->GetXaxis()->GetXmax()); 

              for (int ntick = 1; ntick <= h_wire_baseline_corrected->GetNbinsX(); ntick++)
                h_wvfm_tick_correction->SetBinContent(
                    ntick, 
                    h_wire_baseline_corrected->GetBinContent(ntick+twvfm_tick_correction.back()));

              // finally add to output histograms
              h_summed_wire_info_per_bin.at(thisHit->View()).at(bin_it)->Add(h_wvfm_tick_correction);

            } // !make_sigma_map
          } // if maximum tick in bin
        } // loop bins
      //} // loop planes
    } // loop hits

    sp_x                ->push_back(tsp_x                 );
    sp_x_t0             ->push_back(tsp_x_t0              );
    sp_y                ->push_back(tsp_y                 );
    sp_z                ->push_back(tsp_z                 );
    hit_peak_time       ->push_back(thit_peak_time        );
    hit_peak_time_t0corr->push_back(thit_peak_time_t0corr );
    hit_peak_time_stddev->push_back(thit_peak_time_stddev );
    hit_rms             ->push_back(thit_rms              );
    hit_charge          ->push_back(thit_charge           );
    hit_multiplicity    ->push_back(thit_multiplicity     );
    hit_view            ->push_back(thit_view             );
    hit_channel         ->push_back(thit_channel          );
    sp_x                ->push_back(tsp_x                 );
    sp_y                ->push_back(tsp_y                 );
    sp_z                ->push_back(tsp_z                 );
    sp_x_t0             ->push_back(tsp_x_t0              );
    wvfm_pulse_height   ->push_back(twvfm_pulse_height    );
    wvfm_fit_mean       ->push_back(twvfm_fit_mean        );
    wvfm_fit_sigma      ->push_back(twvfm_fit_sigma       );
    wvfm_fit_chisq      ->push_back(twvfm_fit_chisq       );
    wvfm_tick_correction->push_back(twvfm_tick_correction );
    wvfm_bin_no         ->push_back(twvfm_bin_no          );

  } // loop tracks

  if (fill_tree) difftree->Fill();

} // LArDiffusion::analyze

void diffmod::LArDiffusion::beginJob()
{

  h_summed_wire_info_per_bin.resize(3);
  h_sigma_hists             .resize(3);
  h_wvfm_pulse_height_hists .resize(3);

  difftree = tfs->make<TTree>("difftree", "diffusion tree");

  difftree->Branch("run"                  , &run);
  difftree->Branch("sub_run"              , &sub_run);
  difftree->Branch("event"                , &event);
  difftree->Branch("track_length"         , "std::vector<double>"                  , &track_length);
  difftree->Branch("track_avg_trans_dist" , "std::vector<double>"                  , &track_avg_trans_dist);
  difftree->Branch("track_cos_theta"      , "std::vector<double>"                  , &track_cos_theta);
  difftree->Branch("track_theta_xz"       , "std::vector<double>"                  , &track_theta_xz);
  difftree->Branch("track_theta_yz"       , "std::vector<double>"                  , &track_theta_yz);
  difftree->Branch("track_start_x"        , "std::vector<double>"                  , &track_start_x);
  difftree->Branch("track_start_x_t0corr" , "std::vector<double>"                  , &track_start_x_t0corr);
  difftree->Branch("track_start_y"        , "std::vector<double>"                  , &track_start_y);
  difftree->Branch("track_start_z"        , "std::vector<double>"                  , &track_start_z);
  difftree->Branch("track_end_x"          , "std::vector<double>"                  , &track_end_x);
  difftree->Branch("track_end_x_t0corr"   , "std::vector<double>"                  , &track_end_x_t0corr);
  difftree->Branch("track_end_y"          , "std::vector<double>"                  , &track_end_y);
  difftree->Branch("track_end_z"          , "std::vector<double>"                  , &track_end_z);
  difftree->Branch("track_t0"             , "std::vector<double>"                  , &track_t0);
  difftree->Branch("track_t0_x_shift"     , "std::vector<double>"                  , &track_t0_x_shift);
  difftree->Branch("track_t0_tick_shift"  , "std::vector<double>"                  , &track_t0_tick_shift);
  difftree->Branch("hit_peak_time"        , "std::vector< std::vector< double > >" , &hit_peak_time);
  difftree->Branch("hit_peak_time_t0corr" , "std::vector< std::vector< double > >" , &hit_peak_time_t0corr);
  difftree->Branch("hit_peak_time_stddev" , "std::vector< std::vector< double > >" , &hit_peak_time_stddev);
  difftree->Branch("hit_rms"              , "std::vector< std::vector< double > >" , &hit_rms);
  difftree->Branch("hit_charge"           , "std::vector< std::vector< double > >" , &hit_charge);
  difftree->Branch("hit_multiplicity"     , "std::vector< std::vector< int > >   " , &hit_multiplicity);
  difftree->Branch("hit_view"             , "std::vector< std::vector< int > >   " , &hit_view);
  difftree->Branch("hit_channel"          , "std::vector< std::vector< int > >   " , &hit_channel);
  difftree->Branch("hit_maximum_tick"     , "std::vector< std::vector< int > >   " , &hit_maximum_tick);
  difftree->Branch("sp_x"                 , "std::vector< std::vector< double > >" , &sp_x);
  difftree->Branch("sp_x_t0"              , "std::vector< std::vector< double > >" , &sp_x_t0);
  difftree->Branch("sp_y"                 , "std::vector< std::vector< double > >" , &sp_y);
  difftree->Branch("sp_z"                 , "std::vector< std::vector< double > >" , &sp_z);
  difftree->Branch("wvfm_pulse_height"    , "std::vector< std::vector< double > >" , &wvfm_pulse_height);
  difftree->Branch("wvfm_fit_mean"        , "std::vector< std::vector< double > >" , &wvfm_fit_mean);
  difftree->Branch("wvfm_fit_sigma"       , "std::vector< std::vector< double > >" , &wvfm_fit_sigma);
  difftree->Branch("wvfm_fit_chisq"       , "std::vector< std::vector< double > >" , &wvfm_fit_chisq);
  difftree->Branch("wvfm_tick_correction" , "std::vector< std::vector< double > >" , &wvfm_tick_correction);
  difftree->Branch("wvfm_bin_no"          , "std::vector< std::vector< int > >   " , &wvfm_bin_no);

  std::vector<std::string> folderNames ={
    "plane0",
    "plane1",
    "plane2"
  };

  MF_LOG_VERBATIM("LArDiffusion") << "beginning loop over planes";

  for (size_t ifN = 0; ifN < folderNames.size(); ifN++){

    theseTDs.push_back(tfs->mkdir(folderNames.at(ifN), folderNames.at(ifN)));

    MF_LOG_VERBATIM("LArDiffusion") << "made TFileDirectory " << folderNames.at(ifN);

    h_sigma_v_bin_precut.push_back(theseTDs.back().make<TH2D>(
          ("h_sigma_v_bin_precut"+folderNames.at(ifN)).c_str(), 
          ";Bin no. ; #sigma_{t}^{2} (#mus^{2});", 
          number_time_bins, 0, number_time_bins, 
          100, 0, 10));

    h_sigma_v_bin_postcut.push_back(theseTDs.back().make<TH2D>(
          ("h_sigma_v_bin_postcut"+folderNames.at(ifN)).c_str(), 
          ";Bin no. ; #sigma_{t}^{2} (#mus^{2});", 
          number_time_bins, 0, number_time_bins, 
          100, 0, 10));

    h_wvfm_pulse_height_v_bin_precut.push_back(theseTDs.back().make<TH2D>(
          ("h_wvfm_pulse_height_v_bin_precut"+folderNames.at(ifN)).c_str(), 
          ";Bin no. ; Pulse Height (Arb. Units);", 
          number_time_bins, 0, number_time_bins, 
          100, 0, 20));

    h_wvfm_pulse_height_v_bin_postcut.push_back(theseTDs.back().make<TH2D>(
          ("h_wvfm_pulse_height_v_bin_postcut"+folderNames.at(ifN)).c_str(), 
          ";Bin no. ; Pulse Height (Arb. Units);", 
          number_time_bins, 0, number_time_bins, 
          100, 0, 20));

    h_sigma_v_wvfm_pulse_height_precut.push_back(theseTDs.back().make<TH2D>(
          ("h_sigma_v_wvfm_pulse_height_precut"+folderNames.at(ifN)).c_str(), 
          ";#sigma_{t}^{2} (#mus^{2}); Pulse Height (Arb. Units);", 
          100, 0, 10, 
          100, 0, 20));

    h_sigma_v_wvfm_pulse_height_postcut.push_back(theseTDs.back().make<TH2D>(
          ("h_sigma_v_wvfm_pulse_height_postcut"+folderNames.at(ifN)).c_str(), 
          ";#sigma_{t}^{2} (#mus^{2}); Pulse Height (Arb. Units);", 
          100, 0, 10, 
          100, 0, 20));

    h_track_theta_xz_v_bin.push_back(theseTDs.back().make<TH2D>(
          ("h_track_theta_xz_v_bin"+folderNames.at(ifN)).c_str(), 
          ";Bin no. ; #theta_{xz} (Deg.);", 
          number_time_bins, 0, number_time_bins, 
          100, 0, 20));

    h_track_theta_yz_v_bin.push_back(theseTDs.back().make<TH2D>(
          ("h_track_theta_yz_v_bin"+folderNames.at(ifN)).c_str(), 
          ";Bin no. ; #theta_{yz} (Deg.);", 
          number_time_bins, 0, number_time_bins, 
          250, 0, 50));

    h_sigma_hist_medians.push_back(theseTDs.back().make<TH1D>(
          ("h_sigma_hist_medians"+folderNames.at(ifN)).c_str(), 
          ";Median #sigma per bin;", 
          number_time_bins, 0, number_time_bins));

    h_sigma_hist_maxs.push_back(theseTDs.back().make<TH1D>(
          ("h_sigma_hist_maxs"+folderNames.at(ifN)).c_str(), 
          ";Max #sigma per bin;", 
          number_time_bins, 0, number_time_bins));

    // analysis junk
    if (!make_sigma_map) {
      //h_single_waveform = tfs->make<TH1D>("h_single_waveform", ";Time (ticks); Arb. Units;", 100, 0, 100);

      //h_nWvfmsInBin = theseTDs.back().make<TH1D>(("h_nWvfmsInBin"+folderNames.at(ifN)).c_str(), ";Drift bin; No. Waveforms;", 25, 0, 25);

      //theseTDs.push_back(tfs->mkdir(folderNames.at(ifN), folderNames.at(ifN)));

      //MF_LOG_VERBATIM("LArDiffusion") << "made TFileDirectory " << folderNames.at(ifN);

      h_nWvfmsInBin.push_back(theseTDs.back().make<TH1D>(
            ("h_nWvfmsInBin"+folderNames.at(ifN)).c_str(), 
            ";Bin no. ; No. Waveforms;", 
            number_time_bins, 0, number_time_bins)
      );

      // Import sigmaMap, assuming it already exists
      std::string filePath;

      cet::search_path sp("FW_SEARCH_PATH");
      if( !sp.find_file(sigma_map_file_path, filePath) )
        throw cet::exception("LArDiffusion")
          << "Cannot find FHC numu quantile file "
          << sigma_map_file_path;

      MF_LOG_VERBATIM("LArDiffusion::beginJob") 
        << "Running without producing sigma map. Checking that it exists..."
        << "\nGetting sigma map from file path " << filePath;
      TFile sigmaMap(filePath.c_str(), "READ");

      if (sigmaMap.IsOpen() == false){
        MF_LOG_VERBATIM("LArDiffusion::beginJob")
          << "No sigma map! Run module using run_diffusion_module_sigmamap.fcl first, " 
          << "\nor check that you're in the right directory.";
      }

      for (int i = 0; i < number_time_bins; i++){

        // declare summed waveforms
        std::string histo_name = 
          "summed_waveform_bin_"  +
          std::to_string(i)       +
          "_"                     +
          folderNames.at(ifN);

        h_summed_wire_info_per_bin.at(ifN).push_back(
            theseTDs.back().make<TH1D>(histo_name.c_str(), 
              "", 
              number_ticks_per_bin, 
              waveform_intime_start + (i * number_ticks_per_bin), 
              waveform_intime_start + ((i+1) * number_ticks_per_bin)));

        sigmaDistsPerBin.resize(0);

        std::string sigmaMapHistoName = "DiffusionModule/"+folderNames.at(ifN)+"/h_sigma_"+std::to_string(i)+"_"+folderNames.at(ifN);
        h_sigma_hists.at(ifN).push_back((TH1D*)sigmaMap.Get(sigmaMapHistoName.c_str()));

        std::string pulseHeightHistoName = "DiffusionModule/"+folderNames.at(ifN)+"/h_wvfm_pulse_height_"+std::to_string(i)+"_"+folderNames.at(ifN);
        h_wvfm_pulse_height_hists.at(ifN).push_back((TH1D*)sigmaMap.Get(pulseHeightHistoName.c_str()));

        // Three options for picking out waveforms, either
        // 1) get waveforms around the median
        // 2) get the peak bin value
        // 3) calculate the truncated mean

        // OPTION 1) Calculate medians in each bin
        sigmaMedians      .push_back(_waveform_func.getMedian(h_sigma_hists.at(ifN).at(i)));
        pulseHeightMedians.push_back(_waveform_func.getMedian(h_wvfm_pulse_height_hists.at(ifN).at(i)));

        // OPTION 2) Calculate maximum in each bin
        //int sigmaMaxBin       = h_sigma_hists.at(ifN).at(i)->GetMaximumBin();
        //int pulseHeightMaxBin = h_wvfm_pulse_height_hists.at(ifN).at(i)->GetMaximumBin();

        //sigmaMaxs      .push_back(h_sigma_hists.at(ifN).at(i)->GetXaxis()->GetBinCenter(sigmaMaxBin));
        //pulseHeightMaxs.push_back(h_wvfm_pulse_height_hists.at(ifN).at(i)->GetXaxis()->GetBinCenter(pulseHeightMaxBin));

        // OPTION 3) Take sigma hist and calculate truncated mean 
        /*
        for (int j = 1; j < h_sigma_hists.at(ifN).at(i)->GetNbinsX()+1; j++) {
          if (h_sigma_hists.at(ifN).at(i)->GetBinContent(j) > 0) {
            for (int k = 0; k < h_sigma_hists.at(ifN).at(i)->GetBinContent(j); k++ ) {
              sigmaDistsPerBin.push_back(h_sigma_hists.at(ifN).at(i)->GetXaxis()->GetBinCenter(j) );
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
        */
        //MF_LOG_VERBATIM("LArDiffusion") << sigmaMedians.at(i) << "\t" << sigmaMaxs.at(i) << "\t" << trunc_mean;
      }
    }
    else {
      for (int n = 0; n < number_time_bins; n++) {
       
        // create histograms

        std::string sigmaHistName = "h_sigma_" + std::to_string(n)+"_"+folderNames.at(ifN);
        h_sigma_hists.at(ifN).push_back(
            theseTDs.back().make<TH1D>(sigmaHistName.c_str(), 
              ";#sigma_{t};", 
              250, 0, 10) );
        
        std::string pulseHeightHistName = "h_wvfm_pulse_height_" + std::to_string(n)+"_"+folderNames.at(ifN);
        h_wvfm_pulse_height_hists.at(ifN).push_back(
            theseTDs.back().make<TH1D>(pulseHeightHistName.c_str(), 
              ";Pulse Height;", 
              250, 0, 20) );
      }
    }
  }
}

void diffmod::LArDiffusion::printHistogram(TH1D* h){
  for (int i = 0; i < h->GetNbinsX(); i++){
    MF_LOG_VERBATIM("LArDiffusion") << i << " " << h->GetBinLowEdge(i) << " " << h->GetBinContent(i);
  }
}

void diffmod::LArDiffusion::clearVectors(){
  track_length          ->resize(0);
  track_avg_trans_dist  ->resize(0);
  track_cos_theta       ->resize(0);
  track_theta_xz        ->resize(0);
  track_theta_yz        ->resize(0);
  track_start_x         ->resize(0);
  track_start_x_t0corr  ->resize(0);
  track_start_y         ->resize(0);
  track_start_z         ->resize(0);
  track_end_x           ->resize(0);
  track_end_x_t0corr    ->resize(0);
  track_end_y           ->resize(0);
  track_end_z           ->resize(0);
  track_t0              ->resize(0);
  track_t0_tick_shift   ->resize(0);
  track_t0_x_shift      ->resize(0);
  hit_peak_time         ->resize(0);
  hit_peak_time_t0corr  ->resize(0);
  hit_peak_time_stddev  ->resize(0);
  hit_rms               ->resize(0);
  hit_charge            ->resize(0);
  hit_multiplicity      ->resize(0);
  hit_view              ->resize(0);
  hit_channel           ->resize(0);
  hit_maximum_tick      ->resize(0);
  sp_x                  ->resize(0);
  sp_y                  ->resize(0);
  sp_z                  ->resize(0);
  sp_x_t0               ->resize(0);
  wvfm_pulse_height     ->resize(0);
  wvfm_fit_mean         ->resize(0);
  wvfm_fit_sigma        ->resize(0);
  wvfm_fit_chisq        ->resize(0);
  wvfm_tick_correction  ->resize(0);
  wvfm_bin_no           ->resize(0);

}

DEFINE_ART_MODULE(diffmod::LArDiffusion)
