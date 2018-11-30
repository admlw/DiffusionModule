////////////////////////////////////////////////////////////////////////
// Class:       DiffusionAssnTester
// Plugin Type: analyzer (art v2_05_00)
// File:        DiffusionAssnTester_module.cc
//
// Generated at Tue Oct 17 12:12:32 2017 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
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

#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "TFile.h"
#include "TH1.h"

class DiffusionAssnTester;


class DiffusionAssnTester : public art::EDAnalyzer {
  public:
    explicit DiffusionAssnTester(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    DiffusionAssnTester(DiffusionAssnTester const &) = delete;
    DiffusionAssnTester(DiffusionAssnTester &&) = delete;
    DiffusionAssnTester & operator = (DiffusionAssnTester const &) = delete;
    DiffusionAssnTester & operator = (DiffusionAssnTester &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;

  private:

    // Declare member data here.
    int trueTrackID;
    art::ServiceHandle< art::TFileService > tfs;

    TH1D* correctedMinusTruthStartX;
    TH1D* correctedMinusTruthEndX;
    
    double sxTraj = -9999;
    double syTraj = -9999;
    double szTraj = -9999;
    double exTraj = -9999;
    double eyTraj = -9999;
    double ezTraj = -9999;

    double sx = -9999;
    double sy = -9999;
    double sz = -9999;
    double ex = -9999;
    double ey = -9999;
    double ez = -9999;

    double exCorrn = -9999;
    double sxCorrn = -9999;


};


DiffusionAssnTester::DiffusionAssnTester(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{}

void DiffusionAssnTester::analyze(art::Event const & e)
{

  // get track handle
  art::Handle< std::vector<recob::Track> > trackHandle;
  e.getByLabel("diffsel", trackHandle);

  art::Handle< std::vector<simb::MCParticle> > mcpHandle;
  e.getByLabel("largeant", mcpHandle);

  // get assns
  art::FindManyP< recob::Hit > trackHitAssn(trackHandle, e, "diffsel");
  art::FindManyP< anab::T0 > trackT0Assn(trackHandle, e, "diffsel");

  int i = 0;
  for (auto const& track : (*trackHandle)){

    std::vector< art::Ptr<recob::Hit> > hits = trackHitAssn.at(i);
    std::vector< art::Ptr<anab::T0> > t0s = trackT0Assn.at(i);

    std::cout << ">> Looping over track with ID " << track.ID()
      << "\n  >> Number of associated hits: " << hits.size()
      << "\n  >> Number of associated t0s : " << t0s.size()
      << std::endl;

    sx = track.Vertex().X();
    sy = track.Vertex().Y();
    sz = track.Vertex().Z();
    ex = track.End().X();
    ey = track.End().Y();
    ez = track.End().Z();

    std::cout << ">> RECONSTRUCTED INFORMAITON "
      << "\n  >> Track starting positions: " 
      << "\n\t  >> " << sx << ", " << ex
      << "\n\t  >> " << sy << ", " << ey
      << "\n\t  >> " << sz << ", " << ez
      << std::endl;

    double t0Correction = 0;
    // perform t0 correction for tracks
    if (t0s.size() == 1){

      auto t0 = t0s.at(0);
      t0Correction = t0->Time() * 0.001 * 0.1114;
      trueTrackID = t0->ID();
    }

    sxCorrn = sx - t0Correction;
    exCorrn = ex - t0Correction;

    std::cout << ">> -T0 CORRECTED INFORMATION "
      << "\n  >> Track starting positions: " 
      << "\n\t  >> " << sxCorrn << ", " << exCorrn
      << "\n\t  >> " << sy << ", " << ey
      << "\n\t  >> " << sz << ", " << ez
      << std::endl;


    i++;
  }

  // take a look at mcp information
  for (auto const& mcp : (*mcpHandle)){
    if (mcp.TrackId() != trueTrackID) continue;

    // MC cosmic truth information starts outside of TPC
    // need to track to inside -> get MCTrajectory obj.
    const simb::MCTrajectory&  traj = mcp.Trajectory();
    bool isInTPC = false;
    for (size_t k = 0; k < mcp.NumberTrajectoryPoints(); k++){

      if (((traj.X(k) > 0 && traj.X(k) < 256) &&
            (traj.Y(k) > -111.5 && traj.Y(k) < 111.5) &&
            (traj.Z(k) > 0 && traj.Z(k) < 1036)) && isInTPC == false){

        isInTPC = true;
        sxTraj = traj.X(k);
        syTraj = traj.Y(k);
        szTraj = traj.Z(k);

      }

      if (((traj.X(k) <= 0 || traj.X(k) >= 256) ||
            (traj.Y(k) <= -111.5 || traj.Y(k) >= 111.5) ||
            (traj.Z(k) <= 0 || traj.Z(k) >= 1036)) && isInTPC == true){

        exTraj = traj.X(k);
        eyTraj = traj.Y(k);
        ezTraj = traj.Z(k);
        break;
      }


    }

    std::cout << ">> TRUTH INFORMATION "
      << "\n  >> Track starting positions: " 
      << "\n\t  >> " << sxTraj << ", " << exTraj
      << "\n\t  >> " << syTraj << ", " << eyTraj
      << "\n\t  >> " << szTraj << ", " << ezTraj
      << std::endl;


  }

  correctedMinusTruthStartX->Fill(sxCorrn - sxTraj);
  correctedMinusTruthEndX->Fill(exCorrn - exTraj);

}

void DiffusionAssnTester::beginJob()
{
  // Implementation of optional member function here.

  correctedMinusTruthStartX = tfs->make<TH1D>("correctedMinusTruthStartX", ";t0-Corrected s_{X} - true s_{X};", 40, -10, 10);
  correctedMinusTruthEndX = tfs->make<TH1D>("correctedMinusTruthEndX", ";t0-Corrected e_{X} - true e_{X};", 40, -10, 10);
}

DEFINE_ART_MODULE(DiffusionAssnTester)
