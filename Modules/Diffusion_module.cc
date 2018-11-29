////////////////////////////////////////////////////////////////////////
// Class:       Diffusion
// Plugin Type: analyzer (art v2_11_03)
// File:        Diffusion_module.cc
//
// Generated at Thu Nov 29 09:42:17 2018 by Adam Lister using cetskelgen
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

class Diffusion;


class Diffusion : public art::EDAnalyzer {
public:
  explicit Diffusion(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Diffusion(Diffusion const &) = delete;
  Diffusion(Diffusion &&) = delete;
  Diffusion & operator = (Diffusion const &) = delete;
  Diffusion & operator = (Diffusion &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Declare member data here.

};


Diffusion::Diffusion(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void Diffusion::analyze(art::Event const & e)
{
  // Implementation of required member function here.
}

void Diffusion::beginJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(Diffusion)
