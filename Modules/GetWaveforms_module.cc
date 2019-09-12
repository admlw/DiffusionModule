////////////////////////////////////////////////////////////////////////
// Class:       GetWaveforms
// Plugin Type: analyzer (art v3_01_02)
// File:        GetWaveforms_module.cc
//
// Generated at Wed Sep 11 20:21:30 2019 by Adam Lister using cetskelgen
// from cetlib version v3_05_01.
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
#include "lardataobj/RecoBase/Wire.h"

// art
#include "art/Framework/Services/Optional/TFileService.h"

// root
#include "TH1.h"

class GetWaveforms;


class GetWaveforms : public art::EDAnalyzer {
  public:
    explicit GetWaveforms(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    GetWaveforms(GetWaveforms const&) = delete;
    GetWaveforms(GetWaveforms&&) = delete;
    GetWaveforms& operator=(GetWaveforms const&) = delete;
    GetWaveforms& operator=(GetWaveforms&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

    // Selected optional functions.
    void beginJob() override;

    // Get Cantor pairing number
    // -- Cantor pairing number takes in two
    //    integers and calculates a bijectvie NxN->N mapping
    //    Using this to get a unique identified for the 
    //    map of TH1s
    int getCantorPairingNumber(int a, int b);

  private:

    art::ServiceHandle< art::TFileService > tfs;

    // output
    std::map< int, TH1D* > waveformHistograms;

    //fhicl
    std::vector<int> fEventVector;
    std::vector<int> fChannelVector;
    std::string fWireLabel;

    //other
    int run;
    int subRun;
    int event;

};


GetWaveforms::GetWaveforms(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{

  fWireLabel = p.get< std::string >("WireLabel", "diffsel");
  fEventVector = p.get< std::vector<int> >("EventVector");
  fChannelVector = p.get< std::vector<int> >("ChannelVector");

}

void GetWaveforms::analyze(art::Event const& e)
{

  run     = e.run();
  subRun = e.subRun();
  event   = e.event();

  // loop over supplied event list and find event of interest
  bool isFoundEvent = false;
  for (size_t i = 0; i < fEventVector.size(); i++){
    if (fEventVector.at(i) == event)
      isFoundEvent = true;
  }

  if (isFoundEvent == false) return;
  else {
    MF_LOG_VERBATIM("GetWaveforms::analyze")
      << "-- Found event " 
      << event
      << ", processing.";
  }

  // now get handle to wires
  art::Handle< std::vector<recob::Wire> > wireHandle;
  e.getByLabel(fWireLabel, wireHandle);
  std::vector< art::Ptr< recob::Wire > > wirePtrVector;
  art::fill_ptr_vector(wirePtrVector, wireHandle);

  // loop the wires, then find out if the channel is in the 
  // vector of chosen channels
  for (size_t i = 0; i < wirePtrVector.size(); i++){

    art::Ptr< recob::Wire > thisWire = wirePtrVector.at(i);

    int wireChannel = thisWire->Channel();

    for (size_t j = 0; j < fChannelVector.size(); j++){

      if (fChannelVector.at(j) == wireChannel){

        int cantorNumber = this->getCantorPairingNumber(event, 
                                                        wireChannel);

        MF_LOG_VERBATIM("GetWaveforms::analyze")
          << "-- Found match at channel " << wireChannel
          << ", finding histogram using cantor number: " << cantorNumber;

        // retrieve correct histogram using cantor number
        TH1D* hTmp = (TH1D*)waveformHistograms.find(cantorNumber)->second;

        MF_LOG_VERBATIM("GetWaveforms::analyze")
          << "Found histogram: " << hTmp->GetName();

        for (size_t itck = 1; itck < thisWire->NSignal(); itck++){
          hTmp->SetBinContent(itck, thisWire->Signal().at(itck));
        }
      }
    }
  }
}

void GetWaveforms::beginJob()
{

  // fill vector of TH1s
  for (size_t i = 0; i < fEventVector.size(); i++){

    for (size_t j = 0; j < fChannelVector.size(); j++){

      std::string histoName = 
        "waveform_event_" + 
        std::to_string(fEventVector.at(i)) + 
        "_channel_" + 
        std::to_string(fChannelVector.at(j));

      // get cantor number (unique identifier) to act as the key in the map

      int cantorNumber = this->getCantorPairingNumber(fEventVector.at(i), 
                                                      fChannelVector.at(j));

      MF_LOG_VERBATIM("GetWaveforms::beginJob")
        << "creating histogram with cantorNumber: " << cantorNumber
        << ", with name: " << histoName;

      // assume that there's 9600 ticks in the histograms
      // can't quite remember if that's right but close enough

      waveformHistograms.insert( std::pair< int, TH1D* >(
          cantorNumber, 
          tfs->make<TH1D>(histoName.c_str(), 
                          ";Tick No.;ADC*", 
                          9600, 0, 9600)));      
    }

  }

}

int GetWaveforms::getCantorPairingNumber(int a, int b){

  int cantorNumber = (((a+b)*(a+b+1))/2)+b;
  return cantorNumber;

}

DEFINE_ART_MODULE(GetWaveforms)
