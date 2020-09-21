#include "GetWaveform.h"

TH1D* average(TH1* h, TH1D* hAverage, int binSize, int nTicks){
    int counter = 0;
    double average = 0;
    for (int i = 0; i < nTicks; i= i+binSize){
        average = 0;
        for (int j = 0; j < binSize; j++){
            average = average + h->GetBinContent(i+j);
        }

        average = average/binSize;
        hAverage->SetBinContent(counter, average);
        counter++;
    }
    return hAverage;
}


int main(int argv, char** argc){

    TString wcRawDProducer(argc[3]);
    art::InputTag wcRawDTag {wcRawDProducer};

    int binSize = 10;

    ///////////////////////////////////////////////////////////////////////////

    //
    // input files
    //

    // construct vector of input files from input textfile

    std::vector<std::string> filenames;
    filenames.push_back(argc[1]);

    ///////////////////////////////////////////////////////////////////////////

    //
    // create output file and histograms
    //


    TFile f_output(argc[2], "UPDATE");

    int selectedRun = atoi(argc[4]);
    
    int minimumEvent = atoi(argc[5]);
    int maximumEvent = atoi(argc[6]);

    // construct a list of channel numbers to run over
    int minimumChannel = atoi(argc[7]);
    int maximumChannel = atoi(argc[8]);

    int nTicks = atoi(argc[9]);

    std::vector<int> channelList;
    for (int i = minimumChannel; i <= maximumChannel; i++){

        channelList.push_back(i);

    }
    

    TH1D* waveform = new TH1D("waveform", "", nTicks, 0, nTicks);
    TH1* frqspace = 0;
    TH1D* hAverage = new TH1D("hAverage", "", nTicks/binSize, 0, nTicks);
    TH1D* ADCs = new TH1D("ADCs", "", 4000, -2000, 2000);

    ///////////////////////////////////////////////////////////////////////////

    //
    // begin event loop
    //

    for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()){

        int run = ev.eventAuxiliary().run();
        int event = ev.eventAuxiliary().event();

        if (/*selectedRun != run ||*/ (event < minimumEvent || event > maximumEvent)) continue;

        TString dirFileName = Form("Event_%i", event);
        f_output.mkdir(dirFileName);

        const auto& rawDHandle = ev.getValidHandle< std::vector<raw::RawDigit> >(wcRawDTag);

        for (auto const& rawD : (*rawDHandle)){

            int channel = rawD.Channel();

            if (channel < minimumChannel || channel > maximumChannel) continue;

            TString chanNo = Form("Event_%i_TimeWfm_channel%i", event, channel);
            waveform->SetName(chanNo);
            for (size_t i = 0; i < nTicks; i++){
                waveform->SetBinContent(i, rawD.ADC(i)/*-rawD.GetPedestal()*/);
                ADCs->Fill(rawD.ADC(i)/*-rawD.GetPedestal()*/);
            }

            //std::cout << "Event " << event << " channel " << channel << std::endl;
            //std::cout << "Waveform peak " << waveform->GetMaximum() << std::endl;

            /*
            if (ADCs->GetRMS() > 10.0 || ADCs->GetMaximum() > 0.0 ) {
                std::cout << "Event: " << event << std::endl;
                std::cout << "channel: " << channel << std::endl;
                std::cout << "Peak: " << ADCs->GetMaximum() << std::endl;
                std::cout << "RMS total: " << ADCs->GetRMS() << std::endl;
            }
            */

            double rms = 0;
            if (ADCs->Integral() != 0){
                Double_t par[3];
                double xq = 0.5;
                ADCs->GetQuantiles(1,&par[0],&xq);
                xq = 0.5 - 0.34;
                ADCs->GetQuantiles(1,&par[1],&xq);
                xq = 0.5 + 0.34;
                ADCs->GetQuantiles(1,&par[2],&xq);
                rms = sqrt((pow(par[2]-par[0],2) + pow(par[1]-par[0],2))/2.);
                Double_t mean = par[0];
            }

            TString ADCsName = Form("Event_%i_ADCs_channel%i", event, channel);
            ADCs->SetName(ADCsName);
            //std::cout << "RMS truncated " << rms << std::endl;

            TString frqname = Form("Event_%i_FrequencyWfm_channel%i", event, channel);
            frqspace = (TH1*)waveform->FFT(frqspace, "MAG");
            frqspace->SetName(frqname);
            
            TString hAverageName = Form("Event_%i_FreqAvg_channel%i", event, channel);
            hAverage = average(frqspace, hAverage, 10, nTicks);
            hAverage->SetName(hAverageName);


            f_output.cd(dirFileName);
            waveform->Write();
            frqspace->Write();
            hAverage->Write();
            ADCs->Write();

        }
    }
    f_output.Close();

}
