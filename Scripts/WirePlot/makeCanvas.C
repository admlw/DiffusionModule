void makeCanvas(){

    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 500);
    c1->cd();
    c1->Divide(4,2);

    for (int i = 4128; i <= 4136; i++){
        c1->cd(i-4127);

        TH1D* h = (TH1D*)_file0->Get(Form("Event_4301/Event_4301_TimeWfm_channel%i", i));
        h->GetYaxis()->SetRangeUser(2030, 2060);
            h->Draw();
    }
}
