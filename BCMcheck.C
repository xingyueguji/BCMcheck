void getchargeandtime(int runNum, Double_t &total_charge_bcm1_cut,Double_t &total_charge_bcm2_cut,
Double_t &total_charge_bcm4a_cut,Double_t &total_charge_bcm4c_cut,Double_t &total_time_bcm1_cut,
Double_t &total_time_bcm2_cut,Double_t &total_time_bcm4a_cut,Double_t &total_time_bcm4c_cut,Double_t &total_time)
{
  //Created and Read Rootfile
  
  TString filename=Form("../boiling/ROOTfiles/NPS/SCALERS/nps_replay_scalers_%d_1_-1.root", runNum); 

  TFile *file = new TFile(filename);
  
  TTree *tdata = (TTree*) file->Get("T");    //data TTree
  TTree *tscal = (TTree*) file->Get("TSH");  //Scaler TTree


  //Get Scaler Leafs
  Double_t Scal_evNum;    //event number associated with scaler reads
  tscal->SetBranchAddress("evNumber", &Scal_evNum);
  Double_t  Scal_BCM4A_charge;
  tscal->SetBranchAddress("H.BCM4A.scalerCharge",&Scal_BCM4A_charge);
  Double_t  Scal_BCM4A_current;
  tscal->SetBranchAddress("H.BCM4A.scalerCurrent",&Scal_BCM4A_current);  
  Double_t  Scal_BCM4C_charge;
  tscal->SetBranchAddress("H.BCM4C.scalerCharge",&Scal_BCM4C_charge);
  Double_t  Scal_BCM4C_current;
  tscal->SetBranchAddress("H.BCM4C.scalerCurrent",&Scal_BCM4C_current);
  Double_t  Scal_BCM1_charge;
  tscal->SetBranchAddress("H.BCM1.scalerCharge",&Scal_BCM1_charge);
  Double_t  Scal_BCM1_current;
  tscal->SetBranchAddress("H.BCM1.scalerCurrent",&Scal_BCM1_current); 
  Double_t  Scal_BCM1;
  tscal->SetBranchAddress("H.BCM1.scaler",&Scal_BCM1);
  Double_t Scal_BCM4A;
  tscal->SetBranchAddress("H.BCM4A.scaler",&Scal_BCM4A);
  Double_t  Scal_BCM2_charge;
  tscal->SetBranchAddress("H.BCM2.scalerCharge",&Scal_BCM2_charge);
  Double_t  Scal_BCM2_current;
  tscal->SetBranchAddress("H.BCM2.scalerCurrent",&Scal_BCM2_current); 

  Double_t  Scal_time;
  tscal->SetBranchAddress("H.1MHz.scalerTime",&Scal_time);
  

  //Defive Quantities To Store Previous Reads and cumulative quantities
  Double_t prev_time = 0.;
  Double_t prev_charge_bcm4a = 0.;
  Double_t prev_charge_bcm4c = 0.;
  Double_t prev_charge_bcm1 = 0.;
  Double_t prev_charge_bcm2 = 0.;
  Double_t prev_scaler_bcm1 = 0.;
  Double_t prev_scaler_bcm4a =0;

  Double_t total_charge_bcm4a = 0.;
  Double_t total_charge_bcm4c = 0.;
  Double_t total_charge_bcm1 = 0.;
  Double_t total_charge_bcm2 = 0.;


  //Loop Over Scaler Reads 
  Long64_t scal_entries = tscal->GetEntries();

  Int_t evt_flag_bcm4a[scal_entries];             //Store Flag [0 or 1], to know if scaler read passed current cut (1) or not (0)
  Int_t evt_flag_bcm4c[scal_entries];             //Store Flag [0 or 1], to know if scaler read passed current cut (1) or not (0)
  Int_t evt_flag_bcm1[scal_entries];
  Int_t evt_flag_bcm2[scal_entries]; 

  Int_t scal_evt_num[scal_entries];    //Store Event Associated with Scaler Read
  //cout << "Scaler events"<<scal_entries<<endl; 
  for (int i = 0; i < scal_entries; i++) {
    
    //**NOTE: Each scaler read is associated with as specific event number
    //        as (scaler read 1-> event 1000,  scaler read 2 -> event 2300, ...)
    //        This means events up to 1000 correspond to scaler read 1, ...
   
    tscal->GetEntry(i);
    //Save all no cut quantities.
    total_time = Scal_time;
    total_charge_bcm4a = Scal_BCM4A_charge;
    total_charge_bcm4c = Scal_BCM4C_charge;
    total_charge_bcm1 = Scal_BCM1_charge;
    total_charge_bcm2 = Scal_BCM2_charge;
    //Apply cut
    if(Scal_BCM1_current > 2){
      total_time_bcm1_cut = total_time_bcm1_cut + (Scal_time - prev_time);
	    total_charge_bcm1_cut = total_charge_bcm1_cut + (Scal_BCM1_charge - prev_charge_bcm1);
    }
    if(Scal_BCM2_current > 2){
      total_time_bcm2_cut = total_time_bcm2_cut + (Scal_time - prev_time);
	    total_charge_bcm2_cut = total_charge_bcm2_cut + (Scal_BCM2_charge - prev_charge_bcm2);
    }
    if(Scal_BCM4A_current > 2){
      total_time_bcm4a_cut = total_time_bcm4a_cut + (Scal_time - prev_time);
	    total_charge_bcm4a_cut = total_charge_bcm4a_cut + (Scal_BCM4A_charge - prev_charge_bcm4a);
    }
    if(Scal_BCM4C_current > 2){
      total_time_bcm4c_cut = total_time_bcm4c_cut + (Scal_time - prev_time);
	    total_charge_bcm4c_cut = total_charge_bcm4c_cut + (Scal_BCM4C_charge - prev_charge_bcm4c);
    }
    prev_time = Scal_time;
    prev_charge_bcm4a = Scal_BCM4A_charge;
    prev_charge_bcm4c = Scal_BCM4C_charge;
    prev_charge_bcm1 = Scal_BCM1_charge;
    prev_charge_bcm2 = Scal_BCM2_charge;


  }
}

void BCMcheck(){
  //open the runlist file
  Int_t runNUM;
  string line;
  TString filename = "BCM.dat";
  ifstream ifs;
  ifs.open(filename);

  vector <double> bcm1ratio;
  vector <double> bcm2ratio;
  vector <double> bcm4cratio;
  vector <double> abovetime1;
  vector <double> abovetime2;
  vector <double> abovetime4c;


  TCanvas *c1 = new TCanvas("","",800,600);
  vector<TGraph*> graphArray;
  TGraph *g1;
  TGraph *g2;
  TGraph *g3;

  while (getline(ifs, line)){
    //convert run from string to int
    runNUM = stoi(line);
    cout << runNUM << endl;

  Double_t total_charge_bcm1_cut = 0;
  Double_t total_charge_bcm2_cut = 0;
  Double_t total_charge_bcm4a_cut = 0;
  Double_t total_charge_bcm4c_cut = 0;
  Double_t total_time_bcm1_cut = 0;
  Double_t total_time_bcm2_cut = 0;
  Double_t total_time_bcm4a_cut = 0;
  Double_t total_time_bcm4c_cut = 0;
  Double_t total_time = 0;

    getchargeandtime(runNUM,total_charge_bcm1_cut,total_charge_bcm2_cut,total_charge_bcm4a_cut,total_charge_bcm4c_cut,
    total_time_bcm1_cut,total_time_bcm2_cut,total_time_bcm4a_cut,total_time_bcm4c_cut,total_time);
    bcm1ratio.push_back(total_charge_bcm1_cut/total_charge_bcm4a_cut);
    bcm2ratio.push_back(total_charge_bcm2_cut/total_charge_bcm4a_cut);
    bcm4cratio.push_back(total_charge_bcm4c_cut/total_charge_bcm4a_cut);
    abovetime1.push_back((total_time_bcm1_cut/total_time)*100);
    abovetime2.push_back((total_time_bcm2_cut/total_time)*100);
    abovetime4c.push_back((total_time_bcm4c_cut/total_time)*100);

    g1 = new TGraph(bcm1ratio.size(),&abovetime1[0],&bcm1ratio[0]);
    g1->SetMarkerColor(1);
    g1->SetMarkerStyle(23);
    g1->SetMarkerSize(1);
    g2 = new TGraph(bcm2ratio.size(),&abovetime2[0],&bcm2ratio[0]);
    g2->SetMarkerColor(2);
    g2->SetMarkerStyle(21);
    g2->SetMarkerSize(1);
    g3 = new TGraph(bcm4cratio.size(),&abovetime4c[0],&bcm4cratio[0]);
    g3->SetMarkerColor(3);
    g3->SetMarkerStyle(22);
    g3->SetMarkerSize(1);        

    graphArray.push_back(g1);
    graphArray.push_back(g2);
    graphArray.push_back(g3);
    bcm1ratio.clear();
    bcm2ratio.clear();
    bcm4cratio.clear();
    abovetime1.clear();
    abovetime1.clear();
    abovetime4c.clear();
  }



for (int i = 0; i < graphArray.size(); i++) {
    graphArray[i]->GetYaxis()->SetRangeUser(0.9, 1.02);
    c1->cd();
    if (i == 0) {
        graphArray[i]->GetXaxis()->SetLimits(0, 100);
        graphArray[i]->Draw("AP");  // Draw the first graph without "SAME"
    } else {
        graphArray[i]->Draw("SP");  // Draw subsequent graphs with "SAME"
    }
}
    TLegend *leg=new TLegend(.1,.15,.5,.4);
    leg->AddEntry(graphArray[0],"BCM1/BCM4A","p");
    leg->AddEntry(graphArray[1],"BCM2/BCM4A","p");
    leg->AddEntry(graphArray[2],"BCM4C/BCM4A","p");
    leg->Draw();

    string pdfname;
    c1->Update();
    pdfname = to_string(runNUM) + ".pdf";
    c1->SaveAs("Plot.pdf");


}