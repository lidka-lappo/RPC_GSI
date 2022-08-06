#include <TROOT.h>
#include <TApplication.h>
#include <TStopwatch.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TFile.h>
#include <TBranch.h>
#include <R3BRpcHitData.h>
#include <R3BTofdHitData.h>
#include <R3BLosCalData.h>
#include <R3BLosHitData.h>
#include <R3BEventHeader.h>
#include <R3BCalifaCrystalCalData.h>
#include <R3BCalifaHitData.h>
#include <TClonesArray.h>
#include <iostream>

using namespace std;

struct crystal_info {

  Int_t fCrystalId;
  Float_t fTheta;
  Float_t fPhi;

};

#define TRIG(i) (1 << (i - 1))

double fmod_magic(double a, double b) {
         Int_t c = 2048*5;
	 return fmod(a - b +c +c/2.,c) -c/2.;   
}

void rpc_analysis()
{

  TStopwatch timer;
  timer.Start();
 ///////////////////////////////////////////////////////////
 
  TH2F *Pos_histo = new TH2F("Position","Position",1500,0,1500,41,0,42);
  TH2F *Charge = new TH2F("Charge","Charge",1000,-100,100,41,0,42); 
  TH2F *Charge_scint = new TH2F("Charge_scint","Charge_scint",1000,-0,1000,4,0,5); 
  TH1F *Pos_RPC = new TH1F("Pos_RPC","Pos_RPC",1000,0,1500); 
  TH1F *Pos_NB = new TH1F("Pos_NB","Pos_NB",1000,0,3000); 
  TH2F *Time = new TH2F("Time","Time",1000,-4100,4100,41,0,42); 
  TH2F *Time_scint = new TH2F("Time_scint","Time_scint",1000,-4100,4100,4,0,5); 
  TH2F *pos_corr = new TH2F("Pos_corr","Pos_corr",1000,0,1500,1000,0,3000); 
  TH2F *tofd_vs_rpc = new TH2F("Tofd charge vs Rpc tof","Tofd charge vs Rpc tof",1000,-75,-50,500,0,16); 
  TH2F *tofd_Eloss = new TH2F("Tofd charge vs Tofd tof","Tofd charge vs Tofd tof",500,18,40,500,0,16); 
      
  //time Histos
  
  TH2F *time_charge_corr = new TH2F("time rpc - time NB Vs charge NB","time rpc - time NB Vs charge NB",1000,100,260,1000,0,5); 
  TH2F *time_losVsRpc_charge_corr = new TH2F("time RPC - time LOS Vs charge RPC","time RPC - time LOS Vs charge RPC",500,0,5,500,-100,-50); 
  TH2F *time_losVsNB_charge_corr = new TH2F("time NB - time LOS Vs charge NB","time NB - time LOS Vs charge NB",1000,0,400,1000,-25,-10); 

  TH2F *time_charge_corr_charge_RPC = new TH2F("time rpc - time NB Vs charge RPC","time rpc - time NB Vs charge RPC",500,0,5,1000,-55,-35); 
  
  TH2F *time_los_vs_timeNb = new TH2F("time los vs time NB","time Los vs time NB ",1000,-100,-0,1000,-100,-0); 

  TH1F *time_los_rpc= new TH1F("ToF: Strip 21","ToF: Strip 21",1000,-70,-50);
  TH1F *time_los_rpc_cond= new TH1F("ToF: Strip 21 & charge 4 in TOFD","ToF: Strip 21 & charge 4 in TOFD",1000,-70,-50);
  TH2F *tof_RPC_Vs_Pos_RPC= new TH2F("Time RPC - time LOS vs Position RPC","Time RPC - time LOS vs Position RPC",1000,-100,1500,500,-100,-0);
  TH1F *time_diff_histo = new TH1F("Time RPC - time NB","Time RPC - time NB",1000,-30,-10); 
  TH1F *time_los_NB = new TH1F("Time NB - time LOS","Time NB - time LOS",1000,-30,0);

/////////////////////////////////////////////////////////////////////

  TH2F *phiVsPhiHisto = new TH2F("phiVsPhiHisto","P2P : Phi Vs Phi ",300,-190,190,300,-190,190);
  TH2F *energyAngleHisto= new TH2F("energyAngleHisto","Angle Vs Energy",600,0,90,1000,0,1000);
  TH1F *openingAngleHisto = new TH1F("openingAngleHisto","P2P : Opening Angle",300,0,100);
  TH1F *thetaDistribution = new TH1F("thetaDistribution","Theta Distribution",1000,0,90);
  TH1F *coplanarDistribution = new TH1F("coplanarDistribution","Is Coplanar? ",500,-190,190);

  TH2F *calHitsVsOpAngle = new TH2F("calHitsVsOpAngle","Cal Hits Vs Opening Angle",500,-190,190,100,0,100);
  TH2F *thetaVsThetaHisto = new TH2F("thetaVsThetaHisto","Primaries : Theta Vs Theta",600,0,90,300,0,90);

  TH1F *p2p_charges = new TH1F("p2p_charges","P2P: Charges",800,20,100);

/////////////////////////////////////////////////////////////////////

  std::vector<TString> fileList;

  fileList.push_back("/lustre/r3b/mxarepe/unpacked_data/main_TS_0127_0001.root");
  fileList.push_back("/lustre/r3b/mxarepe/unpacked_data/main_TS_0127_0002.root");
  fileList.push_back("/lustre/r3b/mxarepe/unpacked_data/main_TS_0127_0003.root");
  fileList.push_back("/lustre/r3b/mxarepe/unpacked_data/main_TS_0127_0004.root");
  fileList.push_back("/lustre/r3b/mxarepe/unpacked_data/main_TS_0127_0005.root");
  fileList.push_back("/lustre/r3b/mxarepe/unpacked_data/main_TS_0127_0006.root");


 for(int s = 0 ; s < fileList.size() ; s++){
  TFile *eventFile;
  TTree *eventTree;
 
  cout<<"Reading File : "<<fileList.at(s)<<endl;
  cout<<endl;


  eventFile = TFile::Open(fileList.at(s));
  eventTree = (TTree*)eventFile->Get("evt");
  eventTree->SetBranchStatus("*",0);

  eventTree->SetBranchStatus("R3BRpcHitData.*",1);
  eventTree->SetBranchStatus("EventHeader.*",1);
//  eventTree->SetBranchStatus("CalifaHitData*",1);
//  eventTree->SetBranchStatus("CalifaCrystalCalData*",1);
  eventTree->SetBranchStatus("TofdHit.*",1);

  TClonesArray *hitCA = new TClonesArray("R3BRpcHitData");
  TBranch  *hitBranch = eventTree->GetBranch("R3BRpcHitData");
  hitBranch->SetAddress(&hitCA);

  R3BEventHeader *fEventHeader = new R3BEventHeader;
  eventTree->SetBranchAddress("EventHeader.",&fEventHeader);

  TClonesArray *tofdCA = new TClonesArray("R3BTofdHitData");
  TBranch  *tofdhitBranch = eventTree->GetBranch("TofdHit");
  tofdhitBranch->SetAddress(&tofdCA); 
  
/*  TClonesArray *calCalifaCA = new TClonesArray("R3BCalifaCrystalCalData",5);
  TBranch  *calCalifaBranch = eventTree->GetBranch("CalifaCrystalCalData");
  calCalifaBranch->SetAddress(&calCalifaCA);


  TClonesArray *hitCalifaCA = new TClonesArray("R3BCalifaHitData",5);
  TBranch  *hitCalifaBranch = eventTree->GetBranch("CalifaHitData");
  hitCalifaBranch->SetAddress(&hitCalifaCA);*/

  Int_t nEvents = eventTree->GetEntries();
  Int_t nHits, ntofdHits, nCalifaHits, nCalifaCalHits;
  Int_t c = 2048*5;
  Int_t tpatbin=0;	
  double time_Strip=0;
  double pos_Strip;
  double time_diff;
  double pos_scint;
  double charge;
  double energy_los;
  double tofd_Q=0;	
  double tofd_tof=0; 
  double rpc_tof=0; 
  Int_t fTPat=0;	
  Int_t rpc_strip=0;
  Int_t goodEvent,goodCrystalHit;
  Int_t goodP2P = 0;
  Float_t fEnergy,fTheta2,fTheta1,fPhi1,fPhi2,fEnergy1,fEnergy2,openingAngle;

  for(Int_t t = 0; t< nEvents; t++){
  
   hitCA->Clear();
//   hitCalifaCA->Clear();
//   calCalifaCA->Clear();
   tofdCA->Clear();

	 
   eventTree->GetEvent(t);
   nHits = hitCA->GetEntries();
//   nCalifaHits    = hitCalifaCA->GetEntries();
//   nCalifaCalHits = calCalifaCA->GetEntries();
   ntofdHits= tofdCA->GetEntries();

   fEnergy=0.0;
   goodEvent=0;
   goodCrystalHit=0;

   /*
   On-spill + FOOT. 
   TRIG_LMU_OUT( 1) = BEAM_GATE_AUX and not FOOT_DEAD_AUX and in_los_nrolu;
   TRIG_LMU_OUT( 2) = BEAM_GATE_AUX and not FOOT_DEAD_AUX and in_los_nrolu and in_tofd;
   TRIG_LMU_OUT( 3) = BEAM_GATE_AUX and not FOOT_DEAD_AUX and in_los_nrolu and in_califa_and;
   TRIG_LMU_OUT( 4) = BEAM_GATE_AUX and not FOOT_DEAD_AUX and in_los_nrolu and in_califa_and and not in_califa_veto;
   TRIG_LMU_OUT( 5) = BEAM_GATE_AUX and not FOOT_DEAD_AUX and in_los_nrolu and in_califa_or;
   TRIG_LMU_OUT( 6) = BEAM_GATE_AUX and not FOOT_DEAD_AUX and in_los_nrolu and in_califa_or and not in_califa_veto;
   TRIG_LMU_OUT( 6) = BEAM_GATE_AUX and not FOOT_DEAD_AUX and in_los_nrolu and in_neuland;
  
   On-spill - FOOT. 
   TRIG_LMU_OUT( 7) = BEAM_GATE_AUX and NONFOOT_BONUS_AUX and in_los_nrolu;
   TRIG_LMU_OUT( 8) = BEAM_GATE_AUX and NONFOOT_BONUS_AUX and in_los_nrolu and in_tofd;
   TRIG_LMU_OUT( 9) = BEAM_GATE_AUX and NONFOOT_BONUS_AUX and in_los_nrolu and in_califa_and;
   TRIG_LMU_OUT(10) = BEAM_GATE_AUX and NONFOOT_BONUS_AUX and in_los_nrolu and in_califa_and and not in_califa_veto;
   TRIG_LMU_OUT(11) = BEAM_GATE_AUX and NONFOOT_BONUS_AUX and in_los_nrolu and in_califa_or;
   // TRIG_LMU_OUT(12) = BEAM_GATE_AUX and NONFOOT_BONUS_AUX and in_los_nrolu and in_califa_or and not in_califa_veto;
   TRIG_LMU_OUT(12) = BEAM_GATE_AUX and NONFOOT_BONUS_AUX and in_los_nrolu and in_neuland;
  
   Off-spill. 
   TRIG_LMU_OUT(13) = not BEAM_GATE_AUX and in_califa_or;
   TRIG_LMU_OUT(14) = not BEAM_GATE_AUX and in_neuland;
   TRIG_LMU_OUT(15) = not BEAM_GATE_AUX and in_tofd;
   TRIG_LMU_OUT(16) = not BEAM_GATE_AUX and in_rpc; 
   */
   bool p2p = (TRIG(4) | TRIG(10)) & fEventHeader->GetTpat();
   //if (fEventHeader->GetTpat() & 0x1f){
   if (p2p){
    
    if(nHits == 0){continue;}
   
  /*  for(Int_t z = 0; z < nCalifaCalHits; z++){

     auto hitCalifaCal = (R3BCalifaCrystalCalData*)calCalifaCA->At(z);
     if (0.001*(hitCalifaCal->GetEnergy()) >= 100)
         goodEvent++;
         if(isnan(0.001*(hitCalifaCal->GetEnergy()))){
         continue;
         }
    }

     if(nCalifaHits == 2  &&  goodEvent>=2){
         fPhi1 = TMath::RadToDeg()*((R3BCalifaHitData*)hitCalifaCA->At(0))->GetPhi();
         fPhi2 = TMath::RadToDeg()*((R3BCalifaHitData*)hitCalifaCA->At(1))->GetPhi();

         fTheta1 = TMath::RadToDeg()*((R3BCalifaHitData*)hitCalifaCA->At(0))->GetTheta();
         fTheta2 = TMath::RadToDeg()*((R3BCalifaHitData*)hitCalifaCA->At(1))->GetTheta();

         for (Int_t n =0;n<nCalifaHits;n++)
            fEnergy= fEnergy + 0.001*((R3BCalifaHitData*)hitCalifaCA->At(n))->GetEnergy();

         fEnergy1 = 0.001*((R3BCalifaHitData*)hitCalifaCA->At(0))->GetEnergy();
         fEnergy2 = 0.001*((R3BCalifaHitData*)hitCalifaCA->At(1))->GetEnergy();


         openingAngle = TMath::Sin(TMath::DegToRad()*fTheta1)*TMath::Sin(TMath::DegToRad()*fTheta2)*TMath::Cos(TMath::DegToRad()*fPhi2-TMath::DegToRad()*fPhi1)\
		       	+ TMath::Cos(TMath::DegToRad()*fTheta1)*TMath::Cos(TMath::DegToRad()*fTheta2);
         openingAngle = TMath::RadToDeg()*TMath::ACos(openingAngle);

        if(fEnergy < 560 && TMath::Abs(fPhi1-fPhi2)>165 && TMath::Abs(fPhi1-fPhi2)<195 && nCalifaCalHits<45){


         openingAngleHisto->Fill(openingAngle);

         phiVsPhiHisto->Fill(fPhi1,fPhi2);
         energyAngleHisto->Fill(fTheta1,fEnergy1);
         energyAngleHisto->Fill(fTheta2,fEnergy2);

         thetaVsThetaHisto->Fill(fTheta1,fTheta2);
         goodP2P++;


         }
       }
*/
    for(Int_t i = 0; i< ntofdHits; i++){
     auto hitTofd = (R3BTofdHitData*)tofdCA->At(i);
	tofd_Q = hitTofd->GetEloss();
	tofd_tof = hitTofd->GetTof();
        tofd_Eloss->Fill(tofd_tof,tofd_Q);

    }
 // cout << "HITS!!!!!!!! " <<  nHits << endl;  
    for(Int_t i = 0; i< nHits; i++){
     auto map1 = (R3BRpcHitData*)(hitCA->At(i));
     //RPC data
//cout << " DET ID !!!!!!!! " <<map1->GetDetId() << endl;
     if(map1->GetDetId()==0){
      Pos_histo->Fill(map1->GetPos(),map1->GetChannelId());
      Time->Fill(map1->GetTime(),map1->GetChannelId());
      Charge->Fill(map1->GetCharge(),map1->GetChannelId());
      rpc_strip = map1->GetChannelId();
      if(rpc_strip == 21){
       time_Strip=map1->GetTime();
       charge = map1->GetCharge();
       pos_Strip=map1->GetPos();
       if(tofd_Q < 4.5 &&  tofd_Q > 3.5){
        time_los_rpc->Fill(map1->GetTof());
        time_los_rpc_cond->Fill(map1->GetTof());
       }
       else{	
        time_los_rpc->Fill(map1->GetTof());
       }
       time_losVsRpc_charge_corr->Fill(charge, map1->GetTof());
       tofd_vs_rpc->Fill(map1->GetTof(),tofd_Q);
       tof_RPC_Vs_Pos_RPC->Fill(pos_Strip,map1->GetTof());
       rpc_tof = map1->GetTof();
      }
     }

     //Pmt data 
     if(map1->GetDetId()==1){
      Time_scint->Fill(map1->GetTime(),map1->GetChannelId());
      Charge_scint->Fill(map1->GetCharge(),map1->GetChannelId());
      if(map1->GetChannelId() == 2){
       pos_scint=map1->GetPos();
       if(time_Strip != 0 && rpc_strip == 21){
        time_diff = time_Strip - map1->GetTime() ;
        time_charge_corr->Fill(map1->GetCharge(),time_diff);
        time_charge_corr_charge_RPC->Fill(charge,time_diff);
        time_diff_histo->Fill(time_diff);
        Pos_RPC->Fill(pos_Strip);
        Pos_NB->Fill(pos_scint);
        pos_corr->Fill(pos_Strip,pos_scint);
       }
       time_losVsNB_charge_corr->Fill(map1->GetCharge(),map1->GetTof());
       time_los_NB->Fill(map1->GetTof());
       time_los_vs_timeNb->Fill(map1->GetTof(),rpc_tof);
      }
     }
    }
   }
  }
 }
/*  TCanvas *myFirstCanvas = new TCanvas("myFirstCanvas","myFirstCanvas",2000,2000);
   myFirstCanvas->Divide(2,2);


   myFirstCanvas->cd(1);
   openingAngleHisto->Draw();

   myFirstCanvas->cd(2);
   phiVsPhiHisto->Draw("colz");

   myFirstCanvas->cd(3);
   thetaVsThetaHisto->Draw("colz");

   myFirstCanvas->cd(4);
   energyAngleHisto->Draw("colz");*/

  TCanvas *C1 = new TCanvas("C1","C1",600,800);
  TCanvas *C2 = new TCanvas("C2","C2",600,800);
  TCanvas *C3 = new TCanvas("C3","C3",600,800);
  TCanvas *C4 = new TCanvas("C4","C4",600,800);
  TCanvas *C5 = new TCanvas("C5","C5",600,800);
  TCanvas *C6 = new TCanvas("C6","C6",600,800);
  TCanvas *C8 = new TCanvas("C8","C8",600,800);
  TCanvas *C9 = new TCanvas("C9","C9",600,800);
  TCanvas *C10 = new TCanvas("C10","C10",600,800);
  TCanvas *C13 = new TCanvas("C13","C13",600,800);

  time_los_rpc->GetXaxis()->SetTitle("ToF (ns)");
  time_los_rpc->GetYaxis()->SetTitle("Counts");
		        
  time_los_rpc_cond->GetXaxis()->SetTitle("ToF (ns)");
  time_los_rpc_cond->GetYaxis()->SetTitle("Counts");

  C2->Divide(2,1);
  C4->Divide(2,1);
  C5->Divide(2,1);
//  C8->Divide(3,1);
  C9->Divide(3,1);
  C13->Divide(2,1);
  C1->cd();
  //Pos_histo->Draw("colz");
  time_los_vs_timeNb->Draw("colz");
  C2->cd(1);
  Time->Draw("colz");
  C2->cd(2);
  Time_scint->Draw("colz");
  C3->cd();
  tofd_Eloss->Draw("colz");
  //charge maps
  C4->cd(1);
  Charge->Draw("colz");
  C4->cd(2);
  Charge_scint->Draw("colz");

  C6->cd();
  tofd_vs_rpc->Draw("colz");
  //time difs
  //C8->cd(1);
  //time_diff_histo->Draw("colz");
  //time_los_rpc_cond->Draw("hist"); 
  C8->cd();
  time_los_rpc->Draw("hist");
  //C8->cd(3);
  //time_los_NB->Draw("hist");
  //time and charge correlation
  C5->cd(1);
  time_charge_corr->Draw("colz");
  C5->cd(2);
  time_charge_corr_charge_RPC->Draw("colz");

  C10->cd();
  time_losVsNB_charge_corr->Draw("colz");
 
  C13->cd(1);
  time_losVsRpc_charge_corr->Draw("colz");  
  C13->cd(2);
  tof_RPC_Vs_Pos_RPC->Draw("colz");  

  //pos RPC vs pos NB
  C9->cd(1);
  pos_corr->Draw("colz");
  C9->cd(2);
  Pos_RPC->Draw("colz");
  C9->cd(3);
  Pos_NB->Draw("colz");

  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  std::cout << std::endl << std::endl;
  std::cout << "Macro finished succesfully." << std::endl;
  std::cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << std::endl << std::endl;

}


# ifndef __CINT__
int main(int argc, char **argv){
	TApplication app("app",&argc,argv);
	rpc_analysis();
	app.Run();
	return 0 ;
}
#endif
