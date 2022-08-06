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
      
/////////////////////////////////////////////////////////////////////

  std::vector<TString> fileList;

 // fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22045050208.root");

//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046050208.root");
  fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046060515.root");
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046060952.root");
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046064454.root");
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st2204607926.root");
 //fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046074125.root");
  //fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046073033.root");

//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046071959.root");
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046080126.root");
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046081301.root");
 //fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046082416.root");
  //fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046083458.root");
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
 // eventTree->SetBranchStatus("TofdHit.*",1);

  TClonesArray *hitCA = new TClonesArray("R3BRpcHitData");
  TBranch  *hitBranch = eventTree->GetBranch("R3BRpcHitData");
  hitBranch->SetAddress(&hitCA);

  R3BEventHeader *fEventHeader = new R3BEventHeader;
  eventTree->SetBranchAddress("EventHeader.",&fEventHeader);

  //TClonesArray *tofdCA = new TClonesArray("R3BTofdHitData");
  //TBranch  *tofdhitBranch = eventTree->GetBranch("TofdHit");
  //tofdhitBranch->SetAddress(&tofdCA); 
  
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
//  cout<<"event "<<t<<endl;
   hitCA->Clear();
  // tofdCA->Clear();
if(t%10000==0)
	cout<<t<<endl;	 
   eventTree->GetEvent(t);
   nHits = hitCA->GetEntries();
   //ntofdHits= tofdCA->GetEntries();
//cout<<"Number of Hits: "<<nHits<<endl;
   fEnergy=0.0;
   goodEvent=0;
   goodCrystalHit=0;
//What is p2p?
//  bool p2p = (TRIG(4) | TRIG(10)) & fEventHeader->GetTpat();
  //	cout<<TRIG(4)<<endl;
  // if (p2p){
    
//	cout<<"And here not"<<endl;
    if(nHits == 0){continue;}
   
  //  for(Int_t i = 0; i< ntofdHits; i++){
   //  auto hitTofd = (R3BTofdHitData*)tofdCA->At(i);
//	tofd_Q = hitTofd->GetEloss();
//	tofd_tof = hitTofd->GetTof();

  //  
  //  }
    for(Int_t i = 0; i< nHits; i++){
     auto map1 = (R3BRpcHitData*)(hitCA->At(i));
  //    cout<<"pozycja "<< map1->GetPos()<<endl;
     if(map1->GetDetId()==0){
      Pos_histo->Fill(map1->GetPos(),map1->GetChannelId());
     }

    }
   }
  //}
 }

  TCanvas *C1 = new TCanvas("C1","C1",600,800);


  C1->cd();
  Pos_histo->Draw("colz");
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
