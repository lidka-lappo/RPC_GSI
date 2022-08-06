using namespace std;

struct crystal_info {

  Int_t fCrystalId;
  Float_t fTheta;
  Float_t fPhi;

};


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
  TH1F *strip_histo[42];   

  TH1F *strip_smooth[42];   

  for(int i =0; i <42; i ++)
  {
   	stringstream strs;
	strs << i;
	string tmp = strs.str();
	char *name = (char *) tmp.c_str();
	strip_histo[i] = new TH1F(name,name,150,0,1500);
	strip_smooth[i] = new TH1F(name,name,300,0,1500);

  }

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

  TClonesArray *hitCA = new TClonesArray("R3BRpcHitData");
  TBranch  *hitBranch = eventTree->GetBranch("R3BRpcHitData");
  hitBranch->SetAddress(&hitCA);

  R3BEventHeader *fEventHeader = new R3BEventHeader;
  eventTree->SetBranchAddress("EventHeader.",&fEventHeader);

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
if(t%100000==0)
	cout<<t<<endl;	 
   eventTree->GetEvent(t);
   nHits = hitCA->GetEntries();
   fEnergy=0.0;
   goodEvent=0;
   goodCrystalHit=0;
    if(nHits == 0){continue;}
   
    for(Int_t i = 0; i< nHits; i++){
     auto map1 = (R3BRpcHitData*)(hitCA->At(i));
     if(map1->GetDetId()==0){
      Pos_histo->Fill(map1->GetPos(),map1->GetChannelId());
      strip_histo[map1->GetChannelId()]->Fill(map1->GetPos());
     }

    }
   }
 }
cout<<"kupa"<<endl;
int tmp1, tmp2, tmp3, tmp4, tmp5;
/*for(int i=0; i <42; i++)
{
	cout<<"Kupa "<<i<<endl;
	for(int j=0; j<1500; j++)
	{
	tmp1 = strip_histo[i]->GetBinContent(j-2);
	
	tmp2 = strip_histo[i]->GetBinContent(j-1);
	
	tmp3 = strip_histo[i]->GetBinContent(j);
	
	tmp4 = strip_histo[i]->GetBinContent(j+1);
	
	tmp5 = strip_histo[i]->GetBinContent(j+2);
	strip_smooth[i]->Fill((tmp1+tmp2+tmp3+tmp4+tmp5)/5);
//	cout<<j;
	}

}*/
  TCanvas *C1 = new TCanvas("C1","C1",600,800);
  int n =1;
  TCanvas *C[n];
  for (int i=0; i<n; i++)
  {
   	stringstream strs;
	strs << i;
	string tmp = strs.str();
	char *name = (char *) tmp.c_str();
	C[i] = new TCanvas(name,name,1500,600);
	C[i]->Divide(3,2);
  }
  C1->cd();
  Pos_histo->Draw("colz");
  
for(int i=0; i<n; i++)
{
	for(int j=1; j<7; j++)
	{
		C[0]->cd(j);
		strip_histo[i*6+j]->Draw("");
	}
}

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
