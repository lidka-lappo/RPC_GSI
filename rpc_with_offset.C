using namespace std;
struct crystal_info {

  Int_t fCrystalId;
  Float_t fTheta;
  Float_t fPhi;

};
int nbin = 2500;
int minData = -500;
int maxData = 2000;
int minHist = 0;
int maxHist = 1500;
int offset[41];

void moveHist(TH1F *pos, int shift)
{
	TAxis *a = pos->GetXaxis();
	a->Set(a->GetNbins(), (a->GetXmin()+shift), (a->GetXmax()+shift));
}

void offset_from_file()
{
	ifstream myfile;
	myfile.open("offset.txt", ios::in|ios::out);
	for(int i=0; i<41; i++)
	{
		int tmp;
		myfile>>tmp;
		offset[i]=tmp;
		cout<<offset[i]<<endl;;		
	}
	myfile.close();
}
double fmod_magic(double a, double b) {
         Int_t c = 2048*5;
	 return fmod(a - b +c +c/2.,c) -c/2.;   
}
void rpc_with_offset()
{

  TStopwatch timer;
  timer.Start();
 ///////////////////////////////////////////////////////////
offset_from_file();
  TH2F *Pos_histo = new TH2F("Position","Position",nbin,minData,maxData,41,0,42);


	

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
      int strip = map1->GetChannelId();
      Pos_histo->Fill(offset[strip-1]+map1->GetPos(),strip);
     }

    }
   }
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
	rpc_with_offset();
	app.Run();
	return 0 ;
}
#endif
