using namespace std;
struct crystal_info {

  Int_t fCrystalId;
  Float_t fTheta;
  Float_t fPhi;

};

TH1F *strip_histo[42];   
TF1 *fit[42];
TF1 *fitAll[42];
TH1F *hist_sub[42];
Double_t fitGaussReject(Double_t *x, Double_t *par)
{
 Double_t norm = par[0];
 Double_t mpv =par[1];
 Double_t mean =par[2];
 Double_t sigma = par[3];
 if(x[0]<=250 ||  x[0]>=750)
 {
  return norm*ROOT::Math::gaussian_pdf(x[0], sigma, mean);
 //return norm*TMath::Landau(x[0], mpv, sigma, false);
 }
 else
 {
  TF1::RejectPoint();
  return 0;
 }

}
Double_t fitGaussAll(Double_t *x, Double_t *par)
{
 Double_t norm = par[0];
 Double_t mpv = par[1];
 Double_t mean =par[2];
 Double_t sigma = par[3];
 return norm*ROOT::Math::gaussian_pdf(x[0], sigma, mean);
// return norm*TMath::Landau(x[0], mpv, sigma, false);
}

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

//  TGraph *strips_[42];   


  for(int i =0; i <42; i ++)
  {
   	stringstream strs;
	strs << i;
	string tmp = strs.str();
	char *name = (char *) tmp.c_str();
	strip_histo[i] = new TH1F(name,name,1500,0,1500);
	
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
  TCanvas *C1 = new TCanvas("C1","C1",600,800);
  int n =7;
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
for(int j =0; j<42; j++)
{
   	stringstream strs1, strs2;
	strs1 <<"fit"<< j;
	strs2 <<"fitAll"<< j;
	string tmp11 = strs1.str();
	string tmp22 = strs2.str();
	char *name1 = (char *) tmp11.c_str();
	char *name2 = (char *) tmp22.c_str();
	fit[j] = new TF1(name1, "fitGaussReject", 0.0, 1500.0, 4);
	fit[j]->SetParameters(4000,strip_histo[j]->GetMaximum(), strip_histo[j]->GetMean(), strip_histo[j]->GetRMS());
	strip_histo[j]->Fit(name1, "RQN");

	fitAll[j] = new TF1(name2, "fitGaussAll", 0.0, 1500.0, 4);
	fitAll[j]->SetParameters(fit[j]->GetParameters());
	//strip_histo[4]->Fit("landau");
	//strip_histo[4]->Draw("");
	//fitAll->Draw("SAME");
	//fit->Draw("SAME");
	//  Pos_histo->Draw("colz");

	hist_sub[j]= new TH1F(name1, "h1", 1500, 0, 1500);
	for(int i =0; i<1500; i++)
	{
		hist_sub[j]->SetBinContent(i+1, (fitAll[j]->Eval(i)-strip_histo[j]->GetBinContent(i+1)));
	}
	cout<<"max: "<<hist_sub[j]->GetMaximum()<<endl;;
	//h1->Draw();
}
for(int i=0; i<n; i++)
{
	for(int j=0; j<6; j++)
	{
		C[i]->cd(j+1);
		strip_histo[i*6+j+1]->Draw();
	//	hist_sub[i*6+j+1]->Draw("SAME");
	 	fitAll[i*6+j+1]->Draw("SAME");
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
