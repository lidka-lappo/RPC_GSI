using namespace std;
struct crystal_info {

  Int_t fCrystalId;
  Float_t fTheta;
  Float_t fPhi;

};
int nbin = 1500; //100
int minData = 0;
int maxData = 1500;
int minHist = 0;
int maxHist = 1500;
TH1F *strip_histo_low[41];   
TF1 *fit[42];
TF1 *fitAll[42];
TH1F *hist_sub[42];
int offset[41];

void calc_pik(TH1F **strip_histo, TString fileList)
{
	int newnBin = 47;
	int newMinHist = 795;
int posY=0;
int posX=0;
double tmpMax=0;
for(int j =1; j<42; j++)
{
strip_histo[j]->Scale(1.0/strip_histo[j]->Integral());
//int tmp=strip_histo_low[j]->ShowPeaks(0.1, "nobackground", 0.7); //0.1-0.001
//int tmp=strip_histo[j]->ShowPeaks(0.1, "", 0.5); //0.1-0.001
//cout<<j<<":"<<tmp<<endl;
if(j>3){
	if(j<30)
	{
		if(tmpMax<strip_histo[j]->GetMaximum())
		{
			posY=j;
			//TList *fun = Pos_histo->GetListOfFunctions();
			//TPolyMarker *pm = (TPolyMarker*)fun->FindObject("TPolyMarker");
			//cout<<"(x:y) ("<<pm[0].GetX()<<" : "<<pm[0].GetY()<<")"<<endl;
			posX=strip_histo[j]->GetMaximumBin();
			tmpMax=strip_histo[j]->GetMaximum();
			cout<<" ("<<posX<<" : "<<posY*30<<")"<<endl;

		}
	}
}
}
fstream myfile1;
myfile1.open("dataAll.txt", ios::app);
 myfile1<<"File : "<<fileList<<" ("<<posX<<" : "<<posY*30<<")"<<endl;
myfile1<<"strip: "<< posY <<endl;
myfile1.close();
cout<<"(x:y) ("<<posX<<" : "<<posY<<")"<<endl;
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
		cout<<offset[i]<<endl;
	}
	myfile.close();
}

Double_t fitGaussReject(Double_t *x, Double_t *par)
{
 Double_t norm = par[0];
 Double_t mpv =par[1];
 Double_t mean =par[2];
 Double_t sigma = par[3];
if( x[0]>=800)
 {
  return norm*ROOT::Math::gaussian_pdf(x[0], sigma, mean);
// return norm*TMath::Landau(x[0], mpv, sigma, false);
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
void rpc_pik_searching()
{

  TStopwatch timer;
  timer.Start();
 ///////////////////////////////////////////////////////////
 
  offset_from_file();
  TH2F *Pos_histo = new TH2F("Position","Position",nbin,minData,maxData,41,0,42);

//  TGraph *strips_[42];   



/////////////////////////////////////////////////////////////////////


  std::vector<TString> fileList;
//background
 //fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22045050208.root");
 //fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046083458.root");

//low
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046050208.root");
//
//short
 fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046060515.root");
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046064454.root");
//measures
//down
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046082416.root");
 //fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046060952.root");
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046065513.root");
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046070926.root");
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046071959.root");
//
//

//up
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046081301.root");
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046080126.root");
//  fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046075139.root");
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046074125.root");
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046073033.root");
 for(int s = 0 ; s < fileList.size() ; s++){
  TFile *eventFile;
  TTree *eventTree;
 
TH1F *strip_histo[42];   
  for(int i =0; i <42; i ++)
  {
   	stringstream strs;
   	stringstream strs1;
	strs << i;
	strs1 << i <<"low";
	string tmp = strs.str();
	string tmp1 = strs1.str();
	char *name = (char *) tmp.c_str();
	char *name1 = (char *) tmp1.c_str();
	strip_histo[i] = new TH1F(name,name,nbin,minHist,maxHist);
	strip_histo_low[i] = new TH1F(name1,name1,100,minHist,maxHist);
	
  }

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
      strip_histo[strip]->Fill(offset[strip-1]+map1->GetPos());
   //   strip_histo_low[strip]->Fill(offset[strip-1]+map1->GetPos());
     }

    }
   }

calc_pik(strip_histo, fileList.at(s));
}
//Int_t pik = Pos_histo->ShowPeaks(10,"",0.05);
//cout<<"pik "<<pik<<endl;
//TList *fun = Pos_histo->GetListOfFunctions();
//if(pik!=0)
//{
//TPolyMarker *pm = (TPolyMarker*)fun->FindObject("TPolyMarker");
//cout<<"(x:y) ("<<pm[0].GetX()<<" : "<<pm[0].GetY()<<")"<<endl;
//} 
TCanvas *C1 = new TCanvas("C1","C1",600,800);
/*  int n= 7;
 TCanvas *C[n];
  for (int i=0; i<n; i++)
  {
   	stringstream strs;
	strs << i;
	string tmp = strs.str();
	char *name = (char *) tmp.c_str();
	C[i] = new TCanvas(name,name,1500,600);
	C[i]->Divide(3,2);
  }*/
 C1->cd();
 Pos_histo->Draw("COLZ");
//C[0]->cd(1);
/*
for(int i=0; i<n; i++)
{
	for(int j=0; j<6; j++)
	{
		if(i!=6||j!=5){
			C[i]->cd(j+1);
			strip_histo[i*6+j+1]->Draw();
			//strip_histo_low[i*6+j+1]->Draw();
		}
	}
}*/


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
	rpc_pik_searching();
	app.Run();
	return 0 ;
}
#endif
