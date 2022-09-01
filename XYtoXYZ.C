using namespace std;

int nbin = 1500;
int minData =0;
int maxData = 1500;
int minHist = 0;
int maxHist = 1500;
int offset[41];
Double_t Xx[9] ={};
Double_t Yy[9] ={};
Double_t Zz[9] ={};
Double_t par[3] ={};


void CalculatingPoints(Double_t *posXYZ, Double_t Xs, Double_t Ys)
{
//Calculating the distance to point

double l[4]={0, 0, 0, 0};
//distance of cardinal point measured from frame of RPC
Double_t XR[4] ={1396, 780, 179, 180};
Double_t YR[4] ={1187, 1189, 1146, 696};
for(int i=0; i<4; i++)
{
	 l[i]=sqrt((Xs-XR[i])*(Xs-XR[i])+(Ys-YR[i])*(Ys-YR[i]));
}

//algorithm of calculating x, y, z of any point knowing its the distance from 4 cardinal points 
double A, B, C, D, E, F;

A= 2*(Xx[0]-Xx[1]+par[0]*(Zz[0]-Zz[1]));
B= 2*(Yy[0]-Yy[1]+par[1]*(Zz[0]-Zz[1]));
C=Xx[1]*Xx[1]-Xx[0]*Xx[0]+Yy[1]*Yy[1]-Yy[0]*Yy[0]+Zz[1]*Zz[1]-Zz[0]*Zz[0]-2*par[2]*(Zz[1]-Zz[0])+l[0]*l[0]-l[1]*l[1];

D= 2*(Xx[2]-Xx[3]+par[0]*(Zz[2]-Zz[3]));
E= 2*(Yy[2]-Yy[3]+par[1]*(Zz[2]-Zz[3]));
F=Xx[3]*Xx[3]-Xx[2]*Xx[2]+Yy[3]*Yy[3]-Yy[2]*Yy[2]+Zz[3]*Zz[3]-Zz[2]*Zz[2]-2*par[2]*(Zz[3]-Zz[2])+l[2]*l[2]-l[3]*l[3];

posXYZ[2]=(D*C-F*A)/(E*A-B*D);
posXYZ[1]=-(B*posXYZ[1]+C)/A; 
posXYZ[0]=par[0]*posXYZ[0]+par[1]*posXYZ[1]+par[2];


}

//offset was calculated in rpc_analysis.C
void offset_from_file()
{
	ifstream myfile;
	myfile.open("offset_back.txt", ios::in|ios::out);
	for(int i=0; i<41; i++)
	{
		int tmp;
		myfile>>tmp;
//18 was another ofset to all strips the same, calculated in plot.C
		offset[i]=tmp-18;
	}
	myfile.close();
}
void XYtoXYZ()
{

  offset_from_file();


//Reading cardinal points and fit from file 
int n =9;
ifstream myfile;
myfile.open("pointsfit.txt", ios::in|ios::out);
myfile>>par[0];
myfile>>par[1];
myfile>>par[2];
for(int i=0; i<n; i++)
{
	myfile>>Zz[i];
	myfile>>Xx[i];
	myfile>>Yy[i];
}

myfile.close();
	

/////////////////////////////////////////////////////////////////////


  std::vector<TString> fileList;
//backgroound
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22045050208.root");
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046083458.root");
//little
fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046050208.root");
//
//fallen onece
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046060515.root");
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046064454.root");

//down
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046082416.root");
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046060952.root");
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046065513.root");
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046070926.root");
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046071959.root");


// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046080126.root");
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046081301.root");
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046075139.root");
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046074125.root");
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046073033.root");



//TH2F *Pos_histo = new TH2F("Position","Position",nbin,minData,maxData,41,0,42);
TGraph2D *Pos_histo = new TGraph2D();
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

// for(Int_t t = 0; t< 100000; t++){
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
	Double_t posXYZ[3]={};
	CalculatingPoints(posXYZ, offset[strip-1]+map1->GetPos(), strip*30);
//	cout<<"x,y,z: "<<posXYZ[0]<<", "<<posXYZ[1]<<", "<<posXYZ[2]<<endl;
        Pos_histo->SetPoint(Pos_histo->GetN(), posXYZ[0], posXYZ[1], posXYZ[2]);
     //Pos_histo->Fill(offset[strip-1]+map1->GetPos(),strip);
     }

    }
   }
 }
Pos_histo->SetPoint(Pos_histo->GetN(), 0,0,0);
Pos_histo->SetPoint(Pos_histo->GetN(), -4000,0,0);
TCanvas *C1 =new TCanvas("C1", "C1", 800, 800);
//gStyle-SetOptStat(0);
C1->cd();
Pos_histo->Draw("AP");
}

