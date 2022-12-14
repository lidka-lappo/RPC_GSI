using namespace std;
TF1 *fitAll[42];
TH1F *hist_sub[42];
TH1F *trip_histo[42];
void plot()
{

  std::vector<TString> fileList;
//background
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22045050208.root");
 //fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046083458.root");

//low
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046050208.root");
//
//felt

//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046060515.root");
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046064454.root");

//measures
//down
 fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046082416.root");
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046060952.root");
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046065513.root");
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046070926.root");
//fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046071959.root");


//up
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046080126.root");
// fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046081301.root");
//  fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046075139.root");
 //fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046074125.root");
  //fileList.push_back("/u/land/mxarepe/unpkd_data/GSI_intern/root_files/r3b_st22046073033.root");

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
//      Pos_histo->Fill(map1->GetPos(),map1->GetChannelId());
      strip_histo[map1->GetChannelId()]->Fill(map1->GetPos());
     }

    }
   }
 }
 //TCanvas *C1 = new TCanvas("C1","C1",600,800);
  int n= 7;
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
 //C1->cd();
// strip_histo[6]->Draw();
 //Pos_histo->Draw("colz");
	int newnBin = 700;
	int newMinHist = 800;
for(int j =0; j<42; j++)
{
	//strip_histo[j]->Scale(1.0/strip_histo[j]->Integral());
   	strip_histo[j]->Smooth(1);
	stringstream strs1, strs2;
	strs1 <<"fit"<< j;
	strs2 <<"fitAll"<< j;
	string tmp11 = strs1.str();
	string tmp22 = strs2.str();
	char *name1 = (char *) tmp11.c_str();
	char *name2 = (char *) tmp22.c_str();
	fit[j] = new TF1(name1, "fitGaussReject", minHist, maxHist, 4);
	fit[j]->SetParameters(4000,strip_histo[j]->GetMaximum(), strip_histo[j]->GetMean(), strip_histo[j]->GetRMS());
	strip_histo[j]->Fit(name1, "RQN");

	fitAll[j] = new TF1(name2, "fitGaussAll", minHist, maxHist, 4);
	fitAll[j]->SetParameters(fit[j]->GetParameters());

	hist_sub[j]= new TH1F(name1, "h1", newnBin, newMinHist, maxHist); //47
	for(int i =0; i<newnBin; i++)//47
	{
//	hist_sub[j]->SetBinContent(i+1, (fitAll[j]->Eval(newMinHist+i)));//+54
	hist_sub[j]->SetBinContent(i+1, (fitAll[j]->Eval(newMinHist+i)-strip_histo[j]->GetBinContent(i+(nbin-newnBin)+1)));//+54
	}
	//h1->Draw();

}

	 		

//fitAll[12]->Draw("SAME");
//hist_sub[12]->Draw("");
for(int i=0; i<n; i++)
{
	for(int j=0; j<6; j++)
	{
		if(i!=6||j!=5){
			C[i]->cd(j+1);
			strip_histo[i*6+j+1]->Draw();
			hist_sub[i*6+j+1]->Draw("SAME");
	 		fitAll[i*6+j+1]->Draw("SAME");
		}
	}
}


//	int ref = hist_sub[11]->GetMaximumBin(); //2
	int ref = 0; //2
	cout<<"Ref: " <<ref+(nbin-newnBin+1)<<endl; //+54
fstream myfile;
myfile.open("offset_official.txt", ios::out);
for(int j =1; j<42; j++)
{
//	int nPeaks = hist_sub[j]->ShowPeaks(2,"", 0.7);
	int shift =ref-(hist_sub[j]->GetMaximumBin());
/*	if(abs(shift)>150)
	{
	cout<<j<<" strip to fix, shit:  "<<shift<<endl;
	//hist_sub[j]->GetXaxis()->SetRange(795,900);
	shift=2000;	
	//hist_sub[j]->GetXaxis()->SetRange(795,hmax);
	}*/
/*	if(j==6)
	{	
		hist_sub[j]->GetXaxis()->SetRange(120,700);//4, 46
		
		shift =ref-(hist_sub[j]->GetMaximumBin());

	}*/
	cout<<j<<"wire: "<<((nbin-newnBin+1)+ref+shift)<<endl; //x15
	shift=shift; //x15
	myfile<<shift<<endl;
} 
myfile.close();
timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  std::cout << std::endl << std::endl;
  std::cout << "Macro finished succesfully." << std::endl;
  std::cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << std::endl << std::endl;

}

