void offset_from_file()
{
	ifstream myfile;
	myfile.open("measured.txt", ios::in|ios::out);
	for(int i=0; i<41; i++)
	{
		int tmp;
		myfile>>tmp;
//		offset[i]=tmp-50;
//		cout<<offset[i]<<endl;;		
	}
	myfile.close();
}
void plot()
{
	Double_t xRPC[10] ={};	
	Double_t xMes[10] ={};	
	Double_t yRPC[10] ={};
	Double_t yMes[10] ={};
	
	Double_t EMes[10] ={1,1,1,1,1,1,1,1,1,1};
	//Double_t EyRPC[10] ={21,21,21,21,21,21,21,21,21,21};
	Double_t EyRPC[10] ={9,9,9,9,9,9,9,9,9,9};
//	Double_t EyRPC[10] ={17,17,17,17,17,17,17,17,17,17};
	//Double_t ExRPC[10] ={17,17,17,17,17,17,17,17,17,17};
	Double_t ExRPC[10] ={17,17,17,17,17,17,17,17,17,17};
	ifstream myfile;
	ifstream myfile1;
	ifstream myfile2;
	ifstream myfile3;
	myfile.open("measuredX.txt", ios::in|ios::out);
	myfile1.open("RPCX.txt", ios::in|ios::out);
	myfile2.open("measuredY.txt", ios::in|ios::out);
	myfile3.open("RPCY.txt", ios::in|ios::out);
	for(int i=0; i<10; i++)
	{
		myfile>>xMes[i];
		myfile1>>xRPC[i];
		myfile2>>yMes[i];
		myfile3>>yRPC[i];
		//gr->SetPoint(i+1,measure,RPC);
	}

	myfile.close();

  TGraph *gr = new TGraph(10, xMes, xRPC);
  TF1 *fun = new TF1("fun", "[0]*x+[1]", 0, 1500);
  gr->Fit("fun");
  TCanvas *C1 = new TCanvas("C1","C1",800,800);
  gStyle->SetOptStat(0);
  C1->cd();
  gr->SetMarkerSize(1.5);
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);  
  gr->SetTitle("Measured/RPC");
  gr->GetXaxis()->SetTitle("RPC measured position X");
  gr->GetYaxis()->SetTitle("Sticker position X");
  gr->Draw("AP");
  fun->Draw("SAME");
Double_t *offset = fun->GetParameters();
for(int i=0; i<10; i++)
{
xRPC[i]=xRPC[i]-offset[1];
}

  TGraphErrors *grMes = new TGraphErrors(10, xMes, yMes, EMes, EMes);
  TGraphErrors *grRPC = new TGraphErrors(10, xRPC, yRPC, ExRPC, EyRPC);
  TCanvas *C2 = new TCanvas("C2","C2",800,800);
//  gStyle->SetOptStat(0);
  C2->cd();
  grRPC->SetMarkerSize(0.5);
  grMes->SetMarkerSize(0.5);
  grRPC->SetMarkerColor(4);
  grMes->SetMarkerColor(3);
  grRPC->SetMarkerStyle(21);  
  grMes->SetMarkerStyle(20);  
  grRPC->SetTitle("Measured/RPC");
  grRPC->GetXaxis()->SetTitle("X(mm)");
  grRPC->GetYaxis()->SetTitle("Y(mm)");
  grRPC->Draw("AP");
  grMes->Draw("PSAME");
}


