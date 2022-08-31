void VisPoints()
{
	int n =9;
	Double_t Xx[9] ={};
	Double_t Yy[9] ={};
	Double_t Zz[9] ={};
	
	ifstream myfile;
	myfile.open("points.txt", ios::in|ios::out);
	for(int i=0; i<n; i++)
	{
		myfile>>Xx[i];
		cout<<"x"<<Xx[i];
		myfile>>Yy[i];
		cout<<"y"<<Yy[i];
		myfile>>Zz[i];
		cout<<"z"<<Zz[i]<<endl;
		//gr->SetPoint(i+1,measure,RPC);
	}

	myfile.close();

  TF2 *f2 = new TF2("f2","[0]*x+[1]*y+[2]", -1300,-10,-400, 700); 
  f2->SetParameters(0,0,-3800);
   TGraph2D *gr = new TGraph2D(n, Yy, Zz, Xx);
  gr->Fit("f2");
  TCanvas *C1 = new TCanvas("C1","C1",800,800);
  gStyle->SetOptStat(0);
  C1->cd();
  gr->SetMarkerSize(1.5);
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);  
  gr->SetTitle("Measured/RPC");
  gr->GetXaxis()->SetTitle("RPC measured position X");
  gr->GetYaxis()->SetTitle("Sticker position X");
 // gr->Draw("AP");
  f2->Draw("surf1");
  gr->Draw("SAMEAP");
}


