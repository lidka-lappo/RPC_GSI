void CalculatingDistance(double *l, Double_t Xs, Double_t Ys)
{
	Double_t XR[4] ={1396, 780, 179, 180};
	Double_t YR[4] ={1187, 1189, 1146, 696};
	//ifstream myfile;
	//myfile.open("stickers.txt", ios::in|ios::out);
	//Double_t Xs = 147;
	//Double_t Ys = 165;
	for(int i=0; i<4; i++)
	{
	 l[i]=sqrt((Xs-XR[i])*(Xs-XR[i])+(Ys-YR[i])*(Ys-YR[i]));
	//cout<<"l"<<l[i]<<", ";
	}
//	cout<<endl;
	//myfile.close();

}

void CalculatingPoints(Double_t *posXYZ,Double_t *par, Double_t *Xx, Double_t *Yy, Double_t *Zz, Double_t Xs, Double_t Ys)
{

double A, B, C, D, E, F;
double l[4]={1614, 1202, 984, 535};
CalculatingDistance(l, Xs, Ys);
//double l[4]={1614, 984, 1259, 0};
A= 2*(Xx[0]-Xx[1]+par[0]*(Zz[0]-Zz[1]));
B= 2*(Yy[0]-Yy[1]+par[1]*(Zz[0]-Zz[1]));
C=Xx[1]*Xx[1]-Xx[0]*Xx[0]+Yy[1]*Yy[1]-Yy[0]*Yy[0]+Zz[1]*Zz[1]-Zz[0]*Zz[0]-2*par[2]*(Zz[1]-Zz[0])+l[0]*l[0]-l[1]*l[1];

D= 2*(Xx[2]-Xx[3]+par[0]*(Zz[2]-Zz[3]));
E= 2*(Yy[2]-Yy[3]+par[1]*(Zz[2]-Zz[3]));
F=Xx[3]*Xx[3]-Xx[2]*Xx[2]+Yy[3]*Yy[3]-Yy[2]*Yy[2]+Zz[3]*Zz[3]-Zz[2]*Zz[2]-2*par[2]*(Zz[3]-Zz[2])+l[2]*l[2]-l[3]*l[3];
/*A= 2*(Xx[1]-Xx[3]+par[0]*(Zz[1]-Zz[3]));
B= 2*(Yy[1]-Yy[3]+par[1]*(Zz[1]-Zz[3]));
C=Xx[3]*Xx[3]-Xx[1]*Xx[1]+Yy[3]*Yy[3]-Yy[1]*Yy[1]+Zz[3]*Zz[3]-Zz[1]*Zz[1]-2*par[2]*(Zz[3]-Zz[1])+l[0]*l[0]-l[1]*l[1];

D= 2*(Xx[7]-Xx[9]+par[0]*(Zz[7]-Zz[9]));
E= 2*(Yy[7]-Yy[9]+par[1]*(Zz[7]-Zz[9]));
F=Xx[9]*Xx[9]-Xx[7]*Xx[7]+Yy[9]*Yy[9]-Yy[7]*Yy[7]+Zz[9]*Zz[9]-Zz[7]*Zz[7]-2*par[2]*(Zz[9]-Zz[7])+l[2]*l[2]-l[3]*l[3];
*/

posXYZ[1]=(D*C-F*A)/(E*A-B*D);
posXYZ[0]=-(B*posXYZ[1]+C)/A;
posXYZ[2]=par[0]*posXYZ[0]+par[1]*posXYZ[1]+par[2];


}

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
		myfile>>Zz[i];
		cout<<"z"<<Zz[i];
		myfile>>Xx[i];
		cout<<"x"<<Xx[i];
		myfile>>Yy[i];
		cout<<"y"<<Yy[i]<<endl;
		//gr->SetPoint(i+1,measure,RPC);
	}

	myfile.close();

  TF2 *f2 = new TF2("f2","[0]*x+[1]*y+[2]", -1500,10,-500, 700); 
  f2->SetParameters(0,0,-3800);
   TGraph2D *gr = new TGraph2D(n, Xx, Yy, Zz);
	//gr->RemovePoint(0);
  gr->Fit("f2");
  
 //calculating different points
 //

TGraph2D *gr1 = new TGraph2D();
Double_t *par = f2->GetParameters();
Double_t posXYZ[3]={};
Double_t Xs = 147;
Double_t Ys = 165;
ifstream myfile1;
myfile1.open("stickers.txt", ios::in|ios::out);
for(int i=0; i<10; i++)
{

	myfile1>>Xs;
	myfile1>>Ys;
	CalculatingPoints(posXYZ,par, Xx, Yy, Zz, Xs, Ys);
	cout<<"x,y,z : " <<posXYZ[0]<<", "<<posXYZ[1]<<", "<<posXYZ[2]<<endl;
	gr1->SetPoint(gr1->GetN(), posXYZ[2], posXYZ[0], posXYZ[1]);
}
myfile.close();
TGraph2D *gr2 = new TGraph2D(n, Zz, Xx, Yy);
gr2->SetPoint(gr->GetN(), 0, 0, 0);
TCanvas *C1 = new TCanvas("C1","C1",800,800);
  gStyle->SetOptStat(0);
  C1->cd();
  gr1->SetMarkerSize(1.5);
  gr1->SetMarkerColor(3);
  gr1->SetMarkerStyle(20);  
  gr2->SetMarkerSize(1.5);
  gr2->SetMarkerColor(4);
  gr2->SetMarkerStyle(21);  
  gr2->SetTitle("Measured/RPC");
  gr2->GetXaxis()->SetTitle("RPC measured position X");
  gr2->GetYaxis()->SetTitle("Sticker position X");
 // gr->Draw("AP");
  gr2->Draw("AP");
  //f2->Draw("SAME");
  gr1->Draw("SAMEAP");
}


