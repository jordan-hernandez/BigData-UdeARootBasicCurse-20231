
using namespace std;
using namespace RooFit;

void analisis()
{

 gStyle->SetOptTitle(0);
 //setTDRStyle();

TFile f ("ToyMC_Bc_bins.root", "READ");
TTree *tree = (TTree*) f.Get ("tree");

   Double_t        Nspull;
   Double_t        Nsdif;
   Double_t        Nbpull;
   Double_t        Mupull;
   Double_t        Nsig;
   Double_t        NsigE;
   Double_t        Nbkg;
   Double_t        NbkgE;
   Double_t        mu;
   Double_t        muE;
   Double_t        a;
   Double_t        aE;
   Double_t        sigma;
   Double_t        sigmaE;

   Double_t        llh;
   Double_t        edm;
   Int_t           status;
   Int_t           covQual;
   Int_t           badNll;


tree->SetBranchAddress("Nspull", &Nspull);
tree->SetBranchAddress("Nsdif", &Nsdif);
tree->SetBranchAddress("Nbpull", &Nbpull);
tree->SetBranchAddress("Mupull", &Mupull);
tree->SetBranchAddress("Nsig", &Nsig);
tree->SetBranchAddress("NsigE", &NsigE);
tree->SetBranchAddress("Nbkg", &Nbkg);
tree->SetBranchAddress("NbkgE", &NbkgE);
tree->SetBranchAddress("mu", &mu);
tree->SetBranchAddress("muE", &muE);
tree->SetBranchAddress("a", &a);
tree->SetBranchAddress("aE", &aE);
tree->SetBranchAddress("sigma", &sigma);
tree->SetBranchAddress("sigmaE", &sigmaE);
tree->SetBranchAddress("status", &status);
tree->SetBranchAddress("covQual", &covQual);
tree->SetBranchAddress("badNll", &badNll);
tree->SetBranchAddress("llh", &llh);
tree->SetBranchAddress("edm", &edm);


Long64_t nentries = tree->GetEntries();
cout<<" Entries : "<<nentries<<endl;

Double_t MminP = -3.0;
Double_t MmaxP = 3.0;  

RooRealVar Mu("Mu","Masa",MminP,MmaxP);
RooDataSet dataMu("data","data",RooArgSet(Mu));

RooRealVar Yi("Yi"," B_{c} yield Pull",0, MminP,MmaxP); 
RooDataSet dataYi("dataYi","dataYi",RooArgSet(Yi));

RooRealVar Yib("Yib"," B_{c} yield bkg Pull",0, MminP,MmaxP); 
RooDataSet dataYib("dataYib","dataYib",RooArgSet(Yib)); 


for (int evt=0; evt < nentries; evt++) 
{
  tree->GetEntry(evt);

  Mu= Mupull;
  Yi= Nspull;
  Yib= Nbpull; 


  dataMu.add(RooArgSet(Mu));
  dataYi.add(RooArgSet(Yi));
  dataYib.add(RooArgSet(Yib));

}

cout<<" Entries : "<<nentries<<endl;



RooRealVar meanYi("meanYi","meanYi",0.0,MminP,MmaxP);
RooRealVar widthYi("widthYi","widthYi",1.0,0.0,5.0);
RooGaussian SigYi("SigYi","SignalYi",Yi,meanYi,widthYi); 
 

RooFitResult* fitYi = SigYi.fitTo(dataYi,Minos(kFALSE),Save(kTRUE), NumCPU(4));
fitYi->Print("v");


RooRealVar meanYib("meanYib","meanYib",0.0,MminP,MmaxP);
RooRealVar widthYib("widthYib","widthYib",1.0,0.0,5.0);
RooGaussian SigYib("SigYib","SignalYib",Yib,meanYib,widthYib); 
 

 
RooFitResult* fitYib = SigYib.fitTo(dataYib,Minos(kFALSE),Save(kTRUE), NumCPU(4));
fitYib->Print("v");



RooRealVar meanMu("meanMu","meanMu",0.0,MminP,MmaxP);
RooRealVar widthMu("widthMu","widthMu",1.0,0.0,5.0);
RooGaussian SigMu("SigMu","SignalMu",Mu,meanMu,widthMu); 



RooFitResult* fitMu = SigMu.fitTo(dataMu,Minos(kFALSE),Save(kTRUE), NumCPU(4));
fitMu->Print("v");


// RooPlot *frame = Yi.frame(Title("Fit Mass"));
// dataYi.plotOn(frame);
// SigYi.plotOn(frame,Components(SigYi),LineColor(kRed),LineWidth(2),Name("signal")); 

// TCanvas *c1 = new TCanvas("rf102_dataimport", "rf102_dataimport",50,50,1200,800 );
// frame->Draw();

// c1->Draw();
// c1->Print("Yipull.png");
}