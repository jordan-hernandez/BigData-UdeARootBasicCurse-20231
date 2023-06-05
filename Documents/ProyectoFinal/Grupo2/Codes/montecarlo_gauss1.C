using namespace RooFit;
using namespace std;

void montecarlo_gauss1(Float_t ntotal=2000, Int_t seed=5, Double_t ptl=12, Double_t pth=70)
{
//Modelo


RooRealVar M("M","Masa",0.14, 0.16);                                                                                                                                          
RooRealVar mean("mean", "mean", 0.1455, 0.14, 0.16);
RooRealVar width("width"," Mass width",0.001255, 0, 0.2);
RooGaussian Sig("Sig"," Signal PDF",M,mean,width);

//Background 
RooRealVar c("c","c",0,-1,1);
RooPolynomial bkg1("bkg1","Background",M,RooArgSet(c),0);



//Pesos de Background y señal
RooRealVar Ns("Ns","Ns",0.,1500);
RooRealVar Nb("Nb","Nb",0.,1500); 

//Modelo para la masa
RooAddPdf MassModel("MassModel","MassModel",RooArgList(Sig,bkg1),RooArgList(Ns,Nb));


//------------ Fit procedure -------------------
RooRandom::randomGenerator()->SetSeed(seed);

// esto que sigue solo es para que no imprima en pantalla todo lo que hace del ajuste
RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
RooMsgService::instance().setSilentMode(kTRUE);
RooMsgService::instance().setStreamStatus(1,false);
RooMsgService::instance().getStream(1).removeTopic(Integration) ;  
RooMsgService::instance().getStream(1).removeTopic(Minimization) ;  
RooMsgService::instance().getStream(1).removeTopic(Fitting) ;  
RooMsgService::instance().getStream(1).removeTopic(NumIntegration) ;
RooMsgService::instance().getStream(1).removeTopic(Optimization) ;
RooMsgService::instance().getStream(1).removeTopic(ObjectHandling) ;
RooMsgService::instance().getStream(1).removeTopic(Eval) ;
RooMsgService::instance().Print() ;

//Variables
Double_t Nsig,NsigE,Nbkg,NbkgE,sigma, sigmaE,mu,muE,Nspull, Nsdif, Nbpull, Mupull;
Int_t status,covQual,badNll;

TFile *file = TFile::Open(Form("Datos_MC.root"),"RECREATE");
TTree *tree = new TTree("tree","tree");
file->cd();

tree->Branch("Nspull",&Nspull);
tree->Branch("Nsdif",&Nsdif);
tree->Branch("Nbpull",&Nbpull);    
tree->Branch("Mupull",&Mupull);

tree->Branch("Nsig",&Nsig);
tree->Branch("NsigE",&NsigE);

tree->Branch("Nbkg",&Nbkg);
tree->Branch("NbkgE",&NbkgE);

tree->Branch("sigma",&sigma);
tree->Branch("sigmaE",&sigmaE);

tree->Branch("mu",&mu);
tree->Branch("muE",&muE);

tree->Branch("status",&status);
tree->Branch("covQual",&covQual);
tree->Branch("badNll",&badNll);

Double_t nsi, nsie, nbi, nbie, ci, cie,sigm,sigme ;
Double_t mui, muie;

Int_t nruns=0;
for(Int_t n=0;n<ntotal;n++)
{
// cout<<"entro"<<endl;
mean.setVal(0.145);
mean.setError(0.000055);

width.setVal(0.0012);
c.setVal(0.0001);
Ns.setVal(1500);
Ns.setError(0.995);
Nb.setVal(500);
mui = 0.145;
nsi = 500;


//Generar montecarlo
RooDataSet *dataToy = MassModel.generate(RooArgSet(M), Extended(kTRUE));

//Fit a estos datos
RooFitResult* fitres = MassModel.fitTo(*dataToy,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(6));

//Se almacenan las variables
Nspull = (Ns.getVal()-nsi)/Ns.getError();
Nsdif = (Ns.getVal()-nsi);
Nbpull = (Nb.getVal()-nbi)/Nb.getError();   
Mupull = (mean.getVal()-mui)/mean.getError();

Nsig = Ns.getVal();
NsigE = Ns.getError();

Nbkg = Nb.getVal();
NbkgE = Nb.getError();

sigm = width.getVal();
sigme = width.getError();

mu = mean.getVal();
muE = mean.getError();


status = fitres->status();
covQual = fitres->covQual();
badNll = fitres->numInvalidNLL();

tree->Fill();

nruns++;

if(nruns%10==0)cout<<"Completado: "<<nruns/ntotal*100<<"%"<<endl;//Porcentaje de avance
delete dataToy;
delete fitres;
}
tree->Write();
}