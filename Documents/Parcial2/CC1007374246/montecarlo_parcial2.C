#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TChain.h>

#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooArgusBG.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "RooVoigtian.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooCBShape.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TLatex.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooNumIntConfig.h"
#include "RooRandom.h"

using namespace RooFit;
using namespace std;


void montecarlo_parcial2(Float_t ntotal=500, Int_t seed=5, Double_t ptl=12, Double_t pth=70)
{
//Este código genera Montecarlo para el modelo de masa
gStyle->SetOptTitle(0);
gStyle->SetOptFit(0);
gStyle->SetOptStat(0);



//---- Mass model -------
Double_t supM = 6.5;
Double_t infM = 6.05;

RooRealVar M("M","Masa",6.05,6.5);                                                                                                                                          
RooRealVar mean("mean", "mean", 6.27, 6.25, 6.3);
RooRealVar width("width"," Mass width",0.015, 0.01, 0.04);
RooGaussian Sig("Sig"," Signal PDF",M,mean,width);

//Background 
RooRealVar c("c","c",-10.0,10.0);
RooChebychev bkg1("bkg1","Background",M,c);



//Pesos de Background y señal
RooRealVar Ns("Ns","Ns",0.,2000);
RooRealVar Nb("Nb","Nb",0.,2000);   

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
Double_t Nsig,NsigE,Nbkg,NbkgE,sigma, sigmaE,mu,muE,c2,c2E,Nspull, Nsdif, Nbpull, Mupull;
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

tree->Branch("c2",&c2);
tree->Branch("c2E",&c2E);

tree->Branch("sigma",&sigma);
tree->Branch("sigmaE",&sigmaE);

tree->Branch("mu",&mu);
tree->Branch("muE",&muE);

tree->Branch("status",&status);
tree->Branch("covQual",&covQual);
tree->Branch("badNll",&badNll);

//****************************************
//    Variables de entreada To_montecarlo.txt
//*****************************************
Double_t nsi, nsie, nbi, nbie, ci, cie,sigm,sigme ;
Double_t mui, muie;

ifstream entrada (Form("To_montecarlo.txt"));
if ( !entrada ) 
{
cout << "No se pudo abrir el archivo entrada fit pr bins" << endl;
exit( 1 );
}

entrada>> nsi >> nsie >> nbi >> nbie >> ci >> cie >> mui >> muie >> sigm >> sigme ;
entrada.close();
cout << "Ns, Nb, g2 values: " <<  nsi << " " << nbi << endl;

Int_t nruns=0;
for(Int_t n=0;n<ntotal;n++)
{
//Ns.setConstant(kTRUE);
Ns.setVal(nsi);
Ns.setError(nsie);

Nb.setVal(nbi);
Nb.setError(nbie);

c.setVal(ci);
c.setError(cie);
  
mean.setVal(mui);
mean.setError(muie);

width.setVal(sigm);
width.setError(sigme);

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

c2 = c.getVal();
c2E = c.getError();

status = fitres->status();
covQual = fitres->covQual();
badNll = fitres->numInvalidNLL();

tree->Fill();

nruns++;

if(nruns%100==0)cout<<"Completado: "<<nruns/ntotal*100<<"%"<<endl;//Porcentaje de avance
delete dataToy;
delete fitres;
}

tree->Write();
}//End analysis

