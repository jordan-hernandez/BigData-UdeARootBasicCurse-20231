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
void sinteticos2_MC(Int_t ntotal=1000, Int_t seed=5)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //gStyle->SetErrorX(0);

  Double_t Mmin = 6.05; 
  Double_t Mmax = 6.5;

  RooRealVar Bcm("bcm","bcm",6.05,6.5);
  
  //---- Mass model -------
  Double_t supM = Mmax;
  Double_t infM = Mmin;
  
  RooRealVar a0("a0","a0",-10,10) ;
  RooChebychev bkg("bkg","Background",Bcm,a0) ;
  
  //**** define nominal signal ****
  //gaussians                                                                                                     
  RooRealVar mean("mean","mean of gaussians",6.275,6.05,6.5) ;
  RooRealVar sigma1("sigma1","sigma1 of gaussians",0.01,0.001,0.015) ;
  RooGaussian sig1("sig1","Signal component 1",Bcm,mean,sigma1) ;
    
  //********final PDF ********
  RooRealVar nsig("nsig","number of signal events",0.,1500) ;
  RooRealVar nbkg("nbkg","number of background events",0,1500) ;
  RooAddPdf  model("model","sig+bkg", RooArgList(sig1,bkg), RooArgList(nsig,nbkg)) ;

  //------------ Fit procedure -------------------
  
  //RooRandom::randomGenerator()->SetSeed(3);
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
  
  Double_t llh,edm,Nsig,NsigE,Nbkg,NbkgE,g1,g1E,mu,muE,conscheb,conschebE, Nspull, Nsdif, Nbpull, Mupull;
  Int_t status,covQual,badNll;
  TFile *file = new TFile("ejemplo.root","RECREATE");
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

  tree->Branch("mu",&mu);
  tree->Branch("muE",&muE);

  tree->Branch("g1",&g1);
  tree->Branch("g1E",&g1E);
  
  tree->Branch("conscheb",&conscheb);
  tree->Branch("conschebE",&conschebE);

  tree->Branch("edm",&edm);
  tree->Branch("llh",&llh);
  tree->Branch("status",&status);
  tree->Branch("covQual",&covQual);
  tree->Branch("badNll",&badNll);

  
  //****************************************
  //    Input values from DATA fit in pt bins
  //*****************************************
  Double_t nisgf, nsigef, nbkgf, nbkgef, a0f, a0ef;
  Double_t meantxt, meantxte, g1i, g1ie;
 
  ifstream entrada ("fit.txt");
  entrada.is_open();
  if ( !entrada ) 
    {
      cout << "No se pudo abrir el archivo entrada 'resultTFit'" << endl;
      exit( 1 );
    }
  entrada>> nisgf >> nsigef >> nbkgf >> nbkgef >> a0f >> a0ef >> meantxt >> meantxte >> g1i >> g1ie ;
  entrada.close();
  cout << "nsig, nbkg values: " <<  nisgf << " " << nbkgf << endl;
  //return;
  
  Int_t nruns=0;
  
  for(Int_t n=0;n<ntotal;n++)
  {
    //nsig.setConstant(kTRUE);
    nsig.setVal(nisgf);
    nsig.setError(nsigef);
    
    nbkg.setVal(nbkgf);
    nbkg.setError(nbkgef);

    a0.setVal(a0f);
    a0.setError(a0ef);
                
    mean.setVal(meantxt);
    mean.setError(meantxte);
    
    sigma1.setVal(g1i);
    sigma1.setError(g1ie);
    
    RooDataSet *dataToy = model.generate(RooArgSet(Bcm), Extended(kTRUE));
    
    RooFitResult* fitres = model.fitTo(*dataToy,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));

    Nspull = (nsig.getVal()-nisgf)/nsig.getError();
    Nsdif = (nsig.getVal()-nisgf);
    Nbpull = (nbkg.getVal()-nbkgf)/nbkg.getError();   
    Mupull = (mean.getVal()-meantxt)/mean.getError();

    Nsig = nsig.getVal();
    NsigE = nsig.getError();
    
    Nbkg = nbkg.getVal();
    NbkgE = nbkg.getError();

    g1 = sigma1.getVal();
    g1E = sigma1.getError();

    mu = mean.getVal();
    muE = mean.getError();

    conscheb = a0.getVal();
    conschebE = a0.getError();
    
    llh = fitres->minNll();
    status = fitres->status();
    covQual = fitres->covQual();
    edm = fitres->edm();
    badNll = fitres->numInvalidNLL();

    tree->Fill();

    nruns++;
   /* 
    TCanvas* canv_nominal = CreateCanvasNom("canv_nominal", fitres, dataToy, Bcm, supM, infM, model, sig1, bkg, nsig, nbkg, sigma1, mean,ptl,pth); 
    canv_nominal->Print(Form("mass_BsFit_ptbins_ToyMC_%1.0f_%1.0f.png",ptl,pth));
    canv_nominal->Print(Form("mass_BsFit_ptbins_ToyMC_%1.0f_%1.0f.pdf",ptl,pth));
    */
    
    delete dataToy;
    delete fitres;

  }

 tree->Write();
 file->Write("",TObject::kOverwrite);
 
 cout<<" for end: "<<nruns<<endl;

}//End analysis
