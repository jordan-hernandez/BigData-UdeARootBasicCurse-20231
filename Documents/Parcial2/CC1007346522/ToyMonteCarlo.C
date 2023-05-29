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


void ToyMonteCarlo(Int_t ntotal=1,Int_t seed = 17)
{
  // Ajustes de estilo
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  // Valores límtes para la Masa
  Double_t Mmin = 6.05; 
  Double_t Mmax = 6.5;

  // Variable a usar
  RooRealVar M("M"," M(B_{c}) (GeV)",Mmin,Mmax);
  
  // ---- MassModel ----

  // -Parámetros Señal-
  RooRealVar mean("mean"," Mass mean",6.27,6.25,6.3,"GeV");

  // Gausiana 1
  RooRealVar width("width"," Mass width",0.10,0.001,0.3,"GeV"); // Sigma1
  RooGaussian Sig("Sig"," Signal PDF",M,mean,width);

  // Gausiana 2
  RooRealVar width2("width2"," Mass width2 ",0.20,0.001,10.0,"GeV"); // Sigma2
  RooGaussian Sig2("Sig2"," Signal PDF B",M,mean,width2);

  // -Parámetros Background-
  RooRealVar c("c","c",-10.0,10.0);
  RooExponential Bkg("Bkg","Exp. Background",M,c);

  // Cantidad de datos por cada componente
  RooRealVar Ns("Ns","Ns",0.,2000);
  RooRealVar Nb("Nb","Nb",0.,2000);   
  RooRealVar fs("fs","fs",0.8,0.,1.);

  // Suma de las dos gausianas
  RooAddPdf Sumgaus("sumgau","sumgau",RooArgList(Sig,Sig2),RooArgList(fs));

  // Modelo de masa (2 Gausianas + Exponencial)
  RooAddPdf MassModel("MassModel","MassModel",RooArgList(Sumgaus,Bkg),RooArgList(Ns,Nb));

  // Hacemos que la semilla de generación sea la misma para todos los procesos
  RooRandom::randomGenerator()->SetSeed(seed);

  // Para que la consola no muestre tantos datos no útiles
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
  
  // Definiendo variables a utilizar 
  Double_t llh,edm,Nsig,NsigE,Nbkg,NbkgE,g1,g1E,mu,muE,g2,g2E,exp,expE, FS, FSE, Nspull, Nsdif, Nbpull, Mupull;
  Int_t status,covQual,badNll;

  // Creamos un archivo .root para guardar los datos generados
  TFile *File2 = TFile::Open("ToyMC.root","RECREATE");
  TTree *Tree2 = new TTree("Tree2","Tree2");
  File2->cd();

  // Definimos las Ramas del archivo
  Tree2->Branch("Nspull",&Nspull);
  Tree2->Branch("Nsdif",&Nsdif);
  Tree2->Branch("Nbpull",&Nbpull);    
  Tree2->Branch("Mupull",&Mupull);
  
  Tree2->Branch("Nsig",&Nsig);
  Tree2->Branch("NsigE",&NsigE);
  
  Tree2->Branch("Nbkg",&Nbkg);
  Tree2->Branch("NbkgE",&NbkgE);

  Tree2->Branch("mu",&mu);
  Tree2->Branch("muE",&muE);

  Tree2->Branch("g1",&g1);
  Tree2->Branch("g1E",&g1E);

  Tree2->Branch("g2",&g2);
  Tree2->Branch("g2E",&g2E);

  Tree2->Branch("FS",&FS);
  Tree2->Branch("FSE",&FSE);
  
  Tree2->Branch("exp",&exp);
  Tree2->Branch("expE",&expE);

  Tree2->Branch("edm",&edm);
  Tree2->Branch("llh",&llh);
  Tree2->Branch("status",&status);
  Tree2->Branch("covQual",&covQual);
  Tree2->Branch("badNll",&badNll);

  // Definimos variables temporales para realizar el llenado del archivo
  Double_t nsi, nsie, nbi, nbie, ci, cie, fsi, fsie;
  Double_t mui, muie, g1i, g1ie, g2i, g2ie;
  
  // Leemos los parámetros que encontramos en el Fit de los datos
  ifstream entrada ("TotalFit.txt");
  //entrada.is_open();
  if ( !entrada ) 
    {
      cout << "No se pudo abrir el archivo entrada fit pr bins" << endl;
      exit( 1 );
    }
    // Sacando los valores
  entrada>> nsi >> nsie >> nbi >> nbie >> ci >> cie >> fsi >> fsie >> mui >> muie >> g1i >> g1ie >> g2i >> g2ie ;
  entrada.close();
  cout << "Ns, Nb, g2 values: " <<  nsi << " " << nbi << "  " << g2i << endl;
  
  Int_t nruns=0;
  
  for(Int_t n=0;n<ntotal;n++)
  {
    // Definimos los valores de los parámetros y su error
    Ns.setVal(nsi);
    Ns.setError(nsie);
    
    Nb.setVal(nbi);
    Nb.setError(nbie);

    c.setVal(ci);
    c.setError(cie);
        
    fs.setVal(fsi);
    fs.setError(fsie);
        
    mean.setVal(mui);
    mean.setError(muie);
    
    width.setVal(g1i);
    width.setError(g1ie);
    
    width2.setVal(g2i);
    width2.setError(g2ie);

    // Generamos datos aleatorios
    RooDataSet *dataToy = MassModel.generate(RooArgSet(M), Extended(kTRUE));
    
    // Fiteamos los datos con el modelo encomtrado
    RooFitResult* fitres = MassModel.fitTo(*dataToy,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));

    // Diferencia entre el último fiteo y el primero
    Nspull = (Ns.getVal()-nsi)/Ns.getError();
    Nsdif = (Ns.getVal()-nsi);
    Nbpull = (Nb.getVal()-nbi)/Nb.getError();   
    Mupull = (mean.getVal()-mui)/mean.getError();

    // Llenamos las variables con los valores obtenidos y sus errores
    Nsig = Ns.getVal();
    NsigE = Ns.getError();
    
    Nbkg = Nb.getVal();
    NbkgE = Nb.getError();

    g1 = width.getVal();
    g1E = width.getError();

    mu = mean.getVal();
    muE = mean.getError();

    FS = fs.getVal();
    FSE = fs.getError();
    
    g2 = width2.getVal();
    g2E = width2.getError();

    exp = c.getVal();
    expE = c.getError();
    
    llh = fitres->minNll();
    status = fitres->status();
    covQual = fitres->covQual();
    edm = fitres->edm();
    badNll = fitres->numInvalidNLL();

    // Llenamos el archivo
    Tree2->Fill();

    nruns++;

    // Eliminamos los datos generados y el resultado del fiteo
    delete dataToy;
    delete fitres;

    }

    // Escribimos el archivo
    Tree2->Write();
    
    cout<<"Total: "<<nruns<<endl; // Mensaje de finalizado

}
