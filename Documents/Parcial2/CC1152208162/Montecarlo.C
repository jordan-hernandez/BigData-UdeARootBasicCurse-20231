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


void Montecarlo(Float_t ntotal=500, Int_t seed=5, Double_t ptl=12, Double_t pth=70)
{
    gStyle -> SetOptTitle(0);
    gStyle -> SetOptFit(0);
    gStyle -> SetOptStat(0);

    // Modelo de masa:
    Double_t supM = 6.5;
    Double_t infM = 6.05;

    RooRealVar M("M", "Mass", 6.05, 6.5);                                                                                                                                          
    RooRealVar mean("mean", "Mean", 6.27, 6.25, 6.3);
    RooRealVar width("width"," Mass width", 0.015, 0.01, 0.04);
    RooGaussian Sig("Sig", " Signal PDF", M, mean, width);

    //Background 
    RooRealVar c("c", "c", -10.0, 10.0);
    RooChebychev bg("bg", "Background", M, c);

    //Pesos de Background y señal
    RooRealVar Ns("Ns", "Ns", 0., 2000);
    RooRealVar Nb("Nb", "Nb", 0., 2000);   

    //Modelo para la masa
    RooAddPdf ModeloMasa("ModeloMasa", "ModeloMasa", RooArgList(Sig, bg), RooArgList(Ns, Nb));
    
    // Proceso de fit:
    RooRandom::randomGenerator() -> SetSeed(seed);

    // No imprimir en pantalla lo que se va haciendo:
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
    Double_t Nsig, NsigE, Nbkg, NbkgE, Sigma, SigmaE, Mu, MuE, C2, C2E, Nspull, Nsdif, Nbpull, Mupull;
    Int_t status, covQual, badNll;

    TFile *file = TFile::Open(Form("ToyMontecarlo.root"), "RECREATE");
    TTree *tree = new TTree("tree", "tree");
    file -> cd();

    tree -> Branch("Nspull", &Nspull);
    tree -> Branch("Nsdif", &Nsdif);
    tree -> Branch("Nbpull", &Nbpull);    
    tree -> Branch("Mupull", &Mupull);
    tree -> Branch("Nsig", &Nsig);
    tree -> Branch("NsigE", &NsigE);
    tree -> Branch("Nbkg", &Nbkg);
    tree -> Branch("NbkgE", &NbkgE);
    tree -> Branch("C2", &C2);
    tree -> Branch("C2E", &C2E);
    tree -> Branch("Sigma", &Sigma);
    tree -> Branch("SigmaE", &SigmaE);
    tree -> Branch("Mu", &Mu);
    tree -> Branch("MuE", &MuE);
    tree -> Branch("status", &status);
    tree -> Branch("covQual", &covQual);
    tree -> Branch("badNll", &badNll);

    // Variables de entrada:
    Double_t Nsi, NsiE, Nbi, NbiE, Ci, CiE, Sigm, SigmE;
    Double_t Mui, MuiE;

    ifstream paramfile (Form("fit_param.txt"));
    if (!paramfile) 
    {
        cout << "Error al leer el archivo de entradas" << endl;
        exit(1);
    }

    paramfile >> Nsi >> NsiE >> Nbi >> NbiE >> Ci >> CiE >> Mui >> MuiE >> Sigm >> SigmE ;
    paramfile.close();
    
    cout << "Ns, Nb, g2 values: " <<  Nsi << " " << Nbi << endl;

    Int_t nruns = 0;
    for(Int_t n=0;n<ntotal;n++) {
        Ns.setVal(Nsi);
        Ns.setError(NsiE);

        Nb.setVal(Nbi);
        Nb.setError(NbiE);

        c.setVal(Ci);
        c.setError(CiE);

        mean.setVal(Mui);
        mean.setError(MuiE);

        width.setVal(Sigm);
        width.setError(SigmE);

        //Generar montecarlo:
        RooDataSet *dataToy = ModeloMasa.generate(RooArgSet(M), Extended(kTRUE));

        //Fit a estos datos:
        RooFitResult* Fit = ModeloMasa.fitTo(*dataToy, Extended(), Minos(kFALSE), Save(kTRUE), NumCPU(6));

        //Se almacenan las variables:
        Nspull = (Ns.getVal()-Nsi)/Ns.getError();
        Nsdif = (Ns.getVal()-Nsi);
        Nbpull = (Nb.getVal()-Nbi)/Nb.getError();   
        Mupull = (mean.getVal()-Mui)/mean.getError();

        Nsig = Ns.getVal();
        NsigE = Ns.getError();

        Nbkg = Nb.getVal();
        NbkgE = Nb.getError();

        Sigm = width.getVal();
        SigmE = width.getError();

        Mu = mean.getVal();
        MuE = mean.getError();

        C2 = c.getVal();
        C2E = c.getError();

        status = Fit -> status();
        covQual = Fit -> covQual();
        badNll = Fit -> numInvalidNLL();

        tree -> Fill();

        nruns++;

        delete dataToy;
        delete Fit;
	    }

    tree -> Write();
	} 
