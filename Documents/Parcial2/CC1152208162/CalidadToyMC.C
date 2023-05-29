#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <TLatex.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include <TString.h>
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooGaussModel.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TCut.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "RooHist.h"
#include "RooGenericPdf.h"
#include "RooTruthModel.h"
#include "RooDecay.h"
#include "RooProdPdf.h"
#include "RooEffProd.h"
#include <TRandom3.h>
#include "RooRandom.h"
#include "RooPoisson.h"
#include "RooGamma.h"
#include "RooBernstein.h"
#include "RooGaussian.h"
#include "iostream"
#include "fstream"
#include "RooFitResult.h"
#include "RooCBShape.h"
#include "RooConstVar.h"
#include "TLegend.h"
#include "TPaveStats.h"

using namespace std;
using namespace RooFit;

void CalidadToyMC(Int_t seed = 3, Double_t ptl=12, Double_t pth=70) {
    // Configuración de estilo para la gráfica:
    gStyle -> SetOptTitle(0);

    // Lectura de los datos del Montecarlo:
    TFile f("ToyMontecarlo.root", "READ");
    TTree *tree = (TTree*) f.Get("tree");

    // Declaración de variables para almacenar los resultados del Montecarlo:
    Double_t Nsig, NsigE, Nbkg, NbkgE, Sigma, SigmaE, Mu, MuE, C2, C2E, Nspull, Nsdif, Nbpull, Mupull;
    Int_t status, covQual, badNll;

    // Asignación de las direcciones de memoria de las variables a las ramas del árbol:
    tree -> SetBranchAddress("Nspull", &Nspull);
    tree -> SetBranchAddress("Nsdif", &Nsdif);
    tree -> SetBranchAddress("Nbpull", &Nbpull);
    tree -> SetBranchAddress("Mupull", &Mupull);
    tree -> SetBranchAddress("Nsig", &Nsig);
    tree -> SetBranchAddress("NsigE", &NsigE);
    tree -> SetBranchAddress("Nbkg", &Nbkg);
    tree -> SetBranchAddress("NbkgE", &NbkgE);
    tree -> SetBranchAddress("Mu", &Mu);
    tree -> SetBranchAddress("MuE", &MuE);
    tree -> SetBranchAddress("C2", &C2);
    tree -> SetBranchAddress("C2E", &C2E);
    tree -> SetBranchAddress("Sigma", &Sigma);
    tree -> SetBranchAddress("SigmaE", &SigmaE);
    tree -> SetBranchAddress("status", &status);
    tree -> SetBranchAddress("covQual", &covQual);
    tree -> SetBranchAddress("badNll", &badNll);

    Long64_t entradas = tree -> GetEntries();
    cout << "Entradas: " << entradas << endl;

    // Definición de la variable de masa:
    RooRealVar mu("mu", "Masa", -3.0, 3.0);
    RooDataSet dataMu("data", "data", RooArgSet(mu));

    // Llenado del conjunto de datos con los valores del pull de la masa del Montecarlo:
    for (int evt = 0; evt < tree -> GetEntries(); evt++) {
        tree -> GetEntry(evt);
        mu = Mupull;
        dataMu.add(RooArgSet(mu));
	    }

    // Ajuste a una gaussiana:
    RooRealVar mean("mean", "mean", 0, -1.5, 1.5);
    RooRealVar sigma2("sigma2", "sigma2", 1, 0, 2);
    RooGaussian gauss("gauss", "gauss", mu, mean, sigma2);
    RooFitResult* Fit = gauss.fitTo(dataMu, Extended(), Minos(kFALSE), Save(kTRUE), NumCPU(4)); 
    Fit -> Print("v");

    // Creación del marco para la gráfica:
    RooPlot *frame = mu.frame(Title("Fit de la masa"));
    dataMu.plotOn(frame, MarkerSize(0.8), XErrorSize(0), Name("Data"));
    gauss.plotOn(frame, Components(gauss), LineColor(kMagenta+3), LineWidth(2), Name("Signal"));

    // Creación del lienzo para la figura y la gráfica
    TCanvas *canvas = new TCanvas("canvas", "canvas", 50, 50, 1200, 800);
    frame -> SetYTitle("Events");
    frame -> SetXTitle("B_{c} mass pull");
    frame -> SetTitleSize(0.06, "Y");
    frame -> SetTitleOffset(0.5, "Y");
    frame -> GetYaxis() -> CenterTitle();
    frame -> GetXaxis() -> CenterTitle();
    frame -> GetYaxis() -> SetNdivisions(505, 1);
    frame -> GetXaxis() -> SetNdivisions(505, 1);
    frame -> GetXaxis() -> SetTickLength(0.07);
    frame -> SetLabelSize(0.03, "XY");
    frame -> SetTitleSize(0.05, "Y");
    frame -> SetTitleSize(0.045, "X");
    frame -> SetMinimum(-0.1);
    frame -> Draw();

    // Creación de la leyenda para los parámetros
    TLegend *lpar = new TLegend(0.53, 0.68, 0.8, 0.88);
    lpar -> SetTextSize(0.04);
    lpar -> SetFillColor(0);
    lpar -> SetBorderSize(0);
    lpar -> SetFillStyle(0);
    lpar -> AddEntry("", Form("Mean = %1.3f #pm %1.3f MeV", mean.getVal(), mean.getError()), "");
    lpar -> AddEntry("", Form("Sigma = %1.3f #pm %1.3f MeV", sigma2.getVal(), sigma2.getError()), "");
    lpar -> Draw();

    canvas -> Draw();
    canvas -> Print("pull_montecarlo.png");
	}
