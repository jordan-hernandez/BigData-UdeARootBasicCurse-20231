/*Aca probaremos nuestra clase Simulator creada*/
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include "Simulator.h"
#include "Simulator.C"
{
    using namespace RooFit;
    
    // ---- Creando Modelo ---
    RooRealVar Mass("Mass", "Mass", 6.05, 6.5);

    RooRealVar Mean("Mean", "Mean", 6.2703 ,6.05, 6.5);
    RooRealVar Sigma("Sigma", "Sigma", .01,0.00001, 1);
    RooGaussian Signal("Signal", "Signal", Mass, Mean, Sigma);

    RooRealVar C("C", "C", -10, 10);
    RooExponential Backg("Backg", "Backg", Mass, C);

    RooRealVar Nsig("Nsig", "Nsig", 0, 13000);
    RooRealVar Nbkg("Nbkg", "Nbkg", 0, 13000);
    RooAddPdf Model("Model", "Model", RooArgList(Signal, Backg), RooArgList(Nsig, Nbkg));
    
    // ---- Inicializando Clase ----
    Simulator MySimulation(Model , Mass , 100);

    // ---- Métodos de la Clase ----

    // Plot Principal
    MySimulation.MainPlot()->Draw();

    // Plot MonteCarlo
    MySimulation.McPlot(RooArgSet(Mean,Sigma,Nsig,Nbkg))->Draw();

}