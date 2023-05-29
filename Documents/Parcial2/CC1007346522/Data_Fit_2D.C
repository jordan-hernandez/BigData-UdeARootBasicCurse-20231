#include "ClassTree1.C"

#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>

#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooArgusBG.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"

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

using namespace std;
using namespace RooFit;

int Data_Fit_2D(){

    // Creamos el TTree
    TChain * Chain1 = new TChain("Tree1","");
    Chain1->Add("RootFile1.root");
    TTree *mytree = (TTree*) Chain1;

    // Usando ClassTree1 desempaquetamos el Tree
    ClassTree1 Tree4(mytree);

    cout << "* This tree has " << Tree4.fChain->GetEntries() << " entries.\n\n";

    // Valores Límete para las variables
    Double_t Mmin = 6.05; 
    Double_t Mmax = 6.5;
    Double_t Tmin = 0.3; 
    Double_t Tmax = 2.6;
    Double_t T_errmin = 0.0001; 
    Double_t T_errmax = 0.5; 

    // Variables a usar
    RooRealVar M("M", "Mass B_{c} (GeV)", Mmin, Mmax);
    RooRealVar T("T", "#tau", Tmin, Tmax);
    RooRealVar T_err("T_err", "Tau error", T_errmin,T_errmax);

    // Dataset con valores de M y T
    RooDataSet Data("Data_mass_tau", "Data Masa y Tau", RooArgSet(M,T));

    // Llenando el RooDataSet
    Long64_t nentries = Tree4.fChain->GetEntries();
    Int_t nTen = nentries/10;
    Long64_t nbytes = 0, nb = 0;

    for (int jentry=0; jentry<nentries; jentry++) 
    {
        Long64_t ientry = Tree4.LoadTree(jentry);
        if (ientry < 0) break;
        if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
        if(jentry==nentries-1) cout<<endl;

        nb = Tree4.fChain->GetEntry(ientry);   nbytes += nb;
    
        // No permitir que la razon entre Tau y su error sea menor a 5
        if ((Tree4.Tau/Tree4.Tau_err)<5) continue;

        // Verificando que los valores de M, T y T_error esten dentro de los límites
        if ((Tree4.Tau<Tmin)||(Tree4.Tau>Tmax)) continue;
        if ((Tree4.Tau_err<T_errmin)||(Tree4.Tau_err>T_errmax)) continue;
        if ((Tree4.mass<Mmin)||(Tree4.mass>Mmax)) continue;

        M = Tree4.mass;
        T = Tree4.Tau;

        Data.add(RooArgSet(M,T)); // Añadir al Dataset
    }

    //---- MassModel ----

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

    // ---- Tau Model ----
    RooRealVar alpha("alpha", "alpha", -3.0, 3.0);
    RooExponential TauModel("TauModel", "TauModel", T, alpha);

    // ---- Total Model ----
    RooProdPdf TotalModel("TotalModel", "TotalModel", MassModel, TauModel);   

    // ---- Fitting ----
    RooFitResult* ResultFit = TotalModel.fitTo(Data,Extended(),Minos(kFALSE),Save(kTRUE),ConditionalObservables(T));

    // Construyendo Histograma 2D
    Int_t Nbins = 20;
    TH1 *HistModel = TotalModel.createHistogram("TotalModel", M, Binning(Nbins), YVar(T, Binning(Nbins)));
    HistModel->SetLineColor(38);
    
    // Creando Frame de Masa
    RooPlot *frame = M.frame(Title("Mass"));

    // Graficando Datos y Fiteo de la Masa
    Data.plotOn(frame);
    TotalModel.plotOn(frame, ProjWData(T,Data),LineColor(38),XErrorSize(0)); 

    // Creando Frame de Tiempo de Vida
    RooPlot *frame2 = T.frame(Title("Lifetime"));

    // Graficando Datos y Fiteo de la Masa
    Data.plotOn(frame2);
    TotalModel.plotOn(frame2, ProjWData(M,Data),LineColor(38),XErrorSize(0)); 
    
    // Creando Canvas
    TCanvas *Canvas1 = new TCanvas("Mass-Lifetime Plot", "Mass-Lifetime Plot", 1200, 450);
    Canvas1->Divide(3); // Creando 3 Pads

    // Pad 1
    Canvas1->cd(1);
    gPad->SetLeftMargin(0.20);
    HistModel->SetTitle("PDF Mass-Lifetime");
    HistModel->GetZaxis()->SetTitleOffset(2.5);
    HistModel->GetXaxis()->SetDecimals(1); 
    HistModel->GetYaxis()->SetDecimals(1);
    HistModel->GetZaxis()->SetDecimals(1);
    HistModel->SetTitleOffset(1.8,"X");
    HistModel->SetTitleOffset(1.6,"Y");
    HistModel->Draw("surf");
    Data.Draw("same");
    

    // Pad 2
    Canvas1->cd(2);
    gPad->SetLeftMargin(0.15);
    frame->GetYaxis()->SetTitleOffset(1.6);
    frame->SetTitleOffset(1.0,"X");
    frame->SetTitleOffset(1.2,"Y");
    frame->Draw();

    // Pad 3
    Canvas1->cd(3);
    gPad->SetLeftMargin(0.15);
    frame2->GetYaxis()->SetTitleOffset(1.6);
    frame2->SetTitleOffset(1.0,"X");
    frame2->SetTitleOffset(1.2,"Y");
    frame2->Draw();
    
    // Guardando Plot
    Canvas1->Print("plots/MasayTau_Fiteo.png");
    return 0;
}