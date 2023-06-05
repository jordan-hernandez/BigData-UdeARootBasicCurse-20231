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
#include "RooGaussModel.h"
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



void analisis_MC(Int_t seed = 3, Double_t ptl=12, Double_t pth=70)
{
//Este código realiza una grafica al pull del montecarlo para ver si este es correcto. Notamos la Gaussiana sobre la imagen final.
 gStyle->SetOptTitle(0);

//Lectura de los datos del Montecarlo
TFile f ("Datos_MC.root", "READ");
TTree *tree = (TTree*) f.Get ("tree");

//Variables
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
   Double_t        c2;
   Double_t        c2E;
   Double_t        sigma;
   Double_t        sigmaE;
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
tree->SetBranchAddress("sigma", &sigma);
tree->SetBranchAddress("sigmaE", &sigmaE);
tree->SetBranchAddress("status", &status);
tree->SetBranchAddress("covQual", &covQual);
tree->SetBranchAddress("badNll", &badNll);

Long64_t nentries = tree->GetEntries();
cout<<" Entries : "<<nentries<<endl;

RooRealVar Mu("Mu","Masa",-5.0,5.0);
RooDataSet dataMu("data","data",RooArgSet(Mu));

for (int evt = 0; evt <tree->GetEntries();evt++)
{
tree->GetEntry(evt);
Mu =  Mupull;
dataMu.add(RooArgSet(Mu));
}

//Modelo de la masa
RooRealVar mean("mean", "mean", 0, -2.0, 2.0);
RooRealVar sigma2("sigma2", "sigma2", 1, 0, 2);
RooGaussian gauss("gauss", "gauss", Mu, mean, sigma2);
RooFitResult* fitres = gauss.fitTo(dataMu,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));//Ajuste a Gauss
fitres->Print("v");

RooPlot *frame = Mu.frame(Title("Fit Mass"));
dataMu.plotOn(frame,MarkerSize(0.8),XErrorSize(0),Name("Data"));
gauss.plotOn(frame,Components(gauss),LineColor(kRed),LineWidth(2),Name("signal")); 

//Lienzo para la figura y gráfica
TCanvas *c1 = new TCanvas("rf102_dataimport", "rf102_dataimport",50,50,1200,800 );
frame->SetYTitle("Events"); 
frame->SetXTitle("#DeltaM Pull"); 
frame->SetTitleSize(0.06,"Y");
frame->SetTitleOffset(0.5,"Y");
frame->GetYaxis()->CenterTitle();   
frame->GetXaxis()->CenterTitle();
frame->GetYaxis()->SetNdivisions(505,1);
frame->GetXaxis()->SetNdivisions(505,1);
frame->GetXaxis()->SetTickLength(0.07); 
frame->SetLabelSize(0.03,"XY");
frame->SetTitleSize(0.05,"Y");
frame->SetTitleSize(0.045,"X");
frame->SetMinimum(-0.1); 
frame->Draw();

TLegend *legpar = new TLegend(0.53,0.68,0.8,0.88);
legpar->SetTextSize(0.04);
legpar->SetFillColor(0);
legpar->SetBorderSize(0);
legpar->SetFillStyle(0);
legpar->AddEntry("",Form("Mean = %1.3f #pm %1.3f MeV",mean.getVal(), mean.getError()),"");
legpar->AddEntry("",Form("sigma = %1.3f #pm %1.3f MeV",sigma2.getVal(), sigma2.getError()),"");
legpar->Draw();

c1->Draw();
c1->Print("../plots/Montecarlo_pull1.png");
}