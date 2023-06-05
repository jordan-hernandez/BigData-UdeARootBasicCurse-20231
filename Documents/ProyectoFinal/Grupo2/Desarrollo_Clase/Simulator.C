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
#include "TObject.h"
#include "Simulator.h"

using namespace RooFit;

int x = 800;
int y = 600;

ClassImp(Simulator) 

/*Descripcion*/

// Constructores
Simulator::Simulator() : TObject(){}
Simulator::Simulator( RooAbsPdf& Model_ , RooRealVar & Obs_ , int const& nbin_): TObject(){
    nbin = nbin_;
    this->Obs = &Obs_;
    this->Model = &Model_;
    this->DataSet = Model->generate(*Obs , 10000); 
    this->FitResult = Model->fitTo(*DataSet , Extended(true) , Save(true));
}

// Función Plot Principal:
// Histograma de Datos y Modelo e Histograma Pull
TCanvas* Simulator::MainPlot(){  
   
   // Dibujando Canvas
    TCanvas *c = new TCanvas("c", "c", x, y);
    c->Divide(1, 2 , 0 ,0);

    c->cd(1); gPad->SetRightMargin(0.01);

    // Primer Frame
    RooPlot* frame = Obs->frame();
    DataSet->plotOn(frame);
    Model->plotOn(frame);
    frame->Draw("A");
    frame->GetYaxis()->CenterTitle(); 
    frame->GetYaxis()->SetTitleSize(0.07); 
    frame->GetYaxis()->SetTitleOffset(0.5);
    frame->GetYaxis()->SetLabelSize(0.045);

    c->cd(2);gPad->SetRightMargin(0.01); gPad->SetBottomMargin(0.3);

    // Segundo Frame (Pull)
    RooHist* pullHist = frame->pullHist();
    RooPlot* pullFrame = Obs->frame();
    pullFrame->addPlotable(pullHist, "P");
    pullFrame->GetYaxis()->SetNdivisions(505);
    pullFrame->GetYaxis()->SetTitle("(Data - Fit) /#sigma");
    pullFrame->GetYaxis()->CenterTitle(); 
    pullFrame->GetYaxis()->SetTitleSize(0.07); 
    pullFrame->GetYaxis()->SetTitleOffset(0.5);
    pullFrame->GetYaxis()->SetLabelSize(0.045);
    pullFrame->GetXaxis()->CenterTitle(); 
    pullFrame->GetXaxis()->SetTitleSize(0.07); 
    pullFrame->GetXaxis()->SetTitleOffset(0.9);
    pullFrame->GetXaxis()->SetLabelSize(0.045);
    pullFrame->Draw("A");
    TLine* zeroLine = new TLine(6.05, 0, 6.5, 0);
    zeroLine->SetLineStyle(2);
    zeroLine->Draw("same");
    c->Update();

    return c;
}

// Plot Testeo con MonteCarlo
TCanvas* Simulator::McPlot(const RooArgSet & Par){ 

    int NumPar = Par.getSize();
    TIterator* ParIter = Par.createIterator();

    // ---- La Trampa: ----
    // Usando RooMCStudy:
    RooMCStudy *MC = new RooMCStudy(*Model,*Obs, Binned(false), Silence(true), Extended(true), FitOptions(Save(true), PrintEvalErrors(0))); 
    //MC->generateAndFit(10000); //---- FIT DE 10 HORAS!!----
    MC->generateAndFit(1000); //---Fit "Rapido"---

    gStyle->SetOptStat(0);
    TCanvas *MC_canvas = new TCanvas("Estudio MC", "Estudio MC", x, y);
    MC_canvas->Divide(NumPar , 2);

    // Iteraciones por cada parámetro
    auto var = ParIter->Next();
    int i = 1;
    while (var) {
      
      RooPlot *ParMeanFrame ;
      RooPlot *ParMeanPullFrame ;

      if (i>NumPar){std::cout<<"Pare en la iteracion "<<i<<endl;break;} ;

      ParMeanFrame = MC->plotParam(*(RooRealVar*)(var), Bins(nbin)); //desreferrenciar el puntero
      ParMeanPullFrame = MC->plotPull(*(RooRealVar*)(var), Bins(nbin), FitGauss(true));

      MC_canvas->cd(i); gPad->SetLeftMargin(0.1);
      ParMeanFrame->GetYaxis()->SetTitleOffset(1.4);
      ParMeanFrame->Draw();
      MC_canvas->Update();

      MC_canvas->cd(i+NumPar); gPad->SetLeftMargin(0.1);
      ParMeanPullFrame->GetYaxis()->SetTitleOffset(1.4);
      ParMeanPullFrame->Draw();
      MC_canvas->Update();

      var = ParIter->Next();
      i += 1;

    } return MC_canvas;}

Simulator::~Simulator(){}


