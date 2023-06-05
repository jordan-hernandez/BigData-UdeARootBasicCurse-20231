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

int DoubleGaussiansMC(){

    // Valores Límite para las variables
    Double_t Mmin = 1.8; 
    Double_t Mmax = 1.975;

    // Numero de datos
    Int_t datos = 50000;

    // Variables a usar
    RooRealVar M("M", "Mass (K#pi#pi) (GeV)", Mmin, Mmax);

    // ---- MassModel ----

    // -Parámetros Señal-
    RooRealVar mean("mean"," Mass mean",1.875,1.7,1.9,"GeV");

    // Gausiana 1
    RooRealVar width("width"," Mass width",0.010,0.001,0.012,"GeV"); // Sigma1
    RooGaussian Sig("Sig"," Signal PDF",M,mean,width);

    // Gausiana 2
    RooRealVar width2("width2"," Mass width2 ",0.015,0.001,0.05,"GeV"); // Sigma2
    RooGaussian Sig2("Sig2"," Signal PDF B",M,mean,width2);

    // -Parámetros Background-
    RooRealVar c0("c0","c0",0.0,10000.0);
    RooRealVar c1("c1","c1",0.0,1.0);
    RooPolynomial Bkg("Bkg","Exp. Background",M,RooArgList(c0,c1),3);

    // CPesos de Background y señal
    RooRealVar fs("fs","fs",10,0.,1.);
    RooRealVar Ns("Ns","Ns",0.,500);
    RooRealVar Nb("Nb","Nb",1500.,20000);   

    // Suma de las dos gausianas
    RooAddPdf Sumgaus("sumgau","sumgau",RooArgList(Sig,Sig2),RooArgList(fs));

    // Modelo de masa (2 Gausianas + Exponencial)
    RooAddPdf MassModel("MassModel","MassModel",RooArgList(Sumgaus,Bkg),RooArgList(Ns,Nb));

    // ---- Generación de datos ----
    RooDataSet* Data_M = MassModel.generate(M,datos);

    // ---- Fitting ----
    RooFitResult* ResultFit = MassModel.fitTo(*Data_M,Extended(),Minos(kFALSE),Save(kTRUE));

    // ---- Gráficas ----
    Double_t supM = Mmax;
    Double_t infM = Mmin;

    // Número de bines
    Double_t nbin = ((supM-infM)/0.006)+1;

    // Tamaño del Canvas
    int H = 550;
    int W = 650;

    // Creando Canvas
    TCanvas *Canvas1 = new TCanvas("Canvas_MasaK","Canvas_MasaK",500,50,W,H);
    Canvas1->cd();

    // Margenes del Canvas 
    Canvas1->SetLeftMargin(0.11);
    Canvas1->SetRightMargin(0.01);
    Canvas1->SetTopMargin(0.05);
    Canvas1->SetBottomMargin(0.1);

    // Creando frame
    RooPlot* Mframe = M.frame(infM,supM,nbin);
    // Se dibujan los datos y el ajuste
    Data_M->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0));
    MassModel.plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));

    // Histograma Pull
    RooHist* hpullm2 = Mframe->pullHist() ;

    // Dibujar el frame
    MassModel.plotOn(Mframe,Components(Bkg),LineColor(kBlue),LineWidth(2),LineStyle(kDashed),Name("bkg")); 
    Data_M->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("Data"));
    MassModel.plotOn(Mframe);
    Mframe->SetTitle(""); 

    // Ajustes del frame
    Mframe->SetYTitle("Events / 4 MeV"); 
    Mframe->SetLabelSize(0.03,"XY");
    Mframe->SetTitleSize(0.045,"XY");
    Mframe->GetYaxis()->SetRangeUser(0, 2500);   
    Mframe->GetYaxis()->SetNdivisions(505,1);
    Mframe->GetXaxis()->SetNdivisions(1005,1);
    Mframe->GetXaxis()->SetTickLength(0.03);    
    Mframe->GetXaxis()->SetDecimals(1); 
    Mframe->SetTitleOffset(0.85,"X");
    Mframe->SetTitleOffset(1.1,"Y");
    Mframe->Draw();

    // Legend 1
    TLegend *Legend1 = new TLegend(0.18,0.18,0.38,0.3); 
    Legend1->SetTextSize(0.04);
    Legend1->SetFillColor(0);
    Legend1->SetBorderSize(0);
    Legend1->SetFillStyle(0);
    Legend1->AddEntry(Mframe->findObject("Data")," Data","ep"); 
    Legend1->AddEntry(Mframe->findObject("fittotal")," Fit result","l");
    Legend1->AddEntry(Mframe->findObject("bkg"),"Combinatorial background","l");
    Legend1->Draw();

    // Texto adicional
    TLatex *tex2 = new TLatex(0.15,0.88,"CMS");
    tex2->SetNDC();
    tex2->SetTextFont(60);
    tex2->SetTextSize(0.05); 
    tex2->SetLineWidth(2);
    tex2->Draw("L");

    // Legend 2
    TLegend *Legend2 = new TLegend(0.3,0.85,0.85,0.9);
    Legend2->AddEntry("", "5 < p_{T} < 6 Gev, |#eta| < 2.1","");
    Legend2->SetBorderSize(0);
    Legend2->SetTextSize(0.04);
    Legend2->SetFillStyle(0);
    Legend2->SetMargin(0.1);
    Legend2->Draw();

    // 
    Canvas1->Modified();
    gPad->Update();
    gPad->RedrawAxis();
    TLine l;
    l.DrawLine(gPad->GetUxmax(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymin());
    Canvas1->Update();

    // Guardar Gráfica
    Canvas1->Print("../plots/Plot_Mass_K_2pi.png");

    // Canvas del Pull
    TCanvas *Canvas2 = new TCanvas("Pull_Mass_Kpipi", "Pull_Mass_Kpipi",50,50,1200,800 );

    Canvas2->cd();

    // Frame del Pull
    RooPlot* Mframe2 = M.frame(Title(" ")) ;

    // Dibujar Pull
    Mframe2->addPlotable(hpullm2,"P") ;

    // Ajustes del frame
    Mframe2->SetTitle("Pull K#pi#pi M"); 
    Mframe2->SetYTitle(" (Data-Fit)/#sigma");
    Mframe2->GetYaxis()->SetNdivisions(505,1);
    Mframe2->GetXaxis()->SetNdivisions(505,1);
    Mframe2->GetXaxis()->SetTickLength(0.07);   
    Mframe2->SetTitleOffset(0.5,"Y");
    Mframe2->SetTitleSize(0.04,"Y");
    Mframe2->SetLabelSize(0.04,"XY");
    Mframe2->SetTitleSize(0.04,"X");
    Mframe2->Draw();

    // Linea en cero
    TLine *line1 = new TLine(infM,0.0,supM,0.0);
    line1->SetLineColor(4);
    line1->SetLineWidth(1);
    line1->Draw();

    Canvas2->Update();

    // Guardando Gráfica
    Canvas2->Print("../plots/Pull_K2pi.png");

    return 0;

}