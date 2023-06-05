//credit: https://root.cern/doc/master/rf101__basics_8C.html

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
using namespace RooFit;
using namespace std;

int plot_gauss1()
{
    
    // Valores Límite para las variables
    Double_t Xmin = 1.8; 
    Double_t Xmax = 1.975;

    // Numero de datos
    Int_t datos = 1500;

    // Variables a usar
    RooRealVar x("x", "x", Xmin, Xmax);

    // ---- MassModel ----

    // -Parámetros Señal-
    RooRealVar mean("mean", "mean of gaussian", 0.1455, 0.14, 0.16);
    RooRealVar sigma("sigma", "width of gaussian", 0.0012, 0, 0.2);

    // Gaussiana
    RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);

    // -Parámetros Background-
    RooRealVar c("c","c",0);
    RooPolynomial bkg1("bkg1","Background",x,RooArgSet(c),2);

    //Pesos de Background y señal
    RooRealVar Ns("Ns","Ns",0.,500);
    RooRealVar Nb("Nb","Nb",0.,500);

    //Modelo para la masa
    RooAddPdf MassModel("MassModel","MassModel",RooArgList(gauss,bkg1),RooArgList(Ns,Nb));

    // ---- Generación de Datos ----
    RooDataSet *data = MassModel.generate(x, datos);

    // ---- Fitting ----
    RooFitResult* fitres = MassModel.fitTo(*data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));
    //fitres->Print("v");

    // ---- Gráficas ----
    RooPlot *xframe2 = x.frame();

    // Numero de bines
    Int_t bins = 45;

    data->plotOn(xframe2,Name("data"), MarkerSize(1 ),MarkerStyle(8), DrawOption("P"),Binning(bins),DataError(RooAbsData::SumW2),XErrorSize(0));
    MassModel.plotOn(xframe2,LineColor(kBlue),LineWidth(2),Name("fit"));
    //MassModel.plotOn(xframe2,Components(gauss),LineColor(kBlue+1),LineWidth(3),Name("signal")); 
    MassModel.plotOn(xframe2,Components(bkg1),LineColor(kBlue),LineWidth(4), LineStyle(kDashed) ,Name("bkg1")); 

    //Pull
    RooPlot* Mframe = x.frame(0.14,0.16,((0.16-0.14)/0.0008)+1);
    data->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0));
    MassModel.plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));
    RooHist* hpullm2 = Mframe->pullHist() ;

    // Creando Canvas
    TCanvas *c1 = new TCanvas("rf102_dataimport", "rf102_dataimport",50,50,1200,800 );

    // Ajustando el Canvas
    xframe2->GetXaxis()->SetNdivisions(6);
    xframe2->SetMinimum(-1); 
    xframe2->GetYaxis()->SetRangeUser(0, 250);   
    xframe2->SetTitle(" ");
    xframe2->SetYTitle("Events/0.4 Mev"); 
    xframe2->SetXTitle("#Delta M(Gev)");
    xframe2->Draw();

    // Texto adicional
    TLatex *tex2 = new TLatex(0.12,0.836,"CMS");
    tex2->SetNDC();
    tex2->SetTextFont(60);
    tex2->SetTextSize(0.05); 
    tex2->SetLineWidth(2);
    tex2->Draw("L");

    // Legend 1
    auto legend = new TLegend(1,1.5,.62,.38);
    legend->AddEntry("", "29 nb^{-1}(13 Tev)","");
    legend->SetBorderSize(0);
    legend->SetTextSize(0.05);
    legend->SetFillStyle(0);
    legend->SetMargin(0.1);
    legend->Draw();

    // Legend 2
    auto legend2 = new TLegend(0.5,0.8,.4,0.9);
    legend2->AddEntry("", "16 < p_{T} < 24 Gev, |#eta| < 2.1","");
    legend2->SetBorderSize(0);
    legend2->SetTextSize(0.04);
    legend2->SetFillStyle(0);
    legend2->SetMargin(0.1);
    legend2->Draw();

    //Leyenda: objetos de la figura
    TLegend *leg = new TLegend(0.4,0.49,0.83,0.71); 
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(xframe2->findObject("data")," Data","ep"); 
    leg->AddEntry(xframe2->findObject("signal")," Fit ","l");
    leg->AddEntry(xframe2->findObject("bkg1"),"combinatorial background","l");
    leg->Draw();

    c1->Draw();

    // Guardar la figura
    c1->Print("../plots/Plot_DeltaM.png");

    // Canvas para el Pull
    TCanvas *c2 = new TCanvas("rf102_dataimport", "rf102_dataimport",50,50,1200,800 );
    
    // Dibujando el Pull
    framem2->addPlotable(hpullm2,"P") ;

    // Ajustando el Canvas
    RooPlot* framem2 = x.frame(Title(" ")) ;
    framem2->SetYTitle(" (Data-Fit)/#sigma");   
    framem2->GetXaxis()->CenterTitle();
    framem2->GetYaxis()->CenterTitle();
    framem2->SetXTitle("Pull #Delta M"); 
    framem2->GetYaxis()->SetNdivisions(505,1);
    framem2->GetXaxis()->SetNdivisions(505,1);
    framem2->GetXaxis()->SetTickLength(0.07);   
    framem2->SetTitleOffset(0.35,"Y");
    framem2->SetTitleSize(0.04,"Y");
    framem2->SetLabelSize(0.04,"XY");
    framem2->SetTitleSize(0.04,"X");
    framem2->Draw();

    // Linea en cero
    TLine *line1 = new TLine(infM,0.0,supM,0.0);
    line1->SetLineColor(4);
    line1->SetLineWidth(1);
    line1->Draw();

    c2->Draw();
    
    // Guardar la figura
    c2->Print("../plots/Pull_DeltaM.png");

    return 0;
}
