#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TMathBase.h"
using namespace RooFit;
using namespace TMath;



int fit_parcial2() {
//Este código lee Datos_organizados.root y hace un fit a la variable masa usando un modelo
//gaussiano para la señal y un polinomio de Chevychev para el Background
Float_t mass , tau,taue;

TFile f ("Datos_organizados.root", "READ");
TTree *tree = (TTree*) f.Get ("mytree");
tree->SetBranchAddress("mass",&mass);
tree->SetBranchAddress("tau", &tau);
tree->SetBranchAddress("taue", &taue);

RooRealVar M("M","Masa",6.05,6.5);
RooDataSet data("data","data",RooArgSet(M));
Int_t nentries = tree->GetEntries();
cout<<"Entradas: "<<nentries<<endl;

for (int evt = 0; evt <tree->GetEntries();evt++)
{
tree->GetEvent(evt);
if(tau/taue<5.0) continue;//Condición para el error
if(isnan(taue)) continue;//Si está vacía la entrada que omita ese evento
M = mass;
data.add(RooArgSet(M));
}
data.Print();

tree->ResetBranchAddresses();

// Modelo de la señal
RooRealVar mean("mean", "mean", 6.27, 6.25, 6.3);
RooRealVar sigma("sigma", "sigma", 0.015, 0.01, 0.04);
RooGaussian gauss("gauss", "gauss", M, mean, sigma);

//Modelo del Background
RooRealVar c("c","c",-10.0,10.0);
RooChebychev bkg1("bkg1","Background",M,c);

//Numero de eventos en señal y Background
RooRealVar Ns("Ns","Ns",0.,2000);
RooRealVar Nb("Nb","Nb",0.,2000);   

//Modelo para el fit 
RooAddPdf MassModel("MassModel","MassModel",RooArgList(gauss,bkg1),RooArgList(Ns,Nb));
RooFitResult* fitres = MassModel.fitTo(data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));
data.Print("v"); 
fitres->Print("v");

//Se guardan los parametros del fit en un archivo txt
ofstream salida_TotalFit(Form("To_montecarlo.txt"));
salida_TotalFit.is_open();
salida_TotalFit <<  Ns.getVal() << " " << Ns.getError() << " " << Nb.getVal() << " " << Nb.getError() <<" " <<  c.getVal() << " " 
<< c.getError() << " " << mean.getVal() << " " << mean.getError() << " "<<  sigma.getVal() << " " << sigma.getError();

//Grafica de: Datos, fit, señal, Background
RooPlot *frame = M.frame(Title("Fit Mass"));
data.plotOn(frame,MarkerSize(0.8),XErrorSize(0),Name("Data"));
MassModel.plotOn(frame,LineColor(kBlue),LineWidth(2),Name("fit"));
MassModel.plotOn(frame,Components(gauss),LineColor(kRed),LineWidth(2),Name("signal")); 
MassModel.plotOn(frame,Components(bkg1),LineColor(kGreen+2),LineWidth(2), LineStyle(kDashed) ,Name("bkg1")); 

//Lienzo para las gráficas
TCanvas *c1 = new TCanvas("rf102_dataimport", "rf102_dataimport",50,50,1200,800 );
TPad *pad1 = new TPad("pad1", "padi",0.01,0.411,0.9903769, 0.99 );
pad1->SetLeftMargin(0.09);   
pad1->SetRightMargin(0.019);
pad1->SetTopMargin(0.09);
pad1->SetBottomMargin(0.0);  
TPad *pad2 = new TPad("pad2", "pad2",0.01,0.01,0.9903769,0.41);
pad2->SetLeftMargin(0.09);
pad2->SetRightMargin(0.019);  
pad2->SetTopMargin(0.0);
pad2->SetBottomMargin(0.25);
pad1->Draw();
pad2->Draw();

pad1->cd(); 
frame->SetYTitle("Events"); 
frame->SetTitleSize(0.06,"Y");
frame->SetTitleOffset(0.5,"Y");
frame->GetYaxis()->CenterTitle();   
frame->GetXaxis()->CenterTitle();
frame->GetYaxis()->SetNdivisions(505,1);
frame->GetXaxis()->SetNdivisions(505,1);
frame->GetXaxis()->SetTickLength(0.07); 
frame->SetLabelSize(0.04,"XY");
frame->SetMinimum(-5); 
frame->Draw("A");  

//Leyenda: objetos de la figura
TLegend *leg = new TLegend(0.18,0.58,0.38,0.88); 
leg->SetTextSize(0.06);
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->AddEntry(frame->findObject("Data")," Data","ep"); 
leg->AddEntry(frame->findObject("fit")," Fit result","l");
leg->AddEntry(frame->findObject("bkg1"),"Background.","l");
leg->AddEntry(frame->findObject("signal")," signal","l");
leg->Draw();

//Leyenda: Resultado del fit
TLegend *legpar = new TLegend(0.6,0.68,0.8,0.88);
legpar->SetTextSize(0.06);
legpar->SetFillColor(0);
legpar->SetBorderSize(0);
legpar->SetFillStyle(0);
legpar->AddEntry("",Form("M(B_{c}) = %1.2f #pm %1.2f MeV",mean.getVal()*1000.0, mean.getError()*1000.0 ),"");
legpar->AddEntry("",Form("N_{B_{c}} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
legpar->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
legpar->Draw();

//Pull
RooPlot* Mframe = M.frame(6.05,6.5,((6.5-6.05)/0.008)+1);
data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0));
MassModel.plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));
RooHist* hpullm2 = Mframe->pullHist() ;

pad2->cd();
RooPlot* framem2 = M.frame(Title(" ")) ;
framem2->SetYTitle(" (Data-Fit)/#sigma");
framem2->addPlotable(hpullm2,"P") ;
framem2->GetXaxis()->CenterTitle();
framem2->GetYaxis()->CenterTitle();
framem2->SetXTitle("M(B_{c})[Gev]"); 
framem2->GetYaxis()->SetNdivisions(505,1);
framem2->GetXaxis()->SetNdivisions(505,1);
framem2->GetXaxis()->SetTickLength(0.07);   
framem2->SetTitleOffset(0.35,"Y");
framem2->SetTitleSize(0.08,"Y");
framem2->SetLabelSize(0.06,"XY");
framem2->SetTitleSize(0.08,"X");
framem2->Draw();

//Se almacena el lianzo en la carpeta plots
c1->Draw();
c1->Print("plots/Datos_fit_pull.png");
return 0;
}