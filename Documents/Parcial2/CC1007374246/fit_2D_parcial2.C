
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooGenericPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooAbsCollection.h"
using namespace RooFit;

void fit_2D_parcial2()
{
Float_t mass , tau,taue;

//Lectura de datos de filtrados
TFile f ("Datos_organizados.root", "READ");
TTree *tree = (TTree*) f.Get ("mytree");
tree->SetBranchAddress("mass",&mass);
tree->SetBranchAddress("tau", &tau);
tree->SetBranchAddress("taue", &taue);

RooRealVar M("M","Mass",6.05,6.5);
RooRealVar Tau("Tau","Tau",0.3,2.6);

//Dataset
RooDataSet data("data","data",RooArgSet(M,Tau));
Int_t nentries = tree->GetEntries();
cout<<"Entradas: "<<nentries<<endl;

for (int evt = 0; evt <tree->GetEntries();evt++)
{
tree->GetEvent(evt);
if(isnan(taue)) continue;//Filtro: Datos vacios
if(tau/taue<5.0) continue;//Filtro: condicion sobre tau y taue
if(mass<6.05 || mass>6.5) continue;//Filtro: condicion sobre mass
if(tau<0.3 || tau>2.6) continue;//Filtro: condicion sobre tau 
if(taue<0.0001 || taue>0.5) continue;//Filtro: condicion sobre taue

M = mass;
Tau = tau;
data.add(RooArgSet(M, Tau));
}

data.Print();

//Modelo de masa
//Gaussiana
RooRealVar mean("mean", "mean", 6.27, 6.25, 6.3);
RooRealVar sigma("sigma", "sigma", 0.015, 0.01, 0.04);
RooGaussian gauss("gauss", "gauss", M, mean, sigma);
//Background
RooRealVar c("c","c",-10.0,10.0);
RooChebychev bkg1("bkg1","Background",M,c);
//Pesos
RooRealVar Ns("Ns","Ns",0.,1000);
RooRealVar Nb("Nb","Nb",0.,1000);   
RooAddPdf MassModel("MassModel","MassModel",RooArgList(gauss,bkg1),RooArgList(Ns,Nb));

//Lifetime
RooRealVar c2("c2", "c2", -10,10);
RooExponential lifetime("lifetime", "lifetime", Tau, c2);

//Producto de modelo de masa y modelo de tiempo de vida
RooProdPdf MassandLifetime("massandLifetime", "massandLifetime", MassModel, Conditional(lifetime, Tau));

//Ajuste de los datos al modelo masa * tiempo de vida
RooFitResult* fitMT =MassandLifetime.fitTo(data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4)) ; 
fitMT->Print("v");
RooPlot *frame = Tau.frame(Title("Tau"));

//Lienzo
TCanvas* c1 = new TCanvas("c", "c", 800, 600);
TPad *pad1 = new TPad("p1", "", 0.01,0.05,0.49, 0.99);
pad1->SetLeftMargin(0.09);   
pad1->SetRightMargin(0.019);
pad1->SetTopMargin(0.09);
pad1->Draw();
TPad *pad2 = new TPad("p2", "", 0.51,0.05,0.99,0.99);
pad2->SetLeftMargin(0.09);
pad2->SetRightMargin(0.019);  
pad2->SetTopMargin(0.09);
pad2->Draw();
pad1->cd();
RooPlot* frame1 = M.frame(Title("Tau"));
frame1->SetTitle("Mass B_{c} ");
frame1->SetYTitle("Events"); 
frame1->SetXTitle("M(B_{c})[Gev]");
frame1->GetYaxis()->CenterTitle();   
frame1->GetXaxis()->CenterTitle();
frame1->GetYaxis()->SetNdivisions(505,1);
frame1->GetXaxis()->SetNdivisions(505,1);
frame1->GetXaxis()->SetTickLength(0.07); 
data.plotOn(frame1,Name("Data"),MarkerSize(0.8));
MassandLifetime.plotOn(frame1, LineColor(kGreen +2),LineWidth(2),Name("fit"),XErrorSize(0));
frame1->Draw();

TLegend *leg1 = new TLegend(0.18,0.58,0.38,0.88); 
leg1->SetTextSize(0.04);
leg1->SetFillColor(0);
leg1->SetBorderSize(0);
leg1->SetFillStyle(0);
leg1->AddEntry(frame1->findObject("Data")," Data","ep"); 
leg1->AddEntry(frame1->findObject("fit")," Fit result","l");
leg1->Draw();

// auto mtext = new TLatex();
// mtext->SetTextSize(0.031);
// mtext->SetTextFont(42);
// mtext->DrawLatex(6.3, 45,Form("B_{c Mass} = %1.4f #pm %1.4f GeV", massMean->getVal(), massMean->getError()));

pad2->cd();
RooPlot* frame3 = Tau.frame();
frame3->SetTitle("Life time B_{c}");
frame3->SetXTitle("tau [ps]");
frame3->SetYTitle("Events"); 
frame3->GetYaxis()->CenterTitle();   
frame3->GetXaxis()->CenterTitle();
frame3->GetYaxis()->SetNdivisions(505,1);
frame3->GetXaxis()->SetNdivisions(505,1);
frame3->GetXaxis()->SetTickLength(0.07); 
data.plotOn(frame3,MarkerSize(0.8),XErrorSize(0),Name("Data"));
MassandLifetime.plotOn(frame3,LineWidth(2),Name("fit"), LineColor(kGreen+2),MarkerSize(0.8),XErrorSize(0));
frame3->Draw();

TLegend *leg = new TLegend(0.18,0.58,0.38,0.88); 
leg->SetTextSize(0.04);
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->AddEntry(frame3->findObject("Data")," Data","ep"); 
leg->AddEntry(frame3->findObject("fit")," Fit result","l");
leg->Draw();

c1->Draw();
c1->Print("plots/Fit_M_Tau.png");

// Histograma 2D
TH1 *hh = MassandLifetime.createHistogram(" ", M, Binning(50), YVar(Tau, Binning(50)));
hh->SetLineColor(kBlue);


// Grafica
RooPlot *frame2 = M.frame(Title(" titulo"));
data.plotOn(frame2);
MassandLifetime.plotOn(frame2, ProjWData(data, true));

// Lienzo
TCanvas* c3 = new TCanvas("rf306_condpereventerrors", "rf306_condperventerrors", 1200, 400);

gPad->SetLeftMargin(0.20);
hh->GetXaxis()->SetTitleOffset(2.5);
hh->GetYaxis()->SetTitleOffset(4.5);
hh->GetZaxis()->SetTitleOffset(1.5);
hh->SetTitle("Life time B_{c}");
hh->GetYaxis()->CenterTitle();   
hh->GetXaxis()->CenterTitle();
hh->GetZaxis()->CenterTitle();
hh->GetYaxis()->SetNdivisions(505,1);
hh->GetXaxis()->SetNdivisions(505,1);
hh->GetZaxis()->SetNdivisions(505,1);
hh->GetXaxis()->SetTickLength(0.07); 
hh->SetTitleOffset(0.5,"Y");
hh->Draw("surf");

c3->Draw();
c3->Print("plots/2D.png");

}