#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooExponential.h>

using namespace RooFit;
 
void Twodim()
{
    Float_t bcm,tau,etau;
    //Lectura de datos de filtrados
    TFile f ("partc2.root", "READ");
    TTree *t = (TTree*) f.Get ("mytree");
    t->SetBranchAddress("bcm",&bcm);
    t->SetBranchAddress("tau", &tau);
    t->SetBranchAddress("etau", &etau);

    RooRealVar Bcm("Bcm","Bcm",6.05,6.5);
    RooRealVar Tau("Tau","Tau",0.21,2.6);

    //Dataset
    RooDataSet data("data","data",RooArgSet(Bcm,Tau));

    for (int evt = 0; evt <t->GetEntries();evt++)
    {
    t->GetEvent(evt);
    if(tau/etau<5.0) continue;
    Bcm = bcm;
    Tau = tau;
    data.add(RooArgSet(Bcm, Tau));
    }

    //Modelo de masa
    //Gaussiana
    RooRealVar mean("mean", "mean", 6.27, 6.25, 6.3);
    RooRealVar sigma1("sigma1", "sigma1", 0.015, 0.01, 0.04);
    RooGaussian sig1("sig1", "sig1", Bcm, mean, sigma1);
    //Background
    RooRealVar a0("a0","a0",-10,10);
    RooChebychev bkg("bkg","Background",Bcm,a0);
    //Pesos
    RooRealVar nsig("nsig","nsig",0.,1500);
    RooRealVar nbkg("nbkg","nbkg",0.,1500);   
    RooAddPdf modelM("modelM","modelM",RooArgList(sig1,bkg),RooArgList(nsig,nbkg));

    //Lifetime
    RooRealVar ctau("ctau", "ctau", -5,1);
    RooExponential decay("decay", "decay", Tau, ctau);

    //Producto de modelo de masa y modelo de tiempo de vida
    RooProdPdf BcmTau("BcmTau", "BcmTau", modelM, Conditional(decay, Tau));

    //Ajuste de los datos al modelo masa * tiempo de vida
    BcmTau.fitTo(data,Extended(kTRUE),Save(kTRUE)) ; 

    RooPlot* frame1 = Bcm.frame(Title("Bcm"));

    RooPlot* frame = Tau.frame(Title("Tau"));
    // Histograma 2D
    TH1 *hh_decay = BcmTau.createHistogram("2d", Bcm, YVar(Tau));
    hh_decay->SetLineColor(kViolet);

   // Draw all frames on canvas
   TCanvas *c = new TCanvas("canvas", "canvas", 1800, 600);
   c->Divide(3);
   c->cd(1);
   gPad->SetLeftMargin(0.20);
   hh_decay->GetXaxis()->SetTitleOffset(1.5);
   hh_decay->GetZaxis()->SetTitleOffset(2.1);
   hh_decay->SetTitle("B_{c}^{+} and Decay");
   hh_decay->GetYaxis()->CenterTitle();   
   hh_decay->GetXaxis()->CenterTitle();
   hh_decay->GetZaxis()->CenterTitle();
   hh_decay->Draw("surf");

   c->cd(2);
   gPad->SetLeftMargin(0.15);
   frame->SetTitle("Decay");
   frame->SetXTitle("tau [ps]");
   frame->SetYTitle("Events"); 
   frame->GetYaxis()->CenterTitle();   
   frame->GetXaxis()->CenterTitle();
   data.plotOn(frame,Name("Data"));
   BcmTau.plotOn(frame,LineWidth(2),Name("fit"), LineColor(kBlue));
   frame->Draw();

   TLegend *leg = new TLegend(0.67,0.47,0.87,0.77); 
   leg->SetTextSize(0.045);
   leg->SetBorderSize(0);
   leg->AddEntry(frame->findObject("Data")," Data"); 
   leg->AddEntry(frame->findObject("fit")," Fit result");
   leg->Draw();

   c->cd(3);
   gPad->SetLeftMargin(0.15);
   frame1->SetTitle("B_{c}^{+}");
   frame1->SetYTitle("Events"); 
   frame1->SetXTitle("M(J/\\psi\\pi^{+}) [GeV]");
   frame1->GetYaxis()->CenterTitle();   
   frame1->GetXaxis()->CenterTitle();
   data.plotOn(frame1,Name("Data"));
   BcmTau.plotOn(frame1, LineColor(kRed),Name("fit"));
   frame1->Draw();

   TLegend *leg1 = new TLegend(0.67,0.47,0.87,0.77); 
   leg1->SetTextSize(0.045);
   leg1->SetBorderSize(0);
   leg1->AddEntry(frame1->findObject("Data")," Data"); 
   leg1->AddEntry(frame1->findObject("fit")," Fit result");
   leg1->Draw();
   c->Print("2d.png");

}




