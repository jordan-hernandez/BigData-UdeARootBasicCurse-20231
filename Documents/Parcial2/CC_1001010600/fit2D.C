#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
using namespace RooFit;
void fit2D(){
    Float_t mass, tau, tErr;

    TFile f("BcData.root", "READ");
    TTree *tree = (TTree*) f.Get("tree");

    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("tau", &tau);
    tree->SetBranchAddress("tErr", &tErr);

    RooRealVar MASS("MASS", "Mass (GeV)", 6.05, 6.5); //Observable
    RooRealVar TAU("TAU", "tau (ps)", .3, 2.6); 
    RooRealVar tERR("tERR", "tau_{e} (ps)", .0001, .5); 

    RooDataSet massDS("mdataSet", "mdataSet", RooArgSet(MASS));
    RooDataSet tauDS("tdataSet", "tdataSet", RooArgSet(TAU));
    RooDataSet data("data", "data", RooArgSet(MASS, TAU));

    //Se crea dataset para la masa
    for(int event = 0; event < tree->GetEntries(); event++){
        tree->GetEntry(event);
        if(mass<=6.0 || mass>=6.5) continue;
        if(tau<=.3 || tau>=2.6) continue;
        if(tErr<=.0001 || tErr>=.5) continue;

        if((tau/tErr)<5.0) continue;

        MASS = mass;
        TAU = tau;
        tERR = tErr;
        data.add(RooArgSet(MASS, TAU));
    }

    //Masa
    RooRealVar mean1("mean1", "mean1", 6.27, 6.25, 6.3);
    RooRealVar sigma1("sigma1", "sigma1", .015, 0.005, 0.05);
    RooGaussian sig1("sig1", "sig1", MASS, mean1, sigma1);

    RooRealVar a0("a0", "a0", -10,10);
    RooRealVar a1("a1", "a1", -10,10);

    RooChebychev bkg("bkg", "bkg", MASS, RooArgList(a0, a1));

    //Numero de eventos
    RooRealVar nsig("nsig","nsig",1000,0.,13000);
    RooRealVar nbkg("nbkg","nbkg",2000,0.,13000);
    
    //Pdf total para la masa
    RooAddPdf masa("model","sig+bkg", RooArgList(sig1,bkg), RooArgList(nsig,nbkg));

    //Tiempo de vida
    RooRealVar c("c", "c", -10,10);
    RooExponential lifet("lifet", "lifet", TAU, c);

    RooProdPdf massLifeTime("massLifeTime", "massLifeTime", masa, Conditional(lifet, TAU));

    RooFitResult* result = massLifeTime.fitTo(data, ConditionalObservables(tERR), Save(kTRUE));
    auto pars = result->floatParsFinal();
    RooRealVar *massMean; massMean = (RooRealVar *)pars.at(3);
    RooRealVar *lifetime; lifetime = (RooRealVar *)pars.at(2);


    TCanvas* canvas = new TCanvas("c", "c", 1000, 600);

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
    RooPlot* framem = MASS.frame();
    framem->SetTitle("Mass for B_{c} meson");
    framem->SetXTitle("Mass (GeV)");

    data.plotOn(framem);
    massLifeTime.plotOn(framem, LineColor(kMagenta-3));
    framem->Draw();

    auto mtext = new TLatex();
    mtext->SetTextSize(0.031);
    mtext->SetTextFont(42);
    mtext->DrawLatex(6.3, 45,Form("B_{c Mass} = %1.4f #pm %1.4f GeV", massMean->getVal(), massMean->getError()));

    pad2->cd();
    RooPlot* framet = TAU.frame();
    framet->SetTitle("Life time for B_{c} meson");
    framet->SetXTitle("#tau (ps)");

    data.plotOn(framet);
    massLifeTime.plotOn(framet, LineColor(kMagenta-3));
    framet->Draw();

    auto lftext = new TLatex();
    lftext->SetTextSize(0.031);
    lftext->SetTextFont(42);
    lftext->DrawLatex(1.6, 45,Form("B_{c #tau} = %1.4f #pm %1.4f ps", -1/lifetime->getVal(), -lifetime->getError()/lifetime->getVal()));

    canvas->Draw();
    canvas->SaveAs("./plots/2Dfit.png");

    
}