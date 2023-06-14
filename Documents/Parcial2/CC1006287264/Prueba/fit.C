#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <iostream>
#include <fstream>

using namespace RooFit;
using namespace std;

void fit() {
    Float_t Bcmass, tau, tauerror;
    TFile f("bc2_data.root", "READ");
    TTree *tree = (TTree*) f.Get("tree");
    tree->SetBranchAddress("Bcmass", &Bcmass);
    tree->SetBranchAddress("tau", &tau);
    tree->SetBranchAddress("tauerror", &tauerror);

    RooRealVar M("M", "M", 6.05, 6.5);
    RooDataSet data("dataSet", "dataSet", RooArgSet(M));

    // Se crea dataset para la masa
    for (int event = 0; event < tree->GetEntries(); event++) {
        tree->GetEntry(event);
        if (Bcmass <= 6.0 || Bcmass >= 6.5) continue;
        if ((tau / tauerror) < 5) continue;
        if (tau < 0.3 || tau > 2.6) continue;
        if (tauerror < 0.0001 || tauerror > 0.5) continue;
        M.setVal(Bcmass);
        data.add(RooArgSet(M));
    }

    // Gaussiana de la distribución de la masa
    RooRealVar mean("mean", "mean", 6.2751, 6.25, 6.3);
    RooRealVar sigma("sigma", "sigma", 0.015, 0.01, 0.04);
    RooGaussian gaus("gaus", "gaus", M, mean, sigma);

    // Línea del background
    RooRealVar a0("a0", "a0", -10, 10);
    //RooRealVar a1("a1", "a1", -10, 10);
    //RooChebychev bkg("bkg", "bkg", M, RooArgList(a0, a1));
    RooChebychev bkg("bkg", "bkg", M, a0);
    
    // Número de eventos
    RooRealVar nsig("nsig", "nsig", 0., 1500);
    RooRealVar nbkg("nbkg", "nbkg", 0., 1500);

    RooAddPdf massmodel("massmodel", "massmodel", RooArgList(gaus, bkg), RooArgList(nsig, nbkg));
    RooFitResult* fitres = massmodel.fitTo(data, Extended(kTRUE), Save(kTRUE), NumCPU(4));

    // Guardamos los datos resultado del fit en un archivo para luego usarlos en el montecarlo.
    ofstream myfile;
    myfile.open("fit.txt");
   /* myfile << nsig.getVal() << " " << nsig.getError() << " " << nbkg.getVal() << " " << nbkg.getError()
           << " " << mean.getVal() << " " << mean.getError() << " " << sigma.getVal() << " " << sigma.getError()
           << " " << a0.getVal() << " " << a0.getError() << " " << a1.getVal() << " " << a1.getError();
    myfile.close();
   */
    myfile << nsig.getVal() << " " << nsig.getError() << " " << nbkg.getVal() << " " << nbkg.getError()
           << " " << a0.getVal() << " " << a0.getError() << " " << mean.getVal() << " " << mean.getError()
           << " " << sigma.getVal() << " " << sigma.getError() ;
    myfile.close();
    
    // Gráfico
    TCanvas* c = new TCanvas("c", "c", 1000, 800);

    TPad *pad1 = new TPad("p1", "", 0.01, 0.411, 0.9903769, 0.99);
    pad1->SetLeftMargin(0.09);
    pad1->SetRightMargin(0.019);
    pad1->SetTopMargin(0.09);
    pad1->SetBottomMargin(0.0);
    pad1->Draw();

    TPad *pad2 = new TPad("p2", "", 0.01, 0.01, 0.9903769, 0.41);
    pad2->SetLeftMargin(0.09);
    pad2->SetRightMargin(0.019);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.25);
    pad2->Draw();

    pad1->cd();
    RooPlot* mframe = M.frame(6.05, 6.5, 60);
    data.plotOn(mframe, Name("data"));

    mframe->SetYTitle("Events / 10 MeV");
    mframe->SetLabelSize(0.05, "XY");
    mframe->SetTitleSize(0.08, "X");
    mframe->SetTitleSize(0.05, "Y");
    mframe->GetYaxis()->CenterTitle();

    massmodel.plotOn(mframe, Components(gaus), LineColor(kViolet-6), LineStyle(kDashed), Name("signal"));
    massmodel.plotOn(mframe, Components(bkg), LineColor(kGreen), LineStyle(kDashed), Name("bkg"));
    massmodel.plotOn(mframe, LineColor(kBlue-4), Name("model"));

    TLegend *legMass = new TLegend(0.18, 0.58, 0.38, 0.88);
    legMass->SetTextSize(0.08);
    legMass->SetFillColor(0);
    legMass->SetBorderSize(0);
    legMass->SetFillStyle(0);
    legMass->SetTextSize(0.06);
    legMass->AddEntry(mframe->findObject("data"), "Data", "pe");
    legMass->AddEntry(mframe->findObject("model"), "Fit", "l");
    legMass->AddEntry(mframe->findObject("signal"), "Signal", "l");
    legMass->AddEntry(mframe->findObject("bkg"), "Bkg", "l");

    TLatex *tex1 = new TLatex(0.88, 0.95, "L = 19.7 fb^{-1} (#sqrt{s} = 8 TeV)");
    tex1->SetNDC();
    tex1->SetTextAlign(31);
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.04);
    tex1->SetLineWidth(2);

    TLatex *tex2 = new TLatex(0.15, 0.95, "CMS");
    tex2->SetNDC();
    tex2->SetTextFont(61);
    tex2->SetTextSize(0.04);
    tex2->SetLineWidth(2);

    mframe->SetTitle("");
    mframe->Draw();
    legMass->Draw();

    TLegend *legfit = new TLegend(0.6, 0.68, 0.8, 0.88);
    legfit->SetTextSize(0.05);
    legfit->SetFillColor(0);
    legfit->SetBorderSize(0);
    legfit->SetFillStyle(0);

    legfit->AddEntry("", Form("M(B_{c}^{+})= %1.3f #pm %1.3fMeV", mean.getVal() * 1000, mean.getError() * 1000), "");
    legfit->AddEntry("", Form("#sigma = %1.3f #pm %1.3f MeV", sigma.getVal() * 1000, sigma.getError() * 1000), "");

    tex1->Draw();
    tex2->Draw();
    legfit->Draw();

    // Construcción del pull

    pad2->cd();
    RooHist* massPull = mframe->pullHist();
    RooPlot* frame2 = M.frame(Title(" "));
    frame2->addPlotable(massPull, "P");
    frame2->SetYTitle(" (Data-Fit)/#sigma");
    frame2->GetXaxis()->CenterTitle();
    frame2->GetYaxis()->CenterTitle();
    frame2->SetXTitle("M(J/#psi#pi^{+}) [GeV]");
    frame2->GetYaxis()->SetNdivisions(505, 1);
    frame2->GetXaxis()->SetNdivisions(505, 1);
    frame2->GetXaxis()->SetTickLength(0.07);
    frame2->SetTitleOffset(0.35, "Y");
    frame2->SetTitleSize(0.08, "Y");
    frame2->SetLabelSize(0.06, "XY");
    frame2->SetTitleSize(0.08, "X");
    frame2->Draw();

    c->Draw();
    c->SaveAs("mass_fit.svg");
}






/*
using namespace RooFit;
using namespace std;
void fit(){
    float Bcmass, tau, tauerror;
    TFile f("bc2_data.root", "READ");
    TTree *tree = (TTree*) f.Get("tree");
    tree->SetBranchAddress("Bcmass", &Bcmass);
    tree->SetBranchAddress("tau", &tau);
    tree->SetBranchAddress("tauerror", &tauerror);

    RooRealVar M("M", "M", 6.05, 6.5); 
    RooDataSet data("dataSet", "dataSet", RooArgSet(M));

    //Se crea dataset para la masa
    for(int event = 0; event < tree->GetEntries(); event++){
        tree->GetEntry(event);
        if(Bcmass<=6.0 || Bcmass>=6.5) continue;
        if((tau/tauerror)<5) continue;
        if(tau<0.3 || tau>2.6) continue;
        if(tauerror<0.0001 || tauerror>0.5) continue;
        M.setVal(Bcmass);
        data.add(RooArgSet(M));
    }   
    

    //gaussiana de la distribución de la masa
    RooRealVar mean("mean", "mean", 6.2751,6.25,6.3);
    RooRealVar sigma("sigma", "sigma", 0.015, 0.01, 0.04);
    RooGaussian gaus("gaus", "gaus", M, mean, sigma);

    
    //línea del background
    RooRealVar a0("a0", "a0", -10,10);
    RooRealVar a1("a1", "a1", -10,10);
    RooChebychev bkg("bkg", "bkg", M, RooArgList(a0, a1));

    //número de eventos
    RooRealVar nsig("nsig","nsig",0.,13000);
    RooRealVar nbkg("nbkg","nbkg",0.,13000);


    RooAddPdf massmodel("massmodel","massmodel",RooArgList(gaus,bkg),RooArgList(nsig,nbkg));
    RooFitResult* fitres = massmodel.fitTo(data,Extended(kTRUE),Save(kTRUE), NumCPU(4));

    //guardamos los datos resultado del fit en un archivo para luego usarlos en el montecarlo. 
    ofstream myfile;
    myfile.open ("fit.txt");
    myfile << nsig.getVal() << " " << nsig.getError()<<" " << nbkg.getVal() << " " << nbkg.getError()<<" "<<mean.getVal() << " " << mean.getError()<< " "<<sigma.getVal() << " " << sigma.getError();
    myfile <<" " <<a0.getVal() << " " << a0.getError()<< " " <<a1.getVal() << " " << a1.getError();
    myfile.close();


    //Grafico
    TCanvas* c = new TCanvas("c", "c", 1000, 800);

    TPad *pad1 = new TPad("p1", "", 0.01,0.411,0.9903769, 0.99 );
    pad1->SetLeftMargin(0.09);   
    pad1->SetRightMargin(0.019);
    pad1->SetTopMargin(0.09);
    pad1->SetBottomMargin(0.0);  
    pad1->Draw();

    TPad *pad2 = new TPad("p2", "", 0.01,0.01,0.9903769,0.41);
    pad2->SetLeftMargin(0.09);
    pad2->SetRightMargin(0.019);  
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.25);
    pad2->Draw();

    pad1->cd();
    RooPlot* mframe = M.frame(6.05, 6.5, 60);
    data.plotOn(mframe, Name("data"));

    mframe->SetYTitle("Events / 10 MeV");
    
    mframe->SetLabelSize(0.05,"XY");
    mframe->SetTitleSize(0.08,"X");
    mframe->SetTitleSize(0.05,"Y"); 
    mframe->GetYaxis()->CenterTitle();  
    
    massmodel.plotOn(mframe, Components(gaus), LineColor(kViolet-6), LineStyle(kDashed), Name("signal"));
    massmodel.plotOn(mframe, Components(bkg), LineColor(kGreen), LineStyle(kDashed), Name("bkg"));
    massmodel.plotOn(mframe, LineColor(kBlue-4), Name("model"));

    TLegend *legMass = new TLegend(0.18,0.58,0.38,0.88);
    legMass->SetTextSize(0.08);
    legMass->SetFillColor(0);
    legMass->SetBorderSize(0);
    legMass->SetFillStyle(0);
    legMass->SetTextSize(0.06); 
    legMass->AddEntry(mframe->findObject("data"), "Data", "pe");
    
    legMass->AddEntry(mframe->findObject("model"), "Fit", "l");
    legMass->AddEntry(mframe->findObject("signal"),"Signal","l");
    legMass->AddEntry(mframe->findObject("bkg"),"Bkg","l");

    
    
    
    TLatex *   tex1 = new TLatex(0.88,0.95,"L = 19.7 fb^{-1} (#sqrt{s} = 8 TeV)");
    tex1->SetNDC();
    tex1->SetTextAlign(31);
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.04); 
    tex1->SetLineWidth(2);
  
      TLatex *tex2 = new TLatex(0.15,0.95,"CMS");
      tex2->SetNDC();
      tex2->SetTextFont(61);
      tex2->SetTextSize(0.04); 
      tex2->SetLineWidth(2);
   
  
    mframe->SetTitle("");
    mframe->Draw();
    legMass->Draw();
    TLegend *legfit = new TLegend(0.6,0.68,0.8,0.88); 
    legfit->SetTextSize(0.05);
    legfit->SetFillColor(0);
    legfit->SetBorderSize(0);
    legfit->SetFillStyle(0);
   
    legfit->AddEntry("",Form("M(B_{c}^{+})= %1.3f #pm %1.3fMeV",mean.getVal()*1000,mean.getError()*1000),"");
    legfit->AddEntry("",Form("#sigma = %1.3f #pm %1.3f MeV",sigma.getVal()*1000,sigma.getError()*1000),"");
    
    tex1->Draw();
    tex2->Draw(); 
    legfit->Draw();

    //construcción del pull

    pad2->cd();
    RooHist* massPull = mframe->pullHist();
    RooPlot* frame2 = M.frame(Title(" ")) ;
    frame2->addPlotable(massPull,"P") ;
    frame2->SetYTitle(" (Data-Fit)/#sigma");
    frame2->GetXaxis()->CenterTitle();
    frame2->GetYaxis()->CenterTitle();
    frame2->SetXTitle("M(J/#psi#pi^{+}) [GeV]"); 
    frame2->GetYaxis()->SetNdivisions(505,1);
    frame2->GetXaxis()->SetNdivisions(505,1);
    frame2->GetXaxis()->SetTickLength(0.07);   
    frame2->SetTitleOffset(0.35,"Y");
    frame2->SetTitleSize(0.08,"Y");
    frame2->SetLabelSize(0.06,"XY");
    frame2->SetTitleSize(0.08,"X");
    frame2->Draw();

    c->Draw();
    c->SaveAs("canvas.svg");

}*/
