
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

// Class generated with tree->MakeClass("paramsToyMc");
#include "paramsToyMc.C"

using namespace std;
using namespace RooFit;


TCanvas* CreateCanvas(TString cname, RooFitResult* result, 
RooDataSet data,  RooRealVar M, Double_t supM, Double_t infM,  
RooGaussian MassModel,  RooRealVar bwg1,  RooRealVar bwm1)  
{



    Double_t nbin = ((supM-infM)/0.2);

    int H = 800;
    int W = 1000;
    
    TCanvas *c1 = new TCanvas(cname,cname,50,50,W,H);
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.07);
    c1->SetTopMargin(0.09);
    c1->SetBottomMargin(0.14); 
    
    RooPlot* Mframe = M.frame(infM,supM,nbin);
    data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.5)); 
    //MassModel.plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(3),Name("fittotal"));// el color rojo es 2
    MassModel.plotOn(Mframe);
    
    Mframe->SetYTitle("Events / 0.2 ");                                                                                                                   
    Mframe->SetLabelSize(0.04,"XY");
    Mframe->SetTitleSize(0.05,"XY");
    Mframe->GetYaxis()->CenterTitle();   
    Mframe->GetXaxis()->CenterTitle();   
    Mframe->SetTitleOffset(1.0,"X");
    Mframe->SetTitleOffset(0.9,"Y");
    Mframe->SetTitleSize(0.06,"XY");
    Mframe->SetMaximum(30.0);
    Mframe->Draw();  
    gStyle->SetOptTitle(0/1);
    
    Double_t G = bwg1.getVal();
    Double_t Ge = bwg1.getError();
    //Double_t G2e = bwg1.getPropagatedError(*result);
    
    TLegend *legpar = new TLegend(0.5,0.7,0.8,0.88);
    legpar->SetTextSize(0.04); //text size in pixels                                 
    legpar->SetFillColor(0);
    legpar->SetBorderSize(0);
    legpar->SetFillStyle(0); 
    legpar->AddEntry("",Form("mean = %1.4f #pm %1.4f GeV ",bwm1.getVal(), bwm1.getError()),"");
    legpar->AddEntry("",Form("#sigma_{1} = %1.4f #pm %1.4f GeV",G, Ge),"");
    legpar->Draw();

    TLegend *legMass = new TLegend(0.64,0.57,0.83,0.65);
    legMass->SetTextFont(43); 
    legMass->SetTextSize(20);  
    legMass->SetFillColor(0); 
    legMass->SetBorderSize(0);
    legMass->SetFillStyle(0); 
    legMass->Draw(); 

    TLatex *   tex1 = new TLatex(0.88,0.926,"ToyMC");
    tex1->SetNDC();
    tex1->SetTextAlign(31);
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.04); 
    tex1->SetLineWidth(2);
    
    TLatex *tex2 = new TLatex(0.15,0.926,"CMS");
    tex2->SetNDC();
    tex2->SetTextFont(61);
    tex2->SetTextSize(0.04); 
    tex2->SetLineWidth(2);
    
    TLatex *tex3 = new TLatex(0.24,0.926,"Preliminary");
    tex3->SetNDC();
    tex3->SetTextFont(52);
    tex3->SetTextSize(0.04); 
    tex3->SetLineWidth(2);
    
    tex1->Draw();  
    tex2->Draw();
    tex3->Draw();
    
    c1->Modified();
    return c1;
  
}


void analysisToyMC()
{
    
    // -------------------------------------------
    // READ TTREE 
    // -------------------------------------------

    Int_t seed = 5;

    TChain *ch = new TChain("tree","");
    ch->Add(Form("ToyMC_Bc_%1i.root/tree", seed));  
    TTree *tree = (TTree*) ch;
    paramsToyMc t(tree); 
    Long64_t nentries = t.fChain->GetEntries();

    cout<<" Tree Entries : "<< nentries <<endl;


    Double_t MminP = -3.0;
    Double_t MmaxP = 3.0;  

    // Data
    RooRealVar Yi("Yi"," B_{c}^{+} yield Pull",MminP,MmaxP); 
    RooDataSet dataYi("dataYi","dataYi",RooArgSet(Yi));

    RooRealVar YiBkg("YiBkg"," B_{c}^{+} yield background Pull",MminP,MmaxP); 
    RooDataSet dataYiBkg("dataYiBkg","dataYiBkg",RooArgSet(YiBkg)); 

    RooRealVar Mu("Mu"," B_{c}^{+} mass Pull",MminP,MmaxP); 
    RooDataSet dataMu("dataMu","dataMu",RooArgSet(Mu)); 



    Int_t nTen = nentries/10;
    Int_t k=0;
    Int_t nbytes = 0, nb = 0;

    for(Long64_t jentry=0; jentry<nentries;jentry++)
        {

        Long64_t ientry = t.LoadTree(jentry);
        if (ientry < 0) break;
        nb = t.fChain->GetEntry(jentry);   nbytes += nb;
        if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
        if(jentry==nentries-1) cout<<endl;

        if(t.status!=0)continue;
        if(t.covQual!=3)continue;

        
        Yi = t.Nspull;
        dataYi.add(RooArgSet(Yi));

        YiBkg = t.Nbpull;
        dataYiBkg.add(RooArgSet(YiBkg));

        Mu = t.Mupull;
        dataMu.add(RooArgSet(Mu));
        
        }



  

    // Model Gaussian Signal 
    RooRealVar meanYi("meanYi","meanYi",0.0,MminP,MmaxP);
    RooRealVar sigmaYi("sigmaYi","sigmaYi",1.0,0.0,5.0);
    RooGaussian SignalYi("SignalYi","SignalYi",Yi, meanYi, sigmaYi); 

    // Model Gaussian Background  
    RooRealVar meanYib("meanYib","meanYib",0.0,MminP,MmaxP);
    RooRealVar sigmaYib("sigmaYib","sigmaYib",1.0,0.0,5.0);
    RooGaussian SigYib("SigYib","SignalYib",YiBkg,meanYib,sigmaYib); 

    // Model Gaussian Mass Mean  
    RooRealVar meanMu("meanMu","meanMu",0.0,MminP,MmaxP);
    RooRealVar sigmaMu("sigmaMu","sigmaMu",1.0,0.0,5.0);
    RooGaussian SigMu("SigMu","SignalMu",Mu, meanMu, sigmaMu); 
    

    // -------------------------------------------
    // FITTING AND CREATE CANVAS
    // -------------------------------------------


    // Fit Signal 
    RooFitResult* fitYi = SignalYi.fitTo(dataYi,Minos(kFALSE),Save(kTRUE), NumCPU(4));
    fitYi->Print("v");

    TCanvas* canv_Yipull = CreateCanvas("canv_Yi", fitYi, dataYi, Yi, MmaxP, MminP, SignalYi, sigmaYi, meanYi);
    canv_Yipull->Print(Form("plots/Pull_YiBc_ToyMC_%1i.png",seed));
    
    // Fit Background 
    RooFitResult* fitYib = SigYib.fitTo(dataYiBkg,Minos(kFALSE),Save(kTRUE), NumCPU(4));
    fitYib->Print("v");

    TCanvas* canv_Yibpull = CreateCanvas("canv_Yib", fitYib, dataYiBkg, YiBkg, MmaxP, MminP, SigYib, sigmaYib, meanYib);
    canv_Yibpull->Print(Form("plots/Pull_YiBkg_Bc_ToyMC_%1i.png",seed));

    // Fit Mass mean 
    RooFitResult* fitMu = SigMu.fitTo(dataMu,Minos(kFALSE),Save(kTRUE), NumCPU(4));
    fitMu->Print("v");

    TCanvas* canv_Mupull = CreateCanvas("canv_Mu", fitMu, dataMu, Mu, MmaxP, MminP, SigMu, sigmaMu, meanMu);
    canv_Mupull->Print(Form("plots/Pull_MuBc_ToyMC_%1i.png",seed));



}