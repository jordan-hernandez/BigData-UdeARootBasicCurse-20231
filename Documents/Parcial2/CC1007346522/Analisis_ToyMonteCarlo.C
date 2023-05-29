
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

#include "ClassToyMC.C"

using namespace std;
using namespace RooFit;

TCanvas* CreateCanvasGauss(TString cname, RooFitResult* result, RooDataSet data,  RooRealVar M, Double_t supM, Double_t infM,  RooGaussian MassModel,  RooRealVar bwg1,  RooRealVar bwm1)  
{
    // Número de bines
    Double_t nbin = ((supM-infM)/0.2);

    int H = 600;
    int W = 800;

    // Creando Canvas
    TCanvas *c1 = new TCanvas(cname,cname,50,50,W,H);
    c1->cd() ;  
    c1->SetLeftMargin(0.1);
    c1->SetRightMargin(0.01);
    c1->SetTopMargin(0.09);
    c1->SetBottomMargin(0.1);
    
    // Creando frame
    RooPlot* Mframe = M.frame(infM,supM,nbin);

    // Graficando datos y modelo
    data.plotOn(Mframe,DataError(RooAbsData::SumW2),XErrorSize(0),MarkerSize(1.5)); 
    MassModel.plotOn(Mframe);
    
    // Ajustes Frame
    Mframe->SetYTitle("Events / 0.2 ");                                                                                                                   
    Mframe->SetLabelSize(0.04,"XY");
    Mframe->SetTitleSize(0.05,"XY");
    Mframe->GetYaxis()->CenterTitle();   
    Mframe->GetXaxis()->CenterTitle();   
    Mframe->SetTitleOffset(0.6,"X");
    Mframe->SetTitleOffset(0.6,"Y");
    Mframe->SetTitleSize(0.06,"XY");
    Mframe->Draw();  
    
    // Valores y errores en valor medio y ancho de la gaussiana
    Double_t G = bwg1.getVal();
    Double_t Ge = bwg1.getError();
    
    // Legend Parámetros
    TLegend *legpar = new TLegend(0.5,0.7,0.8,0.88);
    legpar->SetTextSize(0.04); //text size in pixels                                 
    legpar->SetFillColor(0);
    legpar->SetBorderSize(0);
    legpar->SetFillStyle(0); 
    legpar->AddEntry("",Form("mean = %1.4f #pm %1.4f GeV ",bwm1.getVal(), bwm1.getError()),"");
    legpar->AddEntry("",Form("#sigma = %1.4f #pm %1.4f GeV",G, Ge),"");
    legpar->Draw();

    // Legends de la Figura
    TLatex * tex1 = new TLatex(0.88,0.926,"ToyMC");
    tex1->SetNDC();
    tex1->SetTextAlign(31);
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.04); 
    tex1->SetLineWidth(2);
    
    TLatex *tex2 = new TLatex(0.15,0.926,"UdeA");
    tex2->SetNDC();
    tex2->SetTextFont(61);
    tex2->SetTextSize(0.04); 
    tex2->SetLineWidth(2);
    
    
    tex1->Draw();  
    tex2->Draw();
    
    c1->Modified();
    return c1;
  
}


void Analisis_ToyMonteCarlo()
{
    // Ajustes de estilo
    gStyle->SetOptTitle(0);

    // Creamos el TTree
    TChain *Chain2 = new TChain("tree","");
    Chain2->Add("ToyMC.root/Tree2");
    TTree *tree = (TTree*) Chain2;

    // Usando ClassToyMC desempaquetamos el Tree
    ClassToyMC Tree3(tree); 

    // Número de entradas
    Long64_t nentries = Tree3.fChain->GetEntries();
    cout<<" Entries : "<<nentries<<endl;

    // Límites para los pull
    Double_t MminP = -3.0;
    Double_t MmaxP = 3.0;

    // Variables para los pull y correspondientes Datasets
    RooRealVar Yi("Yi"," B_{c} yield Pull",MminP,MmaxP); 
    RooDataSet dataYi("dataYi","dataYi",RooArgSet(Yi));

    RooRealVar Yib("Yib"," B_{c} yield bkg Pull",MminP,MmaxP); 
    RooDataSet dataYib("dataYib","dataYib",RooArgSet(Yib)); 

    RooRealVar Mu("Mu"," B_{c} mass Pull",MminP,MmaxP); 
    RooDataSet dataMu("dataMu","dataMu",RooArgSet(Mu));   

    // Llenando los Datasets
    Int_t nTen = nentries/10;
    Int_t nbytes = 0, nb = 0;
    for(Long64_t jentry=0; jentry<nentries;jentry++)
        {
        Long64_t ientry = Tree3.LoadTree(jentry);
        if (ientry < 0) break;
        nb = Tree3.fChain->GetEntry(jentry);   nbytes += nb;
        if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
        if(jentry==nentries-1) cout<<endl;

        // Seleccionando datos sin status 0 y covQual 3
        if(Tree3.status!=0)continue;
        if(Tree3.covQual!=3)continue;
        
        // Añadiendo datos a los Datasets
        Yi = Tree3.Nspull;
        dataYi.add(RooArgSet(Yi));

        Yib = Tree3.Nbpull;
        dataYib.add(RooArgSet(Yib));

        Mu = Tree3.Mupull;
        dataMu.add(RooArgSet(Mu));
    }

    // Parámetros de las gaussianas
    RooRealVar meanYi("meanYi","meanYi",0.0,MminP,MmaxP);
    RooRealVar widthYi("widthYi","widthYi",1.0,0.0,5.0);
    RooRealVar meanYib("meanYib","meanYib",0.0,MminP,MmaxP);
    RooRealVar widthYib("widthYib","widthYib",1.0,0.0,5.0);
    RooRealVar meanMu("meanMu","meanMu",0.0,MminP,MmaxP);
    RooRealVar widthMu("widthMu","widthMu",1.0,0.0,5.0);

    // Creando Gaussianas
    RooGaussian SigYi("SigYi","SignalYi",Yi,meanYi,widthYi); 
    RooGaussian SigYib("SigYib","SignalYib",Yib,meanYib,widthYib);
    RooGaussian SigMu("SigMu","SignalMu",Mu,meanMu,widthMu); 
    
    // Fiteando los pull
    RooFitResult* fitYi = SigYi.fitTo(dataYi,Minos(kFALSE),Save(kTRUE), NumCPU(4));
    RooFitResult* fitYib = SigYib.fitTo(dataYib,Minos(kFALSE),Save(kTRUE), NumCPU(4));
    RooFitResult* fitMu = SigMu.fitTo(dataMu,Minos(kFALSE),Save(kTRUE), NumCPU(4));

    // Creando Canvas para cada variable
    TCanvas* canv_Yipull = CreateCanvasGauss("canv_Yi", fitYi, dataYi, Yi, MmaxP, MminP, SigYi, widthYi, meanYi);
    canv_Yipull->Print("plots/Pull_YiBs_ToyMC.png");    

    TCanvas* canv_Yibpull = CreateCanvasGauss("canv_Yib", fitYib, dataYib, Yib, MmaxP, MminP, SigYib, widthYib, meanYib);
    canv_Yibpull->Print("plots/Pull_Yib_ToyMC.png"); 

    TCanvas* canv_Mupull = CreateCanvasGauss("canv_Mu", fitMu, dataMu, Mu, MmaxP, MminP, SigMu, widthMu, meanMu);
    canv_Mupull->Print("plots/Pull_Mu_ToyMC.png");
}