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

#include "tree1.C"

using namespace std;
using namespace RooFit;

TCanvas* CreateCanvasGauss(TString cname, RooFitResult* result, RooDataSet data,  RooRealVar M, Double_t supM, Double_t infM,  RooGaussian MassModel,  RooRealVar bwg1,  RooRealVar bwm1)  
{

//   Double_t nbin = ((supM-infM)/0.2);
    Double_t nbin = 50;

  int H = 600;
  int W = 800;
  
  TCanvas *c1 = new TCanvas(cname,cname,50,50,W,H);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.07);
  c1->SetTopMargin(0.09);
  c1->SetBottomMargin(0.14); 
  
  RooPlot* Mframe = M.frame(infM,supM,nbin);
  data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.5)); 
  MassModel.plotOn(Mframe);
 
  Mframe->SetYTitle("Events");                                                                                                                   
  Mframe->SetLabelSize(0.04,"XY");
  Mframe->SetTitleSize(0.05,"XY");
  Mframe->GetYaxis()->CenterTitle();   
  Mframe->GetXaxis()->CenterTitle();   
  Mframe->SetTitleOffset(1.0,"X");
  Mframe->SetTitleOffset(0.9,"Y");
  Mframe->SetTitleSize(0.06,"XY");
  Mframe->SetMaximum(10000.0);
  Mframe->Draw();  
 
  Double_t G = bwg1.getVal();
  Double_t Ge = bwg1.getError();
  //Double_t G2e = bwg1.getPropagatedError(*result);
  
  TLegend *legpar = new TLegend(0.55,0.7,0.8,0.88);  //(0.5,0.7,0.8,0.88)
  legpar->SetTextSize(0.04); //text size in pixels                                 
  legpar->SetFillColor(0);
  legpar->SetBorderSize(0);
  legpar->SetFillStyle(0); 
  legpar->AddEntry("",Form("#mu = %1.4f #pm %1.4f GeV ",bwm1.getVal(), bwm1.getError()),"");
  legpar->AddEntry("",Form("#sigma = %1.4f #pm %1.4f GeV",G, Ge),"");
  legpar->Draw();

  TLegend *legMass = new TLegend(0.64,0.57,0.83,0.65);
  legMass->SetTextFont(43); 
  legMass->SetTextSize(20);  
  legMass->SetFillColor(0); 
  legMass->SetBorderSize(0);
  legMass->SetFillStyle(0); 
//   legMass->SetHeader(Form("%1.1f #leq p_{T}(B_{s}^{0}) < %1.1f GeV ",ptl,pth));
  legMass->Draw(); 

  //TLatex *   tex1 = new TLatex(0.88,0.926,"L = 78.2 fb^{-1} (#sqrt{s} = 13 TeV)");
  TLatex * tex1 = new TLatex(0.88,0.926,"ToyMC");
  tex1->SetNDC();
  tex1->SetTextAlign(31);
  tex1->SetTextFont(42);
  tex1->SetTextSize(0.04); 
  tex1->SetLineWidth(2);
  
  
  tex1->Draw();  

  
  c1->Modified();
  return c1;
  
}


void mz_anali(Int_t seed = 4, Double_t ptl=12, Double_t pth=70)
{

 gStyle->SetOptTitle(0);
 //setTDRStyle();

TChain *ch = new TChain("tree","");
//ch->Add("ToyMC_Bc2s_1.root/tree");
ch->Add("ToyMC.root/tree1");  
TTree *tree = (TTree*) ch;
tree1 t(tree); 
Long64_t nentries = t.fChain->GetEntries();

cout<<" Entries : "<<nentries<<endl;

using namespace RooFit;
using namespace std;

Double_t MminP = -5.0;
Double_t MmaxP = 5.0;  

// RooRealVar Yi("Yi"," B_{c}^{+} yield Pull",MminP,MmaxP); 
// RooDataSet dataYi("dataYi","dataYi",RooArgSet(Yi));

// RooRealVar Yib("Yib"," B_{c}^{+} yield bkg Pull",MminP,MmaxP); 
// RooDataSet dataYib("dataYib","dataYib",RooArgSet(Yib)); 

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
      if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"-"<<flush;
      if(jentry==nentries-1) cout<<endl;

      if(t.status!=0)continue;
      if(t.covQual!=3)continue;

      Mu = t.Mupull;
      dataMu.add(RooArgSet(Mu));

      
    }
 

// RooRealVar meanMu("meanMu","meanMu",0.0,MminP,MmaxP);
RooRealVar meanMu("meanMu","meanMu",0.0,-1.0,1.0);
RooRealVar widthMu("widthMu","widthMu",1.0,0.0,5.0);
RooGaussian SigMu("SigMu","SignalMu",Mu,meanMu,widthMu); 
 
RooFitResult* fitMu = SigMu.fitTo(dataMu,Minos(kFALSE),Save(kTRUE), NumCPU(4));
// fitMu->Print("v");

TCanvas* canv_Mupull = CreateCanvasGauss("canv_Mu", fitMu, dataMu, Mu, MmaxP, MminP, SigMu, widthMu, meanMu);
canv_Mupull->Print("pullMC.png");


 
}