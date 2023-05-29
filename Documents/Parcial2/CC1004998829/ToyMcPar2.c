#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TChain.h>

#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooArgusBG.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "RooVoigtian.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"

#include "RooCBShape.h"

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
#include "RooRandom.h"

using namespace RooFit;
using namespace std;
/*
TCanvas* CreateCanvasNom(TString cname, RooFitResult* result, RooDataSet *datatre, RooRealVar Bcm, Double_t supM, Double_t infM,  RooAddPdf model, RooGaussian sig1, RooChebychev bkg, RooRealVar nsig, RooRealVar nbkg, RooRealVar sigma1, RooRealVar mean, Double_t ptl, Double_t pth)
{

  //Double_t nbin = ((supM-infM)/0.010) + 1;
  Double_t nbin = ((supM-infM)/0.010);
  
  int H = 600;
  int W = 800;
  TCanvas *c1 = new TCanvas("c1","",50,50,W,H);
  c1->cd();
  c1->SetLeftMargin(0.13);
  c1->SetRightMargin(0.02);
  c1->SetTopMargin(0.09);
  c1->SetBottomMargin(0.13); 
  
  RooPlot* Mframe = Bcm.frame(infM,supM,nbin);
  datatre->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0));
  //model.plotOn(Mframe);
  model.plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));
  model.plotOn(Mframe,Components(sig1),LineColor(kRed),LineWidth(2),Name("Signal")); 
  model.plotOn(Mframe,Components(bkg),LineColor(kGreen),LineWidth(2),Name("bkg")); 
  datatre->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("Data"));
  model.plotOn(Mframe);
  
  Mframe->SetYTitle("Events / 10 MeV"); 
  Mframe->SetLabelSize(0.04,"XY");
  Mframe->SetTitleSize(0.05,"XY");
  Mframe->GetYaxis()->CenterTitle();   
  Mframe->GetXaxis()->CenterTitle();
  Mframe->GetYaxis()->SetNdivisions(505,1);
  Mframe->GetXaxis()->SetNdivisions(505,1);   
  Mframe->GetXaxis()->SetDecimals(1); 
  Mframe->SetTitleOffset(0.9,"X");
  Mframe->SetTitleOffset(1.1,"Y");
  Mframe->SetTitleSize(0.06,"XY");
  Mframe->SetMinimum(1.0); 
  Mframe->Draw();
  
  //TLegend *leg = new TLegend(0.15,0.68,0.35,0.88);
  TLegend *leg = new TLegend(0.18,0.68,0.38,0.88); 
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(Mframe->findObject("Data")," Data","ep"); 
  leg->AddEntry(Mframe->findObject("fittotal")," Fit result","l");
  leg->AddEntry(Mframe->findObject("Signal")," Signal","l");
  leg->AddEntry(Mframe->findObject("bkg"),"Comb. backg.","l");
  leg->Draw();
  
  Double_t Mpsi = mean.getVal()*1000.0;
  Double_t MpsiE = mean.getError()*1000.0;
  Double_t G = sigma1.getVal()*1000.0;
  Double_t GE = sigma1.getError()*1000.0;
  
  TLegend *legpar = new TLegend(0.6,0.68,0.8,0.88);
  legpar->SetTextSize(0.035);
  legpar->SetFillColor(0);
  legpar->SetBorderSize(0);
  legpar->SetFillStyle(0);
  legpar->AddEntry("",Form("Bcm(B_{s}^{0}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
  legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
  legpar->AddEntry("",Form("N_{B_{s}^{0}} = %1.2f #pm %1.2f",nsig.getVal(),nsig.getError()),"");
  legpar->AddEntry("",Form("N_{bkg} = %1.2f #pm %1.2f",nbkg.getVal(),nbkg.getError()),"");
  legpar->Draw();
  
  TLegend *legMass = new TLegend(0.64,0.57,0.83,0.65);
  legMass->SetTextFont(43); 
  legMass->SetTextSize(20);  
  legMass->SetFillColor(0); 
  legMass->SetBorderSize(0);
  legMass->SetFillStyle(0); 
  legMass->SetHeader(Form("%1.1f #leq p_{T}(B_{s}^{0}) < %1.1f GeV ",ptl,pth));
  legMass->Draw(); 
  
  //TLatex *   tex1 = new TLatex(0.92,0.926,"61.2 fb^{-1} (13 TeV)");
  TLatex *   tex1 = new TLatex(0.92,0.926,"ToyMC");
  
  tex1->SetNDC();
  tex1->SetTextAlign(31);
  tex1->SetTextFont(42);
  tex1->SetTextSize(0.05); 
  tex1->SetLineWidth(2);
  
  TLatex *tex2 = new TLatex(0.2,0.926,"CMS");
  tex2->SetNDC();
  tex2->SetTextFont(61);
  tex2->SetTextSize(0.05); 
  tex2->SetLineWidth(2);
  
  TLatex *tex3 = new TLatex(0.29,0.926,"Preliminary");
  tex3->SetNDC();
  tex3->SetTextFont(52);
  tex3->SetTextSize(0.05); 
  tex3->SetLineWidth(2);
  
  tex1->Draw();  
  tex2->Draw();
  tex3->Draw();
  
  c1->Modified();
  gPad->Update();
  gPad->RedrawAxis();
  TLine l;
  l.DrawLine(gPad->GetUxmax(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymin());
  return c1; 
 
}
*/
void ToyMcPar2(Int_t ntotal=1010, Int_t seed=5)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //gStyle->SetErrorX(0);

  Double_t Mmin = 6.05; 
  Double_t Mmax = 6.5;

  RooRealVar Bcm("bcm","bcm",6.05,6.5);
  
  //---- Mass model -------
  Double_t supM = Mmax;
  Double_t infM = Mmin;
  
  RooRealVar a0("a0","a0",-10,10) ;
  RooChebychev bkg("bkg","Background",Bcm,a0) ;
  
  //**** define nominal signal ****
  //gaussians                                                                                                     
  RooRealVar mean("mean","mean of gaussians",6.275,6.05,6.5) ;
  RooRealVar sigma1("sigma1","sigma1 of gaussians",0.01,0.001,0.015) ;
  RooGaussian sig1("sig1","Signal component 1",Bcm,mean,sigma1) ;
    
  //********final PDF ********
  RooRealVar nsig("nsig","number of signal events",0.,1500) ;
  RooRealVar nbkg("nbkg","number of background events",0,1500) ;
  RooAddPdf  model("model","sig+bkg", RooArgList(sig1,bkg), RooArgList(nsig,nbkg)) ;

  //------------ Fit procedure -------------------
  
  //RooRandom::randomGenerator()->SetSeed(3);
  RooRandom::randomGenerator()->SetSeed(seed);
  // esto que sigue solo es para que no imprima en pantalla todo lo que hace del ajuste
  
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(1,false);
  RooMsgService::instance().getStream(1).removeTopic(Integration) ;  
  RooMsgService::instance().getStream(1).removeTopic(Minimization) ;  
  RooMsgService::instance().getStream(1).removeTopic(Fitting) ;  
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration) ;
  RooMsgService::instance().getStream(1).removeTopic(Optimization) ;
  RooMsgService::instance().getStream(1).removeTopic(ObjectHandling) ;
  RooMsgService::instance().getStream(1).removeTopic(Eval) ;
  RooMsgService::instance().Print() ;
  
  Double_t llh,edm,Nsig,NsigE,Nbkg,NbkgE,g1,g1E,mu,muE,conscheb,conschebE, Nspull, Nsdif, Nbpull, Mupull;
  Int_t status,covQual,badNll;
  TFile *file = new TFile("ToyMC_Par2.root","RECREATE");
  TTree *tree = new TTree("tree","tree");
  file->cd();

  tree->Branch("Nspull",&Nspull);
  tree->Branch("Nsdif",&Nsdif);
  tree->Branch("Nbpull",&Nbpull);    
  tree->Branch("Mupull",&Mupull);
  
  tree->Branch("Nsig",&Nsig);
  tree->Branch("NsigE",&NsigE);
  
  tree->Branch("Nbkg",&Nbkg);
  tree->Branch("NbkgE",&NbkgE);

  tree->Branch("mu",&mu);
  tree->Branch("muE",&muE);

  tree->Branch("g1",&g1);
  tree->Branch("g1E",&g1E);
  
  tree->Branch("conscheb",&conscheb);
  tree->Branch("conschebE",&conschebE);

  tree->Branch("edm",&edm);
  tree->Branch("llh",&llh);
  tree->Branch("status",&status);
  tree->Branch("covQual",&covQual);
  tree->Branch("badNll",&badNll);

  
  //****************************************
  //    Input values from DATA fit in pt bins
  //*****************************************
  Double_t nisgf, nsigef, nbkgf, nbkgef, a0f, a0ef;
  Double_t meantxt, meantxte, g1i, g1ie;
 
  ifstream entrada ("output_resulTFit_1_2.txt");
  entrada.is_open();
  if ( !entrada ) 
    {
      cout << "No se pudo abrir el archivo entrada 'resultTFit'" << endl;
      exit( 1 );
    }
  entrada>> nisgf >> nsigef >> nbkgf >> nbkgef >> a0f >> a0ef >> meantxt >> meantxte >> g1i >> g1ie ;
  entrada.close();
  cout << "nsig, nbkg values: " <<  nisgf << " " << nbkgf << endl;
  //return;
  
  Int_t nruns=0;
  
  for(Int_t n=0;n<ntotal;n++)
  {
    //nsig.setConstant(kTRUE);
    nsig.setVal(nisgf);
    nsig.setError(nsigef);
    
    nbkg.setVal(nbkgf);
    nbkg.setError(nbkgef);

    a0.setVal(a0f);
    a0.setError(a0ef);
                
    mean.setVal(meantxt);
    mean.setError(meantxte);
    
    sigma1.setVal(g1i);
    sigma1.setError(g1ie);
    
    RooDataSet *dataToy = model.generate(RooArgSet(Bcm), Extended(kTRUE));
    
    RooFitResult* fitres = model.fitTo(*dataToy,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));

    Nspull = (nsig.getVal()-nisgf)/nsig.getError();
    Nsdif = (nsig.getVal()-nisgf);
    Nbpull = (nbkg.getVal()-nbkgf)/nbkg.getError();   
    Mupull = (mean.getVal()-meantxt)/mean.getError();

    Nsig = nsig.getVal();
    NsigE = nsig.getError();
    
    Nbkg = nbkg.getVal();
    NbkgE = nbkg.getError();

    g1 = sigma1.getVal();
    g1E = sigma1.getError();

    mu = mean.getVal();
    muE = mean.getError();

    conscheb = a0.getVal();
    conschebE = a0.getError();
    
    llh = fitres->minNll();
    status = fitres->status();
    covQual = fitres->covQual();
    edm = fitres->edm();
    badNll = fitres->numInvalidNLL();

    tree->Fill();

    nruns++;
   /* 
    TCanvas* canv_nominal = CreateCanvasNom("canv_nominal", fitres, dataToy, Bcm, supM, infM, model, sig1, bkg, nsig, nbkg, sigma1, mean,ptl,pth); 
    canv_nominal->Print(Form("mass_BsFit_ptbins_ToyMC_%1.0f_%1.0f.png",ptl,pth));
    canv_nominal->Print(Form("mass_BsFit_ptbins_ToyMC_%1.0f_%1.0f.pdf",ptl,pth));
    */
    
    delete dataToy;
    delete fitres;

  }

 tree->Write();
 file->Write("",TObject::kOverwrite);
 
 cout<<" for end: "<<nruns<<endl;

}//End analysis