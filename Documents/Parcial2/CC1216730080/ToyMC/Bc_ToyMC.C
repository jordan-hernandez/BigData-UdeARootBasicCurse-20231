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

TCanvas* CreateCanvasNom(TString cname, RooFitResult* result, RooDataSet *data, RooRealVar M, Double_t supM, Double_t infM,  RooAddPdf MassModel, RooAddPdf sumgau, RooExponential bkg1, RooRealVar Ns, RooRealVar Nb, RooRealVar width, RooRealVar width2, RooRealVar fs, RooRealVar mean, Double_t ptl, Double_t pth)
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
  
  RooPlot* Mframe = M.frame(infM,supM,nbin);
  data->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0));
  //MassModel.plotOn(Mframe);
  MassModel.plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));
  MassModel.plotOn(Mframe,Components(sumgau),LineColor(kRed),LineWidth(2),Name("Signal")); 
  MassModel.plotOn(Mframe,Components(bkg1),LineColor(kGreen),LineWidth(2),Name("bkg")); 
  data->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("Data"));
  MassModel.plotOn(Mframe);
  
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
  //Double_t G = width.getVal()*1000.0;
  //Double_t GE = width.getError()*1000.0;
  Double_t G = sqrt( fs.getVal()*width.getVal()*width.getVal() + (1-fs.getVal())*width2.getVal()*width2.getVal() )*1000.0;
  Double_t GE = (1/G)*sqrt( (fs.getVal()*fs.getVal())*(width.getVal()*width.getVal())*(width.getError()*width.getError()) + ((1-fs.getVal())*(1-fs.getVal()))*(width2.getVal()*width2.getVal())*(width2.getError()*width2.getError()) )*1000.0*1000.0;
  
  TLegend *legpar = new TLegend(0.6,0.68,0.8,0.88);
  legpar->SetTextSize(0.035);
  legpar->SetFillColor(0);
  legpar->SetBorderSize(0);
  legpar->SetFillStyle(0);
  legpar->AddEntry("",Form("M(B_{s}^{0}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
  legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
  legpar->AddEntry("",Form("N_{B_{s}^{0}} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
  legpar->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
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

void Bc_ToyMC(Int_t ntotal=1000, Int_t seed=5, Double_t ptl=12, Double_t pth=70)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //gStyle->SetErrorX(0);

  Double_t Mmin = 6.04; 
  Double_t Mmax = 6.50;

  RooRealVar M("M","M(J/#psi #pi^{+}) (GeV)",Mmin,Mmax);
  //RooDataSet data("data","data",RooArgSet(M));
 
  //RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-10) ;
  //RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-10) ;
  
  //---- Mass model -------
  Double_t supM = Mmax;
  Double_t infM = Mmin;
  
  RooRealVar c("c","c",-20.0,0.0);
  RooExponential bkg1("bkg1","Exp. Background",M,c);
  
  //**** define nominal signal ****
  //gaussians                                                                                                                                                   
  //RooRealVar mean("mean"," Mass mean",3.686097);
  RooRealVar mean("mean"," Mass mean",6.275,6.25,6.3,"GeV");
  RooRealVar width("width"," Mass width",0.10,0.01,0.15,"GeV");
  RooGaussian Sig("Sig"," Signal PDF",M,mean,width);
    
  //********final PDF ********
  RooRealVar Ns("Ns","Ns",0.,20000);
  RooRealVar Nb("Nb","Nb",0.,20000); 

  //RooAddPdf MassModel("MassModel","MassModel",RooArgList(Sig,bkg1),RooArgList(Ns,Nb));
  RooAddPdf MassModel("MassModel","MassModel",RooArgList(Sig,bkg1),RooArgList(Ns,Nb));

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
  
  Double_t llh,edm,Nsig,NsigE,Nbkg,NbkgE,g1,g1E,mu,muE,g2,g2E,exp,expE, FS, FSE, Nspull, Nsdif, Nbpull, Mupull;
  Int_t status,covQual,badNll;
  //TFile *file = new TFile("ToyMC_Bs.root","RECREATE");
  //TFile *file = new TFile(Form("ToyMC_Bs_ptbins_%1i_%1.0f_%1.0f.root",seed,ptl,pth),"RECREATE");
  TFile *file = TFile::Open(Form("ToyMC_Bs_ptbins_%1i_%1.0f_%1.0f.root",seed,ptl,pth),"RECREATE");
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

  tree->Branch("g2",&g2);
  tree->Branch("g2E",&g2E);

  tree->Branch("FS",&FS);
  tree->Branch("FSE",&FSE);
  
  tree->Branch("exp",&exp);
  tree->Branch("expE",&expE);

  tree->Branch("edm",&edm);
  tree->Branch("llh",&llh);
  tree->Branch("status",&status);
  tree->Branch("covQual",&covQual);
  tree->Branch("badNll",&badNll);

  //****************************************
  //    Input values from DATA fit in pt bins
  //*****************************************
  Double_t nsi, nsie, nbi, nbie, ci, cie, fsi, fsie;
  Double_t mui, muie, g1i, g1ie, g2i, g2ie;
 
  ifstream entrada (Form("../output_TotalFit_BsFit_nominal_%1.0f_%1.0f.txt",ptl,pth));
  //entrada.is_open();
  if ( !entrada ) 
    {
      cout << "No se pudo abrir el archivo entrada fit pr bins" << endl;
      exit( 1 );
    }
  entrada>> nsi >> nsie >> nbi >> nbie >> ci >> cie >> mui >> muie >> g1i >> g1ie;
  entrada.close();
  cout << "Ns, Nb values: " <<  nsi << " " << nbi << endl;
  //return;
  
  Int_t nruns=0;
  
  for(Int_t n=0;n<ntotal;n++)
  {
    //Ns.setConstant(kTRUE);
    Ns.setVal(nsi);
    Ns.setError(nsie);
    
    Nb.setVal(nbi);
    Nb.setError(nbie);

    c.setVal(ci);
    c.setError(cie);
        
    mean.setVal(mui);
    mean.setError(muie);
    
    width.setVal(g1i);
    width.setError(g1ie);

    //RooRandom::randomGenerator()->SetSeed(n);
    RooDataSet *dataToy = MassModel.generate(RooArgSet(M), Extended(kTRUE));
    //RooDataSet *dataToy = MassModel.generate(RooArgSet(M),Ns.getVal() +  Nb.getVal() );
    //RooDataSet dataToy = MassModel.generate(RooArgSet(M), Ns.getVal() + Nb.getVal() );
    
    RooFitResult* fitres = MassModel.fitTo(*dataToy,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));
    //RooFitResult* fitres = MassModel.fitTo(data,Extended(),Minos(kTRUE),Save(kTRUE), NumCPU(4));
    //dataToy->Print("v");
    //fitres->Print("v");

    Nspull = (Ns.getVal()-nsi)/Ns.getError();
    Nsdif = (Ns.getVal()-nsi);
    Nbpull = (Nb.getVal()-nbi)/Nb.getError();   
    Mupull = (mean.getVal()-mui)/mean.getError();

    Nsig = Ns.getVal();
    NsigE = Ns.getError();
    
    Nbkg = Nb.getVal();
    NbkgE = Nb.getError();

    g1 = width.getVal();
    g1E = width.getError();

    mu = mean.getVal();
    muE = mean.getError();

    exp = c.getVal();
    expE = c.getError();
    
    llh = fitres->minNll();
    status = fitres->status();
    covQual = fitres->covQual();
    edm = fitres->edm();
    badNll = fitres->numInvalidNLL();

    tree->Fill();

    nruns++;
    /*
    TCanvas* canv_nominal = CreateCanvasNom("canv_nominal", fitres, dataToy, M, supM, infM, MassModel, sumgau, bkg1, Ns, Nb, width, width2, fs, mean,ptl,pth); 
    canv_nominal->Print(Form("plots/mass_BsFit_ptbins_ToyMC_%1.0f_%1.0f.png",ptl,pth));
    canv_nominal->Print(Form("plots/mass_BsFit_ptbins_ToyMC_%1.0f_%1.0f.pdf",ptl,pth));
    */
    
    delete dataToy;
    delete fitres;

  }

 tree->Write();
 tree->MakeClass("classtoy");
 //file->Write("",TObject::kOverwrite);
 
 cout<<" for end: "<<nruns<<endl;

}//End analysis