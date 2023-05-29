#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>

#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooArgusBG.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"

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

using namespace RooFit;

void save_resultTfit(ofstream& salida, RooRealVar Ns, RooRealVar Nb, RooRealVar c, RooRealVar fs, RooRealVar mean,  RooRealVar width, RooRealVar width2)
{
  salida <<  Ns.getVal() << " " << Ns.getError() << " " << Nb.getVal() << " " << Nb.getError() <<" " <<  c.getVal() << " " << c.getError() << " "  
	 << fs.getVal() << " " << fs.getError() << " " << mean.getVal() << " " << mean.getError() << " "
	 <<  width.getVal() << " " << width.getError() << " " <<  width2.getVal() << " " << width2.getError()
    ;
  cout << " el archivo se escribio bien" << endl; 
  return;
}


void save_result(ofstream& salida, RooRealVar Ns, RooRealVar Nb)
{
  salida <<  Ns.getVal() << " " << Ns.getError() << " " << Nb.getVal() << " " << Nb.getError();
  cout << " el archivo se escribio bien" << endl; 
  return;
}


int bc2dFit(Double_t ptl=12.0, Double_t pth=70.0)
{
Float_t Bcmass, tau , tauerror;
    
TChain *ch = new TChain("bstree","");
ch->Add("reducetree_Bc2d_AOD_best.root/mytree");

TTree *tree = (TTree*) ch;
//classreduce t(tree);
Long64_t nentries = tree->GetEntries();
cout<<" Entries : "<<nentries<<endl;
/*
    TFile f ("reducetree_Bc2d_AOD_best.root", "READ");
    TTree *tree= (TTree*)f.Get("mytree");
    */
tree->SetBranchAddress ("Bcmass", &Bcmass);
tree->SetBranchAddress ("tau", &tau);
tree->SetBranchAddress ("tauerror", &tauerror);



Double_t Mmin = 6.05; 
Double_t Mmax = 6.5; 

RooRealVar M("M","",Mmin,Mmax);
RooDataSet data("data","data",RooArgSet(M));

Int_t nTen = nentries/10;
Int_t k=0;
Int_t nbytes = 0, nb = 0;

for(Long64_t jentry=0; jentry<nentries;jentry++)
    {
      Long64_t ientry = tree->LoadTree(jentry);
      if (ientry < 0) break;
      nb = tree->GetEntry(jentry);   nbytes += nb;
      if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
      if(jentry==nentries-1) cout<<endl;

      //Mass windows cuts
      if(Bcmass<=Mmin || Bcmass>=Mmax) continue;	
      if((tau/tauerror)<5.0) continue;
      M=Bcmass;    
      data.add(RooArgSet(M));
    }

//---- Mass model -------
 Double_t supM = Mmax;
 Double_t infM = Mmin;

//FIT
// Signal
RooRealVar mean("mean"," Mass mean",6.2751,6.25,6.3,"GeV");
RooRealVar width("width"," Mass width",0.015,0.01,0.04,"GeV");
RooGaussian Sig("Sig"," Signal PDF",M,mean,width);

RooRealVar width2("width"," Mass width",0.020,0.005,0.05,"GeV");
RooGaussian Sig2("Sig2"," Signal 2 PDF",M,mean,width2);

RooRealVar fs("fs","fs",0.8,0.,1.);
RooAddPdf sumgau("sumgau","sumgau",RooArgList(Sig,Sig2),RooArgList(fs));

// Background
RooRealVar c("c","c",-10.0,10.0);
RooChebychev  bkg1("bkg1","Exp. Background",M,c);

RooRealVar Ns("Ns","Ns",0.,2000);
RooRealVar Nb("Nb","Nb",0.,2000);   

RooAddPdf MassModel("MassModel","MassModel",RooArgList(sumgau,bkg1),RooArgList(Ns,Nb));
RooFitResult* fitres = MassModel.fitTo(data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));

data.Print("v"); 
fitres->Print("v");
Double_t nbin = ((Mmax-Mmin)/0.005)+1;

TCanvas *c1 = new TCanvas("basic", "basic",50,50, 800, 600);
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
pad2->SetFillColor(0);
pad2->SetGridx(0);
pad2->SetGridy(0);

pad1->Draw();
pad2->Draw();
pad1->cd(); 

//RooPlot* framemm = M.frame(Title("")) ;
RooPlot* frame = M.frame(Mmin,Mmax,nbin) ;
data.plotOn(frame,DataError(RooAbsData::SumW2),MarkerSize(.65),XErrorSize(0));
MassModel.plotOn(frame,LineColor(kRed),LineWidth(2),Name("Fit")); 

RooHist* hpullm2 = frame->pullHist() ;


//MassModel.plotOn(framemm,Components(sumgau),LineColor(kBlue),LineWidth(2),Name("Signal")); 
MassModel.plotOn(frame,Components(bkg1),LineColor(kGreen),LineWidth(2),Name("bkg")); 
MassModel.plotOn(frame,Components(sumgau),LineColor(kBlue),LineWidth(2),Name("Signal")); 
//framemm->Draw();

frame->SetTitle("Fit to the B_{c} mass");
frame->SetTitleSize(0.8,"");
frame->SetYTitle("Events");        
frame->SetLabelSize(0.06,"XY");
frame->SetTitleSize(0.19,"XY");
frame->GetYaxis()->CenterTitle();   
frame->GetXaxis()->CenterTitle();
frame->GetYaxis()->SetNdivisions(505,1);
frame->GetXaxis()->SetNdivisions(505,1);   
frame->GetXaxis()->SetDecimals(1);
frame->GetXaxis()->SetTickLength(0.06); 
frame->SetTitleSize(0.11,"X");
frame->SetTitleSize(0.1,"Y");  
frame->SetTitleOffset(0.9,"X");
frame->SetTitleOffset(0.42,"Y");
frame->SetMinimum(-2.0);   
frame->Draw();

//TLegend *leg = new TLegend(0.15,0.68,0.35,0.88);
TLegend *leg = new TLegend(0.18,0.58,0.38,0.88); 
leg->SetTextSize(0.06);
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->AddEntry(frame->findObject("Data")," Data","ep"); 
leg->AddEntry(frame->findObject("Fit")," Fit result","l");
leg->AddEntry(frame->findObject("Signal")," Signal","l");
leg->AddEntry(frame->findObject("bkg"),"Comb. backg.","l");
leg->Draw();

Double_t Mpsi = mean.getVal()*1000.0;
Double_t MpsiE = mean.getError()*1000.0;
Double_t G = sqrt( fs.getVal()*width.getVal()*width.getVal() + (1-fs.getVal())*width2.getVal()*width2.getVal() )*1000.0;
Double_t GE = (1/G)*sqrt( (fs.getVal()*fs.getVal())*(width.getVal()*width.getVal())*(width.getError()*width.getError()) + ((1-fs.getVal())*(1-fs.getVal()))*(width2.getVal()*width2.getVal())*(width2.getError()*width2.getError()) )*1000.0*1000.0;

TLegend *legpar = new TLegend(0.6,0.58,0.8,0.88);
legpar->SetTextSize(0.06);
legpar->SetTextFont(42);
legpar->SetFillColor(0);
legpar->SetBorderSize(0);
legpar->SetFillStyle(0);
legpar->AddEntry("",Form("M(B_{c}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
legpar->AddEntry("",Form("N_{B_{c}} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
legpar->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
legpar->Draw();

pad2->cd();
// RooPlot* Mframe = M.frame(Mmin,Mmax,nbin) ;

RooPlot* framem2 = M.frame() ;


framem2->addPlotable(hpullm2,"P") ;
framem2->SetTitle("");
framem2->SetYTitle(" (N^{fit}-N^{true})/#sigma_{fit}");
framem2->SetXTitle("Mass (Gev)");
framem2->SetLabelSize(0.1,"XY");
framem2->SetTitleSize(0.11,"X");
framem2->SetTitleSize(0.135,"Y");  
framem2->GetYaxis()->CenterTitle();   
framem2->GetXaxis()->CenterTitle();
framem2->GetYaxis()->SetNdivisions(505,1);
framem2->GetXaxis()->SetNdivisions(505,1);
framem2->GetXaxis()->SetTickLength(0.07);   
framem2->SetTitleOffset(0.9,"X");
framem2->SetTitleOffset(0.31,"Y");
framem2->SetMaximum(4.9);
framem2->SetMinimum(-4.9);

framem2->Draw();

TLine *line1 = new TLine(Mmin,0.0,Mmax,0.0);
line1->SetLineColor(1);
//line1->SetLineStyle(2);
line1->SetLineWidth(1);
line1->Draw();


TLatex *tex2 = new TLatex(0.15,0.95,"CMS");
tex2->SetNDC();
tex2->SetTextFont(61);
tex2->SetTextSize(0.04); 
tex2->SetLineWidth(2);
tex2->Draw();

c1->Modified();
c1->Print("fit.png");

 //parameters output file
 ofstream salida_TotalFit(Form("plots_ptandeta/output_TotalFit_BcFit_nominal_%1.0f_%1.0f.txt",ptl,pth));
 salida_TotalFit.is_open();
 save_resultTfit(salida_TotalFit, Ns, Nb, c, fs, mean, width, width2);

/*
 //parameters output file
 ofstream *salida_nominal = save_result(Form("plots/output_BsFit_nominal_%1i_%1i.txt",era,vq),mean, Gt, GtE, Ns, Nb);
 ofstream salida_nominal(Form("plots_ptandeta/output_BsFit_nominal_%1.0f_%1.0f.txt",ptl,pth));
 salida_nominal.is_open();
 save_result(salida_nominal, Ns, Nb);// OJO esta es la funciona a la que se le pasa el archivo por referencia
 
 //made canvas
 TCanvas* canv_nominal = CreateCanvasNom("canv_nominal", fitres, data, M, supM, infM, MassModel, sumgau, bkg1, Ns, Nb, width, width2, fs, mean,ptl,pth);
 canv_nominal->Print(Form("plots_ptandeta/mass_BsFit_ptbins_%1.0f_%1.0f_NOpull.png",ptl,pth));
 canv_nominal->Print(Form("plots_ptandeta/mass_BsFit_ptbins_%1.0f_%1.0f_NOpull.pdf",ptl,pth));
 TCanvas* canv_nominalpull = CreateCanvasNomPull("canv_nominalpull", fitres, data, M, supM, infM, MassModel, sumgau, bkg1, Ns, Nb, width, width2, fs, mean,ptl,pth);  
 canv_nominalpull->Print(Form("plots_ptandeta/mass_BsFitpull_ptbins_%1.0f_%1.0f.png",ptl,pth));
 canv_nominalpull->Print(Form("plots_ptandeta/mass_BsFitpull_ptbins_%1.0f_%1.0f.pdf",ptl,pth));
*/
return 0;
}