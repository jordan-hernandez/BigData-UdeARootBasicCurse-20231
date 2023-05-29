#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooDecay.h"
#include "RooLandau.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
using namespace RooFit;
 
void bc_2d()
{
    Float_t Bcmass, tau , tauerror;
    
    TChain *ch = new TChain("bstree","");
    ch->Add("reducetree_Bc2d_AOD_best.root/mytree");

    TTree *tree = (TTree*) ch;
    //classreduce t(tree);
    Long64_t nentries = tree->GetEntries();
    cout<<" Entries : "<<nentries<<endl;

    tree->SetBranchAddress ("Bcmass", &Bcmass);
    tree->SetBranchAddress ("tau", &tau);
    tree->SetBranchAddress ("tauerror", &tauerror);
   // B - p h y s i c s   p d f   w i t h   p e r - e v e n t  G a u s s i a n   r e s o l u t i o n
   // ----------------------------------------------------------------------------------------------
    Float_t Mmin = 6.05; 
    Float_t Mmax = 6.5;
    Float_t Tmin = 0.3; 
    Float_t Tmax = 2.6; 
    Float_t Terrmin = 0.0001; 
    Float_t Terrmax = 0.5;

   // Observables
    RooRealVar M("M", "M", Mmin, Mmax);
    RooRealVar Tau("Tau", "Tau", Tmin, Tmax);

    RooDataSet data_M("data_M","data_M",RooArgSet(M));
    RooDataSet data_tau("data_tau","data_tau",RooArgSet(Tau));

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
        if(tau<=Tmin || tau>=Tmax) continue;	
        if(tauerror<=Terrmin || tauerror>=Terrmax) continue;	
        if((tau/tauerror)<5.0) continue;
        M = Bcmass;
        Tau=tau;

        data_M.add(RooArgSet(M));
        data_tau.add(RooArgSet(Tau));

        }

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
   RooRealVar c0("c","c",-10.0,10.0);
   RooChebychev  bkg1("bkg1","Exp. Background",M,c0);

   RooRealVar Ns("Ns","Ns",0.,3000);
   RooRealVar Nb("Nb","Nb",0.,3000);   

   RooAddPdf MassModel("MassModel","MassModel",RooArgList(sumgau,bkg1),RooArgList(Ns,Nb));
   RooFitResult* fitres = MassModel.fitTo(data_M,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));

   //Lifetime
   RooRealVar cte("c2", "c2", -1.85, -10,10);
   RooExponential lifetime("lifetime", "lifetime", Tau, cte);
   
   // Conditional
   RooProdPdf massandlifetime("massandlifetime", "massandlifetime", MassModel, Conditional(lifetime, Tau));

   TH1 *hh = massandlifetime.createHistogram("hh", M, Binning(50), YVar(Tau, Binning(50)));
   hh->SetLineColor(kBlue);
 
   Double_t nbin = ((Mmax-Mmin)/0.005)+1;
   RooPlot* frame = M.frame(Mmin,Mmax,nbin);
   data_M.plotOn(frame,DataError(RooAbsData::SumW2),MarkerSize(.65),XErrorSize(0));
   //MassModel.plotOn(frame,LineColor(kRed),LineWidth(2),Name("Fit"));
   massandlifetime.plotOn(frame,LineColor(kBlue),LineWidth(2),Name("Fit"));
 
   RooPlot *frame2 = Tau.frame(Title("Lifetime model"));
   data_tau.plotOn(frame2);
   lifetime.plotOn(frame2);
 
   // Draw all frames on canvas
   TCanvas *c = new TCanvas("massandlifetime", "massandlifetime", 1200, 400);
   c->Divide(3);
   c->cd(1);
   gPad->SetLeftMargin(0.20);
   hh->GetZaxis()->SetTitleOffset(2.5);
   hh->Draw("surf");
   c->cd(2);
   gPad->SetLeftMargin(0.15);
   frame->GetYaxis()->SetTitleOffset(1.6);
   frame->Draw();
   c->cd(3);
   gPad->SetLeftMargin(0.15);
   frame2->GetYaxis()->SetTitleOffset(1.6);
   frame2->Draw();
}