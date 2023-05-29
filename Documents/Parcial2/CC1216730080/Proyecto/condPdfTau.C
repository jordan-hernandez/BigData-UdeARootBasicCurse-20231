#include "../myTreeClass.C"

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
using namespace std;

void condPdfTau() {
    Float_t Mmax = 6.50;
    Float_t Mmin = 6.04;
    Float_t tmin = 0.3;
    Float_t tmax = 2.6;
    Float_t terrmin = 0.0001;
    Float_t terrmax = 0.5;

    TChain *ch = new TChain("bctree","");
    ch->Add("../Bc_data.root/tree");

    TTree *tree = (TTree*) ch;
    myTreeClass myt(tree);
    Long64_t nentries = myt.fChain->GetEntries();
    cout<<" Entries : "<<nentries<<endl;

    //------------------------------
    RooRealVar M("M","M(J/#psi #pi^{+}) (GeV)",Mmin,Mmax);
    RooRealVar Tau("Tau","#tau(J/#psi #pi^{+}) (ps)",tmin,tmax);
    RooRealVar TauErr("TauError","#Delta#tau(J/#psi #pi^{+}) (ps)",terrmin,terrmax);
    RooDataSet data("data","Filtered data Lifetime",RooArgSet(Tau));
    RooDataSet dataErr("dataErr","Filtered data LifetimeError",RooArgSet(TauErr));

    Int_t nTen = nentries/10;
    Int_t k=0;
    Int_t nbytes = 0, nb = 0;
    for(Long64_t jentry=0; jentry<nentries;jentry++){
        Long64_t ientry = myt.LoadTree(jentry);
        if (ientry < 0) break;
        nb = myt.fChain->GetEntry(jentry);   nbytes += nb;
        if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
        if(jentry==nentries-1) cout<<endl;

        //Mass windows cuts
        if(myt.Bcmass<=Mmin || myt.Bcmass>=Mmax) continue;

        //Lifetime windows cuts
        if(myt.Tau<=tmin || myt.Tau>=tmax) continue;

        //Lifetime error window cuts	
        if(myt.TauError<=terrmin || myt.TauError>=terrmax) continue;

        //Filtro para lifetime.
        if((myt.Tau/myt.TauError)<5.0)continue;

        Tau = myt.Tau;
        TauErr = myt.TauError;    
        data.add(RooArgSet(Tau));
        dataErr.add(RooArgSet(TauErr));

    }

    cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
    data.Print("v");

    RooRealVar t("t", "t", tmin, tmax);
    //RooRealVar TauE("Tauerr", "per-event error on tau", terrmin, terrmax);

    //---------- Gaussian resolution model scaled by per-event error ------------
    RooRealVar bias("bias", "bias", 0, -10., 10.);
    RooRealVar sigma("sigma", "per-event error scale factor", 1, 0.0, 10);
    RooGaussModel gm("gm1", "gauss model scaled bt per-event error", Tau, bias, sigma, TauErr);

    //---------- Decay (x) gauss(tau|tauErr) --------------
    RooRealVar tau("tau", "tau", 0.507);
    RooDecay decay_gm("decay_gm", "decay", Tau, tau, gm, RooDecay::SingleSided);

    //---------- Fit conditional decay --------------
    decay_gm.fitTo(data, ConditionalObservables(TauErr));

    // P l o t   c o n d i t i o n a l   d e c a y _ d m ( d t | d t e r r )
    // ---------------------------------------------------------------------
 
    // Make two-dimensional plot of conditional pdf in (dt,dterr)
    TH1 *hh_decay = decay_gm.createHistogram("hh_decay", Tau, Binning(50), YVar(TauErr, Binning(50)));
    hh_decay->SetLineColor(kBlue);
    
    // Plot decay_gm(dt|dterr) at various values of dterr
    RooPlot *frame = Tau.frame(Title("Slices of decay(#tau|#Delta#tau_{err}) at various #Delta#tau_{err}"));
    for (Int_t ibin = 0; ibin < 100; ibin += 20) {
        TauErr.setBin(ibin);
        decay_gm.plotOn(frame, Normalization(5.));
    }
    
    // Make projection of data an dt
    RooPlot *frame2 = Tau.frame(Title("Projection of decay(#tau|#Delta#tau_{err}) on #tau"));
    data.plotOn(frame2);
    
    // Make projection of decay(dt|dterr) on dt.
    //
    // Instead of integrating out dterr, make a weighted average of curves
    // at values dterr_i as given in the external dataset.
    // (The true argument bins the data before projection to speed up the process)
    decay_gm.plotOn(frame2, ProjWData(dataErr, true));
    
    // Draw all frames on canvas
    TCanvas *c = new TCanvas("rf306_condpereventerrors", "rf306_condperventerrors", 1200, 400);
    c->Divide(3);
    c->cd(1);
    gPad->SetLeftMargin(0.20);
    hh_decay->GetZaxis()->SetTitleOffset(2.5);
    hh_decay->Draw("surf");
    c->cd(2);
    gPad->SetLeftMargin(0.15);
    frame->GetYaxis()->SetTitleOffset(1.6);
    frame->Draw();
    c->cd(3);
    gPad->SetLeftMargin(0.15);
    frame2->GetYaxis()->SetTitleOffset(1.6);
    frame2->Draw();
    c->Print("conditional_with_per_event_errors.png");
}