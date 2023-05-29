     
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

void fit2d(){
    Float_t mass, time, error;

    TFile f("masayT.root", "READ");
    TTree *tree = (TTree*)f.Get("tree");

    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("time", &time);
    tree->SetBranchAddress("error", &error);


    //Masa del Bc: 6.2751GeV
    Double_t Mmin = 6.05; 
    Double_t Mmax = 6.5; 

    Double_t taumin = 0.3; 
    Double_t taumax = 2.6; 

    Double_t errmin = 0.0001; 
    Double_t errmax = 0.4; 

    RooRealVar M("M","M(B_{c}^{+}) (GeV)",Mmin,Mmax);
    RooRealVar tau("tau","Lifetime (ps)",taumin,taumax);
    // RooRealVar err("err","Lifetime error",errmin,errmax);

    RooDataSet dataMt("dataMt","dataMt",RooArgSet(M,tau));
    // RooDataSet dataerr("dataerr","dataerr",RooArgSet(err));


    Long64_t nentries = tree->GetEntries();
    // cout<<" Entries : "<<nentries<<endl;

    for (int evt=0; evt < nentries; evt++) 
    {
        tree->GetEvent (evt);

        if (error==0 || isnan(error)) continue;
        if (time<taumin || time>taumax) continue;
        if (error<errmin || error>errmax) continue;

        if((time/error)<5.0) continue;

        M=mass;
        tau=time;   
        dataMt.add(RooArgSet(M,tau));
        // err=error;
        // dataerr.add(RooArgSet(err));
        
    }
    // dataMt.Print("v");


    RooRealVar mean("mean"," Mass mean",6.2751,6.25,6.3,"GeV");
    RooRealVar width("width"," Mass width",0.015,0.01,0.04,"GeV");
    RooGaussian gauss("Sig"," Signal PDF",M,mean,width);


    RooRealVar a("a","a",-10,10.0);
    RooChebychev bkg("bkg","Background",M,a);

    RooRealVar Ns("Ns","Ns",0.,2000);
    RooRealVar Nb("Nb","Nb",0.,2000); 

    RooAddPdf MassModel("MassModel","MassModel",RooArgList(gauss,bkg),RooArgList(Ns,Nb));


    RooRealVar tmed("tmed", "tmed", -1, -3,2);
    RooExponential decay("decay","decay",tau,tmed);

    RooProdPdf pdf2d("pdf2d", "pdf2d", MassModel, Conditional(decay, tau));

    RooFitResult* fitres = pdf2d.fitTo(dataMt,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4), 
                                        ConditionalObservables(decay,tau));

    fitres->Print("v");

    TH1 *hh_decay = pdf2d.createHistogram("hh_decay", M, Binning(50), YVar(tau, Binning(50)));
    // hh_decay->SetLineColor(kBlue);

    //Gráfica

    RooPlot *frame1 = tau.frame(Title("Lifetime fit"));
    dataMt.plotOn(frame1);
    pdf2d.plotOn(frame1);

    RooPlot *frame2 = M.frame(Title("Mass fit"));
    dataMt.plotOn(frame2);
    pdf2d.plotOn(frame2);



    TCanvas *c = new TCanvas("c", "c", 1200, 400);
    c->Divide(3);

    c->cd(1);
    gPad->SetLeftMargin(0.20);
    hh_decay->SetTitleOffset(2.5,"Z");
    hh_decay->SetTitleOffset(2,"Y");
    hh_decay->SetTitleOffset(2,"X");
    hh_decay->SetTitle("Mass-Lifetime relation");
    hh_decay->SetLabelSize(0.03,"X");
    hh_decay->GetZaxis()->SetNdivisions(6,1);
    hh_decay->Draw("surf2");   

    c->cd(2);
    gPad->SetLeftMargin(0.15);
    
    frame1->GetYaxis()->SetTitleOffset(1.6);
    frame1->Draw();

    TLegend *leg1 = new TLegend(0.48,0.75,0.8,0.88); 
    leg1->SetTextSize(0.04);
    leg1->SetFillColor(0);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->AddEntry("","Halftime paramater:","");
    leg1->AddEntry("",Form("#tau = %1.3f #pm %1.3f ps",-1/tmed.getVal(),1/tmed.getError()),"");
    leg1->Draw();

    c->cd(3);
    gPad->SetLeftMargin(0.15);
    frame2->GetYaxis()->SetTitleOffset(1.6);
    frame2->Draw();

    TLegend *leg2 = new TLegend(0.52,0.75,0.8,0.88); 
    leg2->SetTextSize(0.03);
    leg2->SetFillColor(0);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry("",Form("M(B_{c}^{+})= %1.1f #pm %1.1f MeV",mean.getVal()*1000,mean.getError()*1000),"");
    leg2->AddEntry("",Form("#sigma = %1.1f #pm %1.1f MeV",width.getVal()*1000,width.getError()*1000),"");
    leg2->AddEntry("",Form("N_{B_{c}^{+}} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
    leg2->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
    leg2->Draw();

    c->Draw();
    c->Print("fit2d.png");
}