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


// -------------------------------------------
// PLOT DATA & FIT
// -------------------------------------------
TCanvas* CreateCanvas(TString cname, RooFitResult* result, RooDataSet *data, 
RooRealVar M, Double_t Mmax, Double_t Mmin,  RooAddPdf MassModel, RooAddPdf sumGauss, 
RooChebychev bkg, RooRealVar Ns, RooRealVar Nb, RooRealVar sigma1, RooRealVar sigma2, 
RooRealVar fraction, RooRealVar meanGauss)
{


    int H = 800;
    int W = 1000;

    TCanvas *c1 = new TCanvas();

    c1->cd();
    c1->SetLeftMargin(0.13);
    c1->SetRightMargin(0.02);
    c1->SetTopMargin(0.09);
    c1->SetBottomMargin(0.13); 

    Double_t supM = Mmax;
    Double_t infM = Mmin;
    Double_t nbin = ((supM - infM) / 0.010);

    RooPlot* Mframe = M.frame(infM,supM,nbin);

    // Data
    data->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0));

    // Model
    MassModel.plotOn(Mframe);
    MassModel.plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));
    MassModel.plotOn(Mframe, Components(sumGauss),LineColor(6),LineWidth(2),LineStyle(kDashed),Name("signal")); 
    
    // Background 
    MassModel.plotOn(Mframe,Components(bkg),LineColor(kGreen),LineWidth(2),Name("bkg"));
    
    data->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("data"));
    MassModel.plotOn(Mframe);

    Mframe->SetYTitle("Events / 10 MeV");
    Mframe->SetLabelSize(0.04, "XY");
    Mframe->SetTitleSize(0.05, "XY");
    Mframe->GetYaxis()->CenterTitle();
    Mframe->GetXaxis()->CenterTitle();
    Mframe->GetYaxis()->SetNdivisions(506, 1);
    Mframe->GetXaxis()->SetNdivisions(510, 1);
    Mframe->GetXaxis()->SetDecimals(1);
    Mframe->SetTitleOffset(0.9, "X");
    Mframe->SetTitleOffset(0.8, "Y");
    Mframe->SetTitleSize(0.06, "XY");
    Mframe->SetMinimum(0.01);
    Mframe->SetMaximum(600);
    Mframe->Draw();
    gStyle->SetOptTitle(0/1);


    // Legend 1
    TLegend *legend1 = new TLegend(0.18,0.65,0.38,0.88); 
    legend1->SetTextSize(0.04);
    legend1->SetFillColor(0);
    legend1->SetBorderSize(0);
    legend1->SetFillStyle(0);
    legend1->AddEntry(Mframe->findObject("data"), "Data", "ep"); 
    legend1->AddEntry(Mframe->findObject("fittotal"), "Fit", "l");
    legend1->AddEntry(Mframe->findObject("signal"), "B_{c}^{+} Signal", "l");
    legend1->AddEntry(Mframe->findObject("bkg"), "Comb. background.", "l");
    legend1->Draw();


    // Legend 2
    Double_t Mpsi = meanGauss.getVal()*1000.0;
    Double_t MpsiE = meanGauss.getError()*1000.0;
    Double_t G = sqrt( fraction.getVal()*sigma1.getVal()*sigma1.getVal() + (1-fraction.getVal())*sigma2.getVal()*sigma2.getVal() )*1000.0;
    Double_t GE = (1/G)*sqrt( (fraction.getVal()*fraction.getVal())*(sigma1.getVal()*sigma1.getVal())*(sigma1.getError()*sigma1.getError()) + ((1-fraction.getVal())*(1-fraction.getVal()))*(sigma2.getVal()*sigma2.getVal())*(sigma2.getError()*sigma2.getError()) )*1000.0*1000.0;
    
    TLegend *legend2 = new TLegend(0.6,0.6,0.8,0.88);
    legend2->SetTextSize(0.04);
    legend2->SetFillColor(0);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->AddEntry("",Form("M(B_{c}^{+}) = %1.2f #pm %1.2f MeV", Mpsi, MpsiE),"");
    legend2->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV", G, GE),"");
    legend2->AddEntry("",Form("N_{B_{c}^{+}} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
    legend2->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
    legend2->Draw();

    return c1;


}


// -------------------------------------------
// TOY MONTECARLO
// -------------------------------------------
void toyMC()
{

    Int_t seed=5;
    Int_t n_total=500;


    // -------------------------------------------
    // CREATE MODEL 
    // -------------------------------------------

    Double_t Mmin = 6.05;
    Double_t Mmax = 6.5;
    RooRealVar M("M", " M(J/#psi #pi^{+}) [GeV]", Mmin, Mmax); // masa del Bc 
  
    // Create two Gaussian PDFs for the signal 
    RooRealVar meanGauss("#muB"," Mass mean Bc", 6.275, 6.225, 6.300,"GeV");
    RooRealVar sigma1("sigma1"," Mass width", 0.010, 0.001, 0.015, "GeV");
    RooRealVar sigma2("sigma2"," Mass width", 0.020, 0.015, 0.05, "GeV");
    
    RooGaussian signal1("signal1", "Signal 1", M, meanGauss, sigma1);
    RooGaussian signal2("signal2", "Signal 2", M, meanGauss, sigma2);
    
    // Build Chebychev polynomial pdf for background 
    RooRealVar a0("a0", "a0", -10., 10.);
    RooRealVar a1("a1", "a1", -10., 10.);
    RooChebychev bkg("bkg","Background", M, RooArgSet(a0,a1)); //es ortonormal


    // Final PDF: signal + background 
    RooRealVar Ns("Ns","Ns", 0.,12747);
    RooRealVar Nb("Nb","Nb", 0.,12747);   
    RooRealVar fraction("fraction","fraction", 0.8, 0., 1.);

    // Sum Gaussians 
    RooAddPdf sumGauss("sumGauss","sumGauss", RooArgList(signal1, signal2), RooArgList(fraction));
    
    // Model
    RooAddPdf MassModel("MassModel","MassModel", RooArgList(sumGauss, bkg), RooArgList(Ns,Nb));


    //------------ Fit procedure -------------------
    
    RooRandom::randomGenerator()->SetSeed(seed);
    
    // Verbose 0:
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    RooMsgService::instance().setSilentMode(kTRUE);
    RooMsgService::instance().setStreamStatus(1,false);
    RooMsgService::instance().getStream(1).removeTopic(Integration);  
    RooMsgService::instance().getStream(1).removeTopic(Minimization);  
    RooMsgService::instance().getStream(1).removeTopic(Fitting);  
    RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
    RooMsgService::instance().getStream(1).removeTopic(Optimization);
    RooMsgService::instance().getStream(1).removeTopic(ObjectHandling);
    RooMsgService::instance().getStream(1).removeTopic(Eval);
    RooMsgService::instance().Print();


    // -------------------------------------------
    // TTREE PARAMS
    // -------------------------------------------

    Double_t llh,edm, Nsig,NsigE, Nbkg,NbkgE, G1,G1E, Mu,MuE, G2,G2E, A0,A0E, A1,A1E, FS,FSE, 
    Nspull, Nsdif, Nbpull, Mupull;
    Int_t status,covQual,badNll;
  
  
    TFile *file = TFile::Open(Form("ToyMC_Bc_%1i.root", seed),"RECREATE");
    TTree *tree = new TTree("tree","tree");
    file->cd();

    tree->Branch("Nspull",&Nspull);
    tree->Branch("Nsdif",&Nsdif);
    tree->Branch("Nbpull",&Nbpull);    
    tree->Branch("Mupull",&Mupull);

    // ----------------------
    
    tree->Branch("Nsig",&Nsig);
    tree->Branch("NsigE",&NsigE);
    
    tree->Branch("Nbkg",&Nbkg);
    tree->Branch("NbkgE",&NbkgE);

    tree->Branch("Mu",&Mu);
    tree->Branch("MuE",&MuE);

    tree->Branch("G1",&G1);
    tree->Branch("G1E",&G1E);

    tree->Branch("G2",&G2);
    tree->Branch("G2E",&G2E);

    tree->Branch("FS",&FS);
    tree->Branch("FSE",&FSE);
    
    tree->Branch("A0",&A0);
    tree->Branch("A0E",&A0E);

    tree->Branch("A1",&A1);
    tree->Branch("A1E",&A1E);

    // -----------------------

    tree->Branch("edm",&edm);
    tree->Branch("llh",&llh);
    tree->Branch("status",&status);
    tree->Branch("covQual",&covQual);
    tree->Branch("badNll",&badNll);


    // -------------------------------------------
    // READ PARAMS FILE FROM fitBc
    // -------------------------------------------
    Double_t nSignal, nSignalE, nBkg, nBkgE, a0i, a0ie,  a1i, a1ie, fracSignal, fracSignalE;
    Double_t mui, muie, g1i, g1ie, g2i, g2ie;


    ifstream inputFile (Form("Param_Fit_Bc.txt"));
    if ( !inputFile ) 
    {
        cout << "No se pudo abrir el archivo" << endl;
        exit( 1 );
    }

    inputFile >> nSignal >> nSignalE >> nBkg >> nBkgE >> a0i >> a0ie >> a1i >> a1ie >> fracSignal >> fracSignalE >> mui >> muie >> g1i >> g1ie >> g2i >> g2ie ;
    inputFile.close();
    cout << "Ns, Nb, g2 values: " <<  nSignal << " " << nBkg << "  " << g2i << endl;



    // -------------------------------------------
    // FITTING
    // -------------------------------------------

  
    for(Int_t n=0; n<n_total; n++)
    {
        //Ns.setConstant(kTRUE);
        Ns.setVal(nSignal);
        Ns.setError(nSignalE);
        
        Nb.setVal(nBkg);
        Nb.setError(nBkgE);

        a0.setVal(a0i);
        a0.setError(a0ie);

        a1.setVal(a1i);
        a1.setError(a1ie);
            
        fraction.setVal(fracSignal);
        fraction.setError(fracSignalE);
            
        meanGauss.setVal(mui);
        meanGauss.setError(muie);
        
        sigma1.setVal(g1i);
        sigma1.setError(g1ie);
        
        sigma2.setVal(g2i);
        sigma2.setError(g2ie);


        // Generate data with the Mass model
        RooDataSet *dataToy = MassModel.generate(RooArgSet(M), Extended(kTRUE));
 
        // Fit data 
        RooFitResult* fitting = MassModel.fitTo(*dataToy,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));


        // Params and errors 
        Nspull = (Ns.getVal()-nSignal)/Ns.getError();
        Nsdif = (Ns.getVal()-nSignal);
        Nbpull = (Nb.getVal()-nBkg)/Nb.getError();   
        Mupull = (meanGauss.getVal()-mui)/meanGauss.getError();

        Nsig = Ns.getVal();
        NsigE = Ns.getError();
        
        Nbkg = Nb.getVal();
        NbkgE = Nb.getError();

        G1 = sigma1.getVal();
        G1E = sigma1.getError();

        Mu = meanGauss.getVal();
        MuE = meanGauss.getError();

        FS = fraction.getVal();
        FSE = fraction.getError();
        
        G2 = sigma2.getVal();
        G2E = sigma2.getError();

        A0 = a0.getVal();
        A0E = a0.getError();

        A1 = a1.getVal();
        A1E = a1.getError();

        // -----------------
        
        llh = fitting->minNll();
        status = fitting->status();
        covQual = fitting->covQual();
        edm = fitting->edm();
        badNll = fitting->numInvalidNLL();

        // Fill tree with results 
        tree->Fill();
        
        // Create one canvas with the toy montecarlo
        if(n==0)
        {
        TCanvas* canvas_toy_mc = CreateCanvas("canvas_toy_mc", fitting, 
        dataToy, M, Mmax, Mmin, MassModel, sumGauss, bkg, Ns, Nb, 
        sigma1, sigma2, fraction, meanGauss); 

        canvas_toy_mc->Print(Form("plots/mass_BcFit_ToyMC.png"));
        }

        delete dataToy;
        delete fitting;

        cout << "N runs: " << n << endl;

    }

    tree->Write();
    

}


