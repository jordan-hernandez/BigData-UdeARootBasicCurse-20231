#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
using namespace RooFit;
void toyMonteCarlo(){
    Float_t a0f, a0fErr, a1f, a1fErr, meanf, meanfErr, sigmaf, sigmafErr, nbkgf, nbkgfErr, nsigf, nsigfErr;

    ifstream entrada("fitData.txt");
    if(!entrada) 
    {
      cout << "No se pudo abrir el archivo entrada fit pr bins" << endl;
      exit( 1 );
    }
    entrada >> a0f >> a0fErr >> a1f >> a1fErr >> meanf >> meanfErr >> sigmaf >> sigmafErr >> nbkgf >> nbkgfErr >> nsigf >> nsigfErr;
    entrada.close();

    RooRealVar m("m", "m", 6.05, 6.5);

    //Señal
    RooRealVar mean("mean", "mean", 6.26, 6.28); 
    RooRealVar sigma("sigma", "sigma", 0.005, 0.05); 
    RooGaussian sig("sig", "sig", m, mean, sigma);

    //Background
    RooRealVar a0("a0", "a0", -1,1); 
    RooRealVar a1("a1", "a1", -1,1); a1.setVal(a1f); 
    RooChebychev bkg("bkg", "bkg", m, RooArgList(a0, a1));

    //Numero de eventos
    RooRealVar nsig("nsig","nsig", 400,800); 
    RooRealVar nbkg("nbkg","nbkg", 400,800); 
    
    //Pdf total
    RooAddPdf model("model","sig+bkg", RooArgList(sig,bkg), RooArgList(nsig,nbkg));

    //DataSets
    RooRealVar massPull("massPUll", "massPull", -6, 6);
    RooRealVar nsigPull("nsigPull", "nsigPull", -6, 6);
    RooRealVar nbkgPull("nbkgPull", "nbkgPull", -6, 6);

    RooDataSet mpDset("massPUll", "massPull", RooArgSet(massPull));
    RooDataSet nsigDset("nsigPULL", "nsigPULL", RooArgSet(nsigPull));
    RooDataSet nbkgDset("nbkgPULL", "nbkgPULL", RooArgSet(nbkgPull));

    RooRandom::randomGenerator()->SetSeed(84329746);
    for(int i=0; i<200; i++){

      mean.setVal(meanf); mean.setError(meanfErr);
      sigma.setVal(sigmaf); sigma.setError(sigmafErr);
      a0.setVal(a0f); a0.setError(a0fErr);
      nsig.setVal(nsigf); nsig.setError(nsigfErr);
      a1.setVal(a1f); a1.setError(a1fErr);
      nbkg.setVal(nbkgf); nbkg.setError(nbkgfErr);
      
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

      RooDataSet* toyData = model.generate(RooArgSet(m), Extended(kTRUE));
      RooFitResult* result = model.fitTo(*toyData, Extended(kTRUE), Minos(kFALSE), Save(kTRUE));

      massPull.setVal((mean.getVal() - meanf) / mean.getError());
      mpDset.add(RooArgSet(massPull));

      nsigPull.setVal((nsig.getVal() - nsigf) / nsig.getError());
      nsigDset.add(RooArgSet(nsigPull));

      nbkgPull.setVal((nbkg.getVal() - nbkgf) / nbkg.getError());
      nbkgDset.add(RooArgSet(nbkgPull));
    }

    //Fit Pull masa
    RooRealVar masPullMean("masPullMean", "masPullMean", 0, -.5, .5);
    RooRealVar massPullSigma("massPullSigma", "massPullSigma", 1, 0.5, 1.5);
    RooGaussian massPullGauss("massPullGauss", "massPullGauss", massPull, masPullMean, massPullSigma);
    massPullGauss.fitTo(mpDset, Extended(kTRUE), Save(kTRUE));

    //Fit Pull nsignal
    RooRealVar sigPullMean("sigPullMean", "sigPullMean", 0, -.5, .5);
    RooRealVar sigPullSigma("sigPullSigma", "sigPullSigma", 1, 0.5, 1.5);
    RooGaussian sigPullGauss("sigPullGauss", "sigPullGauss", nsigPull, sigPullMean, sigPullSigma);
    sigPullGauss.fitTo(nsigDset, Extended(kTRUE), Save(kTRUE));

    //Fit Pull nsignal
    RooRealVar bkgPullMean("bkgPullMean", "bkgPullMean", 0, -.5, .5);
    RooRealVar bkgPullSigma("bkgPullSigma", "bkgPullSigma", 1, 0.5, 1.5);
    RooGaussian bkgPullGauss("bkgPullGauss", "bkgPullGauss", nbkgPull, bkgPullMean, bkgPullSigma);
    bkgPullGauss.fitTo(nbkgDset, Extended(kTRUE), Save(kTRUE));

    TCanvas* c2 = new TCanvas("c2", "c2", 1400, 600);
    TPad *pad1 = new TPad("p1", "", 0.01,0.05,0.32, 0.99); pad1->Draw();
    TPad *pad2 = new TPad("p2", "", 0.33,0.05,0.65, 0.99); pad2->Draw();
    TPad *pad3 = new TPad("p3", "", 0.67,0.05,0.99, 0.99); pad3->Draw();

    pad1->cd();
    
    RooPlot* mframe2 = massPull.frame(-6, 6, 30);
    mpDset.plotOn(mframe2);
    massPullGauss.plotOn(mframe2, LineColor(kBlue-9));
    
    mframe2->SetTitle("Mass pull distribution");
    mframe2->Draw();

    auto mtext = new TLatex();
    mtext->SetTextSize(0.04);
    mtext->SetTextFont(42);
    mtext->DrawLatex(0.6, 35,Form("#mu = %1.4f #pm %1.4f", masPullMean.getVal(), masPullMean.getError()));
    mtext->DrawLatex(0.6, 32,Form("#sigma = %1.4f #pm %1.4f", massPullSigma.getVal(), massPullSigma.getError()));

    pad2->cd();
    RooPlot* sframe = nsigPull.frame(-6, 6, 30);
    nsigDset.plotOn(sframe);
    sigPullGauss.plotOn(sframe, LineColor(kBlue-9));
    sframe->SetTitle("Signal event pull distribution");
    sframe->Draw();

    auto stext = new TLatex();
    stext->SetTextSize(0.04);
    stext->SetTextFont(42);
    stext->DrawLatex(0.6, 45,Form("#mu = %1.4f #pm %1.4f", sigPullMean.getVal(), sigPullMean.getError()));
    stext->DrawLatex(0.6, 43,Form("#sigma = %1.4f #pm %1.4f", sigPullSigma.getVal(), sigPullSigma.getError()));

    pad3->cd();
    RooPlot* bframe = nbkgPull.frame(-6, 6, 30);
    nbkgDset.plotOn(bframe);
    bkgPullGauss.plotOn(bframe, LineColor(kBlue-9));
    bframe->SetTitle("Bakcground event pull distribution");
    bframe->Draw();

    auto btext = new TLatex();
    btext->SetTextSize(0.04);
    btext->SetTextFont(42);
    btext->DrawLatex(0.6, 40,Form("#mu = %1.4f #pm %1.4f", bkgPullMean.getVal(), bkgPullMean.getError()));
    btext->DrawLatex(0.6, 38,Form("#sigma = %1.4f #pm %1.4f", bkgPullSigma.getVal(), bkgPullSigma.getError()));


    c2->Draw();
    c2->SaveAs("./plots/monteCarlo.png");
}