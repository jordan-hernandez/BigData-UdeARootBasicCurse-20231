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

TCanvas* CreateCanvasMass(TString cname, RooFitResult* result, RooDataSet data, RooRealVar M, 
                        Double_t supM, Double_t infM, RooProdPdf MassTauModel, RooRealVar Ns_m, RooRealVar Nb_m, 
                        RooRealVar width, RooRealVar mean, Double_t ptl, Double_t pth)
{
    Double_t nbin = 40;
    int H = 600;
    int W = 800;

    TCanvas *c1 = new TCanvas(cname,cname,50,50,W,H);
    
    c1->SetLeftMargin(0.1);
    c1->SetRightMargin(0.01);
    c1->SetTopMargin(0.09);
    c1->SetBottomMargin(0.1);
    gPad->SetLogy();

    RooPlot* frameM = M.frame(6.05, 6.5, nbin);
    data.plotOn(frameM, Name("Data"));
    MassTauModel.plotOn(frameM, DrawOption("L"),LineColor(kRed),LineWidth(2), Name("Mass component"));

    frameM->SetYTitle(Form("Events / %f MeV", (6.5-6.05)/nbin)); 
    frameM->SetLabelSize(0.02,"XY");
    frameM->SetTitleSize(0.03,"XY");
    frameM->GetYaxis()->CenterTitle();   
    frameM->GetXaxis()->CenterTitle();
    frameM->GetYaxis()->SetNdivisions(505,1);
    frameM->GetXaxis()->SetNdivisions(505,1);
    frameM->GetXaxis()->SetTickLength(0.06);    
    frameM->GetXaxis()->SetDecimals(1); 
    frameM->SetTitleOffset(0.8,"X");
    frameM->SetTitleOffset(1.0,"Y");
    frameM->SetMinimum(1.0); 
    frameM->Draw();

    TLegend *leg1 = new TLegend(0.12,0.7,0.25,0.9); 
    leg1->SetTextSize(0.04);
    leg1->SetFillColor(0);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->AddEntry(frameM->findObject("Data")," Data","ep"); 
    leg1->AddEntry(frameM->findObject("Mass component"),"2D Fit result mass","l");
    leg1->Draw();

    Double_t Mpsi = mean.getVal()*1000.0;
    Double_t MpsiE = mean.getError()*1000.0;
    Double_t G = width.getVal()*1000.0;
    Double_t GE = width.getError()*1000.0;

    TLegend *legpar1 = new TLegend(0.65,0.7,0.8,0.9);
    legpar1->SetTextSize(0.03);
    legpar1->SetTextFont(42);
    legpar1->SetFillColor(0);
    legpar1->SetBorderSize(0);
    legpar1->SetFillStyle(0);
    legpar1->AddEntry("",Form("M(B_{c}^{+}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
    legpar1->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
    legpar1->AddEntry("",Form("N_{B_{c}^{+}} = %1.0f #pm %1.0f",Ns_m.getVal(),Ns_m.getError()),"");
    legpar1->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb_m.getVal(),Nb_m.getError()),"");
    legpar1->Draw();

    c1->Modified();
    return c1;
}

TCanvas* CreateCanvasTau(TString cname, RooFitResult* result, RooDataSet data, RooRealVar Tau, 
                        Double_t supT, Double_t infT, RooProdPdf MassTauModel, RooRealVar Ns_t, RooRealVar Nb_t, 
                        RooRealVar tauc, RooRealVar bkgc, Double_t ptl, Double_t pth)
{
    Double_t nbin = 40;
    int H = 600;
    int W = 800;

    TCanvas *c1 = new TCanvas(cname,cname,50,50,W,H);
    c1->cd();
    c1->SetLeftMargin(0.1);
    c1->SetRightMargin(0.01);
    c1->SetTopMargin(0.09);
    c1->SetBottomMargin(0.1);
    gPad->SetLogy();

    RooPlot* frameT = Tau.frame();
    data.plotOn(frameT, Name("Data"));
    MassTauModel.plotOn(frameT, DrawOption("L"),LineColor(kBlue),LineWidth(2), Name("Lifetime component"));

    frameT->SetYTitle(Form("Events / %f", (2.6-0.3)/nbin)); 
    frameT->SetLabelSize(0.02,"XY");
    frameT->SetTitleSize(0.03,"XY");
    frameT->GetYaxis()->CenterTitle();   
    frameT->GetXaxis()->CenterTitle();
    frameT->GetYaxis()->SetNdivisions(505,1);
    frameT->GetXaxis()->SetNdivisions(505,1);
    frameT->GetXaxis()->SetTickLength(0.06);    
    frameT->GetXaxis()->SetDecimals(1); 
    frameT->SetTitleOffset(0.8,"X");
    frameT->SetTitleOffset(1.0,"Y");
    frameT->SetMinimum(1.0); 
    frameT->Draw();

    TLegend *leg2 = new TLegend(0.65,0.5,0.8,0.69); 
    leg2->SetTextSize(0.04);
    leg2->SetFillColor(0);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(frameT->findObject("Data")," Data","ep"); 
    leg2->AddEntry(frameT->findObject("Lifetime component"),"2D Fit result lifetime","l");
    leg2->Draw();

    Double_t Tpsi = -1/tauc.getVal();
    Double_t TpsiE = -Tpsi*tauc.getError()/tauc.getVal();
    Double_t Tbkg = -1/bkgc.getVal();
    Double_t TbkgE = -Tbkg*bkgc.getError()/bkgc.getVal();

    TLegend *legpar2 = new TLegend(0.55,0.7,0.7,0.9);
    legpar2->SetTextSize(0.03);
    legpar2->SetTextFont(42);
    legpar2->SetFillColor(0);
    legpar2->SetBorderSize(0);
    legpar2->SetFillStyle(0);
    legpar2->AddEntry("",Form("#tau(B_{c}^{+}) = %1.4f #pm %1.4f ps",Tpsi,TpsiE),"");
    legpar2->AddEntry("",Form("#tau_{bkg} = %1.4f #pm %1.4f ps",Tbkg,TbkgE),"");
    legpar2->AddEntry("",Form("N_{B_{c}^{+}} = %1.0f #pm %1.0f",Ns_t.getVal(),Ns_t.getError()),"");
    legpar2->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb_t.getVal(),Nb_t.getError()),"");
    legpar2->Draw();

    c1->Modified();
    gPad->Update();
    gPad->RedrawAxis();
    return c1;
}

// TCanvas* CreateCanvasPlot(TString cname, RooFitResult* result, RooDataSet data, RooRealVar M, 
//                         Double_t supM, Double_t infM, RooRealVar Tau, 
//                         Double_t supT, Double_t infT,  RooProdPdf MassTauModel, RooAddPdf Masssignal, 
//                         RooRealVar Ns_m, RooRealVar Nb_m, RooRealVar width, 
//                         RooRealVar mean, RooAddPdf Tausignal, RooRealVar Ns_t, RooRealVar Nb_t, 
//                         RooRealVar tauc, RooRealVar bkgc, Double_t ptl, Double_t pth)
// {
//     Double_t nbin = 40;
//     int H = 600;
//     int W = 800;

//     TCanvas *c1 = new TCanvas(cname,cname,50,50,W,H);
//     c1->cd();
//     c1->SetLeftMargin(0.005);
//     c1->SetRightMargin(0.01);
//     c1->SetTopMargin(0.09);
//     c1->SetBottomMargin(0.1);

//     TPad *pad1 = new TPad("pad1", "padi", 0.01, 0.1, 0.49, 0.99);
//     //TPad *pad1 = new TPad("pad1", "padi",0.01,0.411,0.9903769, 0.99 );
//     pad1->SetLeftMargin(0.09);   
//     pad1->SetRightMargin(0.019);
//     pad1->SetTopMargin(0.09);
//     pad1->SetBottomMargin(0.0);  

//     TPad *pad2 = new TPad("pad2", "pad2", 0.51, 0.1, 0.99, 0.99);
//     //TPad *pad2 = new TPad("pad2", "pad2",0.01,0.01,0.9903769,0.41);
//     pad2->SetLeftMargin(0.09);
//     pad2->SetRightMargin(0.019);  
//     pad2->SetTopMargin(0.0);
//     pad2->SetBottomMargin(0.25);
//     //pad2->SetTickx(0);
//     pad2->SetFillColor(0);
//     pad2->SetGridx(0);
//     pad2->SetGridy(0);

//     pad1->Draw();
//     pad2->Draw();

//     //El el pad izquierdo se plotea la masa.
//     pad1->cd();

//     RooPlot* frameM = M.frame();
//     data.plotOn(frameM, Name("Data"));
//     MassTauModel.plotOn(frameM, DrawOption("L"),LineColor(kRed),LineWidth(2), Name("Mass component"));

//     frameM->SetYTitle("Events / 2.5 MeV"); 
//     frameM->SetLabelSize(0.02,"XY");
//     frameM->SetTitleSize(0.03,"XY");
//     frameM->GetYaxis()->CenterTitle();   
//     frameM->GetXaxis()->CenterTitle();
//     frameM->GetYaxis()->SetNdivisions(505,1);
//     frameM->GetXaxis()->SetNdivisions(505,1);
//     frameM->GetXaxis()->SetTickLength(0.06);    
//     frameM->GetXaxis()->SetDecimals(1); 
//     frameM->SetTitleOffset(0.8,"X");
//     frameM->SetTitleOffset(1.0,"Y");
//     frameM->SetMinimum(1.0); 
//     frameM->Draw();

//     TLegend *leg1 = new TLegend(0.12,0.7,0.25,0.9); 
//     leg1->SetTextSize(0.04);
//     leg1->SetFillColor(0);
//     leg1->SetBorderSize(0);
//     leg1->SetFillStyle(0);
//     leg1->AddEntry(frameM->findObject("Data")," Data","ep"); 
//     leg1->AddEntry(frameM->findObject("Mass component"),"2D Fit result mass","l");
//     leg1->Draw();

//     Double_t Mpsi = mean.getVal()*1000.0;
//     Double_t MpsiE = mean.getError()*1000.0;
//     Double_t G = width.getVal()*1000.0;
//     Double_t GE = width.getError()*1000.0;

//     TLegend *legpar1 = new TLegend(0.55,0.7,0.7,0.9);
//     legpar1->SetTextSize(0.03);
//     legpar1->SetTextFont(42);
//     legpar1->SetFillColor(0);
//     legpar1->SetBorderSize(0);
//     legpar1->SetFillStyle(0);
//     legpar1->AddEntry("",Form("M(B_{c}^{+}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
//     legpar1->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
//     legpar1->AddEntry("",Form("N_{B_{c}^{+}} = %1.0f #pm %1.0f",Ns_m.getVal(),Ns_m.getError()),"");
//     legpar1->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb_m.getVal(),Nb_m.getError()),"");
//     legpar1->Draw();

//     //El el pad derecho se plotea la masa.
//     pad2->cd();

//     RooPlot* frameT = Tau.frame();
//     data.plotOn(frameT, Name("Data"));
//     MassTauModel.plotOn(frameT, DrawOption("L"),LineColor(kBlue),LineWidth(2), Name("Lifetime component"));

//     frameT->SetYTitle("Events / 2.5 MeV"); 
//     frameT->SetLabelSize(0.02,"XY");
//     frameT->SetTitleSize(0.03,"XY");
//     frameT->GetYaxis()->CenterTitle();   
//     frameT->GetXaxis()->CenterTitle();
//     frameT->GetYaxis()->SetNdivisions(505,1);
//     frameT->GetXaxis()->SetNdivisions(505,1);
//     frameT->GetXaxis()->SetTickLength(0.06);    
//     frameT->GetXaxis()->SetDecimals(1); 
//     frameT->SetTitleOffset(0.8,"X");
//     frameT->SetTitleOffset(1.0,"Y");
//     frameT->SetMinimum(1.0); 
//     frameT->Draw();

//     TLegend *leg2 = new TLegend(0.55,0.5,0.7,0.69); 
//     leg2->SetTextSize(0.04);
//     leg2->SetFillColor(0);
//     leg2->SetBorderSize(0);
//     leg2->SetFillStyle(0);
//     leg2->AddEntry(frameT->findObject("Data")," Data","ep"); 
//     leg2->AddEntry(frameT->findObject("Lifetime component"),"2D Fit result lifetime","l");
//     leg2->Draw();

//     Double_t Tpsi = 1/tauc.getVal();
//     Double_t TpsiE = Tpsi*tauc.getError()/tauc.getVal();
//     Double_t Tbkg = 1/bkgc.getVal();
//     Double_t TbkgE = Tbkg*bkgc.getError()/bkgc.getVal();

//     TLegend *legpar2 = new TLegend(0.55,0.7,0.7,0.9);
//     legpar2->SetTextSize(0.03);
//     legpar2->SetTextFont(42);
//     legpar2->SetFillColor(0);
//     legpar2->SetBorderSize(0);
//     legpar2->SetFillStyle(0);
//     legpar2->AddEntry("",Form("#tau(B_{c}^{+}) = %1.2f #pm %1.2f MeV",Tpsi,TpsiE),"");
//     legpar2->AddEntry("",Form("#tau_{bkg} = %1.2f #pm %1.2f MeV",Tbkg,TbkgE),"");
//     legpar2->AddEntry("",Form("N_{B_{c}^{+}} = %1.0f #pm %1.0f",Ns_t.getVal(),Ns_t.getError()),"");
//     legpar2->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb_t.getVal(),Nb_t.getError()),"");
//     legpar2->Draw();

//     c1->cd();

//     TLatex *tex3 = new TLatex(0.22,0.95,"Preliminary");
//     tex3->SetNDC();
//     tex3->SetTextFont(52);
//     tex3->SetTextSize(0.04); 
//     tex3->SetLineWidth(2);
//     tex3->Draw();

//     c1->Modified();
//     gPad->Update();
//     gPad->RedrawAxis();
//     gPad->SetLogy();

//     return c1;
// }

void Bcmass_Tau(Double_t ptl=12.0, Double_t pth=70.0){

    Float_t Mmax = 6.50;
    Float_t Mmin = 6.04;
    Float_t tmin = 0.3;
    Float_t tmax = 2.6;
    Float_t terrmin = 0.0001;
    Float_t terrmax = 0.5;

    TChain *ch = new TChain("bctree","");
    ch->Add("../Bc_data.root/tree");

    TTree *tree = (TTree*) ch;
    myTreeClass t(tree);
    Long64_t nentries = t.fChain->GetEntries();
    cout<<" Entries : "<<nentries<<endl;

    //------------------------------
    RooRealVar M("M","M(J/#psi #pi^{+}) (MeV)",Mmin,Mmax);
    RooRealVar Tau("Tau","#tau(J/#psi #pi^{+}) (ps)",tmin,tmax);
    RooRealVar TauErr("TauError","#Delta#tau(J/#psi #pi^{+}) (ps)",terrmin,terrmax);
    RooDataSet data("data","Filtered data Bc mass/lifetime",RooArgSet(M, Tau, TauErr));

    Int_t nTen = nentries/10;
    Int_t k=0;
    Int_t nbytes = 0, nb = 0;
    for(Long64_t jentry=0; jentry<nentries;jentry++){
        Long64_t ientry = t.LoadTree(jentry);
        if (ientry < 0) break;
        nb = t.fChain->GetEntry(jentry);   nbytes += nb;
        if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
        if(jentry==nentries-1) cout<<endl;

        //Mass windows cuts
        if(t.Bcmass<=Mmin || t.Bcmass>=Mmax) continue;

        //Lifetime windows cuts
        if(t.Tau<=tmin || t.Tau>=tmax) continue;

        //Lifetime error window cuts	
        if(t.TauError<=terrmin || t.TauError>=terrmax) continue;

        //Filtro para lifetime.
        if((t.Tau/t.TauError)<5.0)continue;

        M = t.Bcmass;
        Tau = t.Tau;
        TauErr = t.TauError;    
        data.add(RooArgSet(M, Tau, TauErr));

    }

    cout<<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<endl;
    data.Print("v");

    ////////////////////////////////////////
    //---- C O M B I N E D  M O D E L ----//
    ////////////////////////////////////////
    Double_t supM = Mmax;
    Double_t infM = Mmin;
    Double_t supT = tmax;
    Double_t infT = tmin;
    Double_t supTerr = terrmax;
    Double_t infTerr = terrmin;

    //-------------- Mass Model -------------------
    //Modelo para el background: exponencial.
    RooRealVar m_c_bkg("m_c_bkg","Mass bkg c",-20.0,0.0);
    RooExponential Mbkg("Mbkg","Exp. Background",M,m_c_bkg);

    //Modelo para la señal: gaussiana.
    RooRealVar mean("mean"," Mass mean",6.275,6.25,6.3,"MeV");
    RooRealVar width("width"," Mass width",0.10,0.01,0.15,"MeV");
    RooGaussian Sig("Sig"," Signal PDF",M,mean,width);

    //Modelo combinado: exponencial + gaussiana.
    RooRealVar Ns_m("Ns_mass","Ns mass",0.,nentries);
    RooRealVar Nb_m("Nb_mass","Nb mass",0.,nentries);   

    RooAddPdf MassModel("MassModel","Mass Model",RooArgList(Sig,Mbkg),RooArgList(Ns_m,Nb_m));

    //-------------- Lifetime Model -------------------
    //Modelo para el background: exponencial.
    RooRealVar lt_c_bkg("lt_c_bkg","Lifetime bkg c",-5.0,5.0, "ps^-1");
    RooExponential Tbkg("Tbkg","Exp. Lifetime bkg",Tau,lt_c_bkg);

    //Modelo para la señal: exponencial.
    RooRealVar lt_c("lt_c","Lifetime signal c",-5.0,5.0,"ps^-1");
    RooExponential TSig("TSig","Lifetime signal",Tau,lt_c);

    //Modelo combinado: exponencial + exponencial.
    RooRealVar Ns_t("Ns_lt","Ns lifetime",0.,nentries);
    RooRealVar Nb_t("Nb_lt","Nb lifetime",0.,nentries);

    RooAddPdf TauModel("TauModel","Tau Model",RooArgList(TSig,Tbkg),RooArgList(Ns_t,Nb_t));

    //-------------- Total Mass/Tau Model -------------------
    RooProdPdf MassTauModel("MassTauModel", "Mass/Tau 2D Model", RooArgSet(MassModel, TauModel));

    ////////////////////////////////////
    //----F I T  P R O C E D U R E----//
    ////////////////////////////////////

    Ns_m.setVal(nentries/2);
    Nb_m.setVal(nentries/2);
    Ns_t.setVal(nentries/2);
    Nb_t.setVal(nentries/2);

    RooFitResult* fitres = MassTauModel.fitTo(data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));
    fitres->Print("v");

    //TCanvas *c1 = new TCanvas("c1","",50,50,W,H);
    //gPad->SetLogy();

    TCanvas* canv_nominalpull = CreateCanvasMass("canv_nominalpull", fitres, data, M, supM, infM, 
                                                MassTauModel, Ns_m, Nb_m, width, 
                                                mean, ptl, pth);
    canv_nominalpull->Print(Form("Bc2DFit_masscomponent_%1.0f_%1.0f.png",ptl,pth));
    TCanvas* canv_nominalpull2 = CreateCanvasTau("canv_nominalpull", fitres, data, Tau, supT, infT, 
                                                MassTauModel, Ns_t, Nb_t, lt_c, 
                                                lt_c_bkg, ptl, pth);
    canv_nominalpull2->Print(Form("Bc2DFit_lifetimecomponent_%1.0f_%1.0f.png",ptl,pth));

    // RooPlot* frameM = M.frame();
    // data.plotOn(frameM);
    // MassTauModel.plotOn(frameM);
    // frameM->Draw();

    // RooPlot* frameT = Tau.frame();
    // data.plotOn(frameT);
    // MassTauModel.plotOn(frameT);
    // frameT->Draw();
}