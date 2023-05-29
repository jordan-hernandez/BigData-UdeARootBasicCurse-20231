#include<iostream>
 
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

{
    using namespace RooFit;
    /*Lectura TTree------------------------------------------------------------------------------------------------------*/
    TFile file("Data_b_Hadron.root", "READ");
    if (!file.IsOpen()) {
        std::cout << "Error: No se puedo abrir el TFile!!" << std::endl;
        exit(1) ;
    }

    TTree *tree = (TTree*) file.Get("Data_tree");
    if (!tree) {
        std::cout << "Error: No se pudo obtener el TTree del TFile!" << std::endl;
        file.Close();
        exit(1);}
    
    // Construimos nuestro objeto DataSet para la masa y el tau
    RooRealVar Mass("Mass", "Mass", 6.05, 6.5);
    RooRealVar Tau("Tau", "Tau", 0.3, 2.6); //rangos sugeridos
    RooRealVar eTau("eTau", "eTau", 0.0001, 0.5); 
    RooDataSet Data("Data", "Data", tree, RooArgSet(Mass , Tau , eTau ));

   // Crear un conjunto de datos que contenga las variables Mass y Tau
    RooDataSet* Data_Mass_Tau = static_cast<RooDataSet*>(Data.reduce(RooArgSet(Mass, Tau)));

    // Crear un conjunto de datos que contenga las variables Tau y eTau 
    RooDataSet* Data_Tau_eTau = static_cast<RooDataSet*>(Data.reduce(RooArgSet(Tau, eTau)));
    //construccion del modelo Masa,Lifetime:----------------------------------------------------------------------------
    
    //Masa:
    RooRealVar Mean("Mean", "Mean", 6, 6.6);
    RooRealVar Sigma("Sigma", "Sigma", 0.001, 1);
    RooGaussian Signal("Signal", "Signal", Mass, Mean, Sigma);

    RooRealVar C("C", "C", -10, 10);
    RooExponential Backg("Backg", "Backg", Mass, C);

    RooRealVar Nsig("Nsig", "Nsig", 0, 13000);
    RooRealVar Nbkg("Nbkg", "Nbkg", 0, 13000);
    RooAddPdf MassModel("MassModel", "MassModel", RooArgList(Signal, Backg), RooArgList(Nsig, Nbkg));

    //Life Time model 
    RooRealVar T("T", "T", 1,-10, 10);
    RooExponential LTModel("LTModel", "LTModel", Tau, T);

    //Producto de las dos pdf
    RooProdPdf Model2D("Model2D", "Model2D", MassModel, Conditional(LTModel , Tau));
 
    /*Fiteo y Graficación-----------------------------------------------------------------------------------------------*/
    Int_t Nbins =  100; //Número de bins para la graficación
    RooFitResult* FitResult2DModel = Model2D.fitTo(*Data_Mass_Tau,Save(true),Extended(true)); 
    
    RooPlot* Massframe = Mass.frame(Title("Mass"));
    RooPlot* LTframe = Tau.frame(Title("Life Time"));

    Data_Mass_Tau->plotOn(Massframe);
    Model2D.plotOn(Massframe, ProjWData(Tau,*Data_Mass_Tau),LineColor(38));

    Data_Mass_Tau->plotOn(LTframe);
    Model2D.plotOn(LTframe, ProjWData(Mass,*Data_Mass_Tau),LineColor(38));

    /*Grafico Masa  LifeTime----------------------------------------------------------------------------------------------------------*/
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    TH1 *Hist2D = Model2D.createHistogram("Hist2D", Mass, Binning(Nbins), YVar(Tau, Binning(Nbins)));
    Hist2D->SetTitle("Mass and LifeTime");
    Hist2D->GetZaxis()->SetTitleOffset(2.5);
    Hist2D->GetXaxis()->SetTitleOffset(2.);
    Hist2D->GetYaxis()->SetTitleOffset(2.);
    Hist2D->GetXaxis()->SetTitle("Mass[GeV]");
    Hist2D->GetYaxis()->SetTitle("Tau[ps]");
    Hist2D->Draw("SURF3");
    Data_Mass_Tau->Draw("same");

    c1->Draw();
    c1->SaveAs("Mass-LifeTime.png");
    delete c1;

    /*Grafico Masa y lifetime--------------------------------------------------------------------------------------------------------------------*/
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    c2->Divide(2);

    c2->cd(1); 
    Massframe ->GetYaxis()->SetTitle((std::string(Massframe->GetYaxis()->GetTitle()) + "[GeV]").c_str());
    Massframe ->GetYaxis()->CenterTitle(); 
    //Massframe ->GetYaxis()->SetTitleSize(0.07); 
    Massframe ->GetYaxis()->SetTitleOffset(1.3);
    Massframe ->GetYaxis()->SetLabelSize(0.045);
    Massframe ->SetXTitle("M(J/#psi #pi^{+})[Gev]");
    Massframe->Draw();
    
    c2->cd(2); 
    LTframe->GetYaxis()->SetTitle((std::string(LTframe->GetYaxis()->GetTitle()) + "[ps]").c_str());
    LTframe->GetYaxis()->CenterTitle(); 
    LTframe->GetYaxis()->SetTitleOffset(1.3);
    LTframe ->GetYaxis()->SetLabelSize(0.045);
    LTframe ->SetXTitle("Tau[ps]");
    LTframe->Draw();
    c2->Update();

    c2->Draw();
    c2->SaveAs("MassProj.png");

    delete c2;
    delete FitResult2DModel;
    delete Massframe;
    delete LTframe ;

    /*Per- errors------------------------------------------------------------------------------------------------------*/

    TCanvas *c3 = new TCanvas("c3", "c3", 800, 600);
    //decay (no adimite extended ML)
    RooRealVar bias("bias", "bias", 0, -10, 10);
    RooRealVar sigma2("sigma2", "sigma2", 1, 0.1, 10);
    RooGaussModel gm("gm", "gm", Tau, bias, sigma2, eTau);
    RooDecay decay("decay", "decay", Tau, T, gm, RooDecay::DoubleSided);

    RooFitResult* FitResult2Ddecay = decay.fitTo(*Data_Tau_eTau, Save(true),ConditionalObservables(eTau));
    TH1 *hh_decay = decay.createHistogram("decay", Tau, Binning(50), YVar(eTau, Binning(50)));
    hh_decay->SetLineColor(kBlue);

    gPad->SetLeftMargin(0.20);
    hh_decay->GetZaxis()->SetTitleOffset(2.5);
    hh_decay->GetYaxis()->SetTitleOffset(2.5);
    hh_decay->GetXaxis()->SetTitleOffset(2.5);
    hh_decay->Draw("SURF3");
    c3->Draw();
    c3->SaveAs("decay(dt|dterr).png");



    delete hh_decay;
    delete FitResult2Ddecay;
    delete c3;
    

    }

