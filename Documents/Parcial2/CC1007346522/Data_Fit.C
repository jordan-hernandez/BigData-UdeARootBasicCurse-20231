#include "ClassTree1.C"

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

using namespace std;
using namespace RooFit;

// Funcion para guardar los parametros del fit
void save_resultTfit(ofstream& salida, RooRealVar Ns, RooRealVar Nb, RooRealVar c, RooRealVar fs, RooRealVar mean,  RooRealVar width, RooRealVar width2)
{ 
    salida <<  Ns.getVal() << " " << Ns.getError() << " " << Nb.getVal() << " " << Nb.getError() <<" "
        <<  c.getVal() << " " << c.getError() << " "
        << fs.getVal() << " " << fs.getError() << " " << mean.getVal() << " " << mean.getError() << " "
        <<  width.getVal() << " " << width.getError() << " " <<  width2.getVal() << " " << width2.getError();
    cout << " el archivo se escribio bien" << endl; 
    return;
}

// Funcion para crear figura
TCanvas* CreateCanvasNomPull(TString cname, RooFitResult* result, RooDataSet data, RooRealVar M, Double_t supM, Double_t infM,  RooAddPdf MassModel, RooAddPdf sumgau, RooExponential bkg1, RooRealVar Ns, RooRealVar Nb, RooRealVar width, RooRealVar width2, RooRealVar fs, RooRealVar mean) 
{   
    // Número de bines
    Double_t nbin = ((supM-infM)/0.005)+1;

    int H = 600;
    int W = 800;

    // Creando Canvas
    TCanvas *c1 = new TCanvas(cname,cname,50,50,W,H);
    c1->cd() ;  
    c1->SetLeftMargin(0.005);
    c1->SetRightMargin(0.01);
    c1->SetTopMargin(0.09);
    c1->SetBottomMargin(0.1);

    // Creando pads
    TPad *pad1 = new TPad("pad1", "pad1",0.01,0.411,0.9903769,0.99);
    pad1->SetLeftMargin(0.09);   
    pad1->SetRightMargin(0.019);
    pad1->SetTopMargin(0.09);
    pad1->SetBottomMargin(0.0);  

    TPad *pad2 = new TPad("pad2", "pad2",0.01,0.01,0.9903769,0.41);
    pad2->SetLeftMargin(0.09);
    pad2->SetRightMargin(0.019);  
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.25);
    pad2->SetTickx(0);
    pad2->SetFillColor(0);
    pad2->SetGridx(0);
    pad2->SetGridy(0);

    pad1->Draw();
    pad2->Draw();
    pad1->cd(); 

    // Creando frame
    RooPlot* Mframe = M.frame(infM,supM,nbin);
    // Se dibujan los datos y el ajuste
    data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0));
    MassModel.plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));
    
    RooHist* hpullm2 = Mframe->pullHist() ; // Crear Histograma Pull
    
    // Dibujar el frame
    MassModel.plotOn(Mframe,Components(sumgau),LineColor(kRed),LineWidth(2),Name("Signal")); 
    MassModel.plotOn(Mframe,Components(bkg1),LineColor(kGreen),LineWidth(2),Name("bkg")); 
    data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("Data"));
    MassModel.plotOn(Mframe);
    Mframe->SetTitle(""); 
    
    // Ajustes del frame
    Mframe->SetYTitle("Events / 2.5 MeV"); 
    Mframe->SetLabelSize(0.07,"XY");
    Mframe->SetTitleSize(0.08,"XY");
    Mframe->GetYaxis()->CenterTitle();   
    Mframe->GetXaxis()->CenterTitle();
    Mframe->GetYaxis()->SetNdivisions(505,1);
    Mframe->GetXaxis()->SetNdivisions(505,1);
    Mframe->GetXaxis()->SetTickLength(0.0);    
    Mframe->GetXaxis()->SetDecimals(1); 
    Mframe->SetTitleOffset(0.8,"X");
    Mframe->SetTitleOffset(0.6,"Y");
    Mframe->SetMinimum(0.5); 
    Mframe->Draw();
    
    // Legend 1
    TLegend *leg = new TLegend(0.18,0.58,0.38,0.88); 
    leg->SetTextSize(0.06);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(Mframe->findObject("Data")," Data","ep"); 
    leg->AddEntry(Mframe->findObject("fittotal")," Fit result","l");
    leg->AddEntry(Mframe->findObject("Signal")," Signal","l");
    leg->AddEntry(Mframe->findObject("bkg"),"Comb. backg.","l");
    leg->Draw();
    
    // Valores y errores en valor medio y ancho de la señal
    Double_t Mpsi = mean.getVal()*1000.0;
    Double_t MpsiE = mean.getError()*1000.0;

    Double_t G = sqrt( fs.getVal()*width.getVal()*width.getVal() + (1-fs.getVal())*width2.getVal()*width2.getVal() )*1000.0;
    Double_t GE = (1/G)*sqrt( (fs.getVal()*fs.getVal())*(width.getVal()*width.getVal())*(width.getError()*width.getError()) + ((1-fs.getVal())*(1-fs.getVal()))*(width2.getVal()*width2.getVal())*(width2.getError()*width2.getError()) )*1000.0*1000.0;
    
    // Legend Parámetros
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
    
    // Pad para el Histograma Pull
    pad2->cd();
    
    RooPlot* framem2 = M.frame(infM,supM,nbin) ;
    framem2->addPlotable(hpullm2,"P") ;
    
    // Ajustes del frame
    framem2->SetTitle("");
    framem2->SetYTitle(" (Data-Fit)/#sigma");
    framem2->SetLabelSize(0.1,"XY");
    framem2->SetTitleSize(0.13,"X");
    framem2->SetTitleSize(0.11,"Y");  
    framem2->GetYaxis()->CenterTitle();   
    framem2->GetXaxis()->CenterTitle();
    framem2->GetYaxis()->SetNdivisions(505,1);
    framem2->GetXaxis()->SetNdivisions(505,1);
    framem2->GetXaxis()->SetTickLength(0.07);   
    framem2->SetTitleOffset(0.9,"X");
    framem2->SetTitleOffset(0.4,"Y");
    framem2->SetMaximum(4.9);
    framem2->SetMinimum(-4.9);
    
    framem2->Draw();

    // Linea en cero
    TLine *line1 = new TLine(infM,0.0,supM,0.0);
    line1->SetLineColor(4);
    line1->SetLineWidth(1);
    line1->Draw();
    
    c1->cd();

    TLatex *   tex1 = new TLatex(0.88,0.95,"L = 19.7 fb^{-1} (#sqrt{s} = 8 TeV)");
    tex1->SetNDC();
    tex1->SetTextAlign(31);
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.04); 
    tex1->SetLineWidth(2);
    
    TLatex *tex2 = new TLatex(0.15,0.95,"UdeA");
    tex2->SetNDC();
    tex2->SetTextFont(61);
    tex2->SetTextSize(0.04); 
    tex2->SetLineWidth(2);    

    tex1->Draw();  
    tex2->Draw();

    c1->Modified();
    gPad->Update();
    gPad->RedrawAxis();
    TLine l;
    l.DrawLine(gPad->GetUxmax(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymin());
    c1->Update();
    return c1; 
  
}

int Data_Fit(){

    // Creamos el TTree
    TChain * Chain1 = new TChain("Tree1","");
    Chain1->Add("RootFile1.root");
    TTree *mytree = (TTree*) Chain1;

    // Usando ClassTree1 desempaquetamos el Tree
    ClassTree1 Tree1(mytree);

    cout << "* This tree has " << Tree1.fChain->GetEntries() << " entries.\n\n";

    // Valores Límete para las variables
    Double_t Mmin = 6.05; 
    Double_t Mmax = 6.5;
    Double_t Tmin = 0.3; 
    Double_t Tmax = 2.6;
    Double_t T_errmin = 0.0001; 
    Double_t T_errmax = 0.5;

    // Variables a usar
    RooRealVar M("M", "Mass B_{c} (GeV)", Mmin, Mmax);
    RooRealVar T("T", "#tau", Tmin, Tmax);
    RooRealVar T_err("T_err", "Tau error", T_errmin,T_errmax);

    // Dataset con valores de M
    RooDataSet Data_M("Data", "Data", RooArgSet(M));

    // Llenando el RooDataSet
    Long64_t nentries = Tree1.fChain->GetEntries();
    Int_t nTen = nentries/10;
    Long64_t nbytes = 0, nb = 0;

    for (int jentry=0; jentry<nentries; jentry++) 
    {
        Long64_t ientry = Tree1.LoadTree(jentry);
        if (ientry < 0) break;
        if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
        if(jentry==nentries-1) cout<<endl;

        nb = Tree1.fChain->GetEntry(ientry);   nbytes += nb;
        
        // No permitir que la razon entre Tau y su error sea menor a 5
        if ((Tree1.Tau/Tree1.Tau_err)<5) continue;

        // Verificar que la masa esté en el rango establecido
        if ((Tree1.mass<Mmin)||(Tree1.mass>Mmax)) continue;
        
        M = Tree1.mass;
        Data_M.add(RooArgSet(M)); // Añadir al Dataset
 
    }

    // ---- MassModel ----

    // -Parámetros Señal-
    RooRealVar mean("mean"," Mass mean",6.27,6.25,6.3,"GeV");

    // Gausiana 1
    RooRealVar width("width"," Mass width",0.10,0.001,0.3,"GeV"); // Sigma1
    RooGaussian Sig("Sig"," Signal PDF",M,mean,width);

    // Gausiana 2
    RooRealVar width2("width2"," Mass width2 ",0.20,0.001,10.0,"GeV"); // Sigma2
    RooGaussian Sig2("Sig2"," Signal PDF B",M,mean,width2);

    // -Parámetros Background-
    RooRealVar c("c","c",-10.0,10.0);
    RooExponential Bkg("Bkg","Exp. Background",M,c);

    // Cantidad de datos por cada componente
    RooRealVar Ns("Ns","Ns",0.,2000);
    RooRealVar Nb("Nb","Nb",0.,2000);   
    RooRealVar fs("fs","fs",0.8,0.,1.);

    // Suma de las dos gausianas
    RooAddPdf Sumgaus("sumgau","sumgau",RooArgList(Sig,Sig2),RooArgList(fs));

    // Modelo de masa (2 Gausianas + Exponencial)
    RooAddPdf MassModel("MassModel","MassModel",RooArgList(Sumgaus,Bkg),RooArgList(Ns,Nb));

    // ---- Fitting ----
    RooFitResult* ResultFit = MassModel.fitTo(Data_M,Extended(),Minos(kFALSE),Save(kTRUE));

    // Guardando los parámetros encontrados
    ofstream salida_TotalFit("TotalFit.txt");
    salida_TotalFit.is_open();
    save_resultTfit(salida_TotalFit, Ns, Nb, c, fs, mean, width, width2);
 
    // Hacer Gráfica
    Double_t supM = Mmax;
    Double_t infM = Mmin;
    TCanvas* Canvas1 = CreateCanvasNomPull("Canvas_MasaBc", ResultFit, Data_M, M, supM, infM, MassModel, Sumgaus, Bkg, Ns, Nb, width, width2, fs, mean);  
    Canvas1->Print("plots/MasaBc_FiteoyPull.png");

    return 0;
}