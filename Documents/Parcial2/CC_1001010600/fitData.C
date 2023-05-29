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
void fitData(){
    Float_t mass, tau, tErr;

    TFile f("BcData.root", "READ");
    TTree *tree = (TTree*) f.Get("tree");

    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("tau", &tau);
    tree->SetBranchAddress("tErr", &tErr);

    RooRealVar MASS("MASS", "Mass (GeV)", 6.05, 6.5); //Observable
    RooDataSet massDS("dataSet", "dataSet", RooArgSet(MASS));

    //Se crea dataset para la masa
    for(int event = 0; event < tree->GetEntries(); event++){
        tree->GetEntry(event);
        if(mass<=6.0 || mass>=6.5) continue;
        if((tau/tErr)<5.0) continue;

        MASS = mass;
        massDS.add(RooArgSet(MASS));
    }

    //Señal
    RooRealVar mean1("mean1", "mean1", 6.27, 6.25, 6.3);
    RooRealVar sigma1("sigma1", "sigma1", .015, 0.005, 0.05);
    RooGaussian sig1("sig1", "sig1", MASS, mean1, sigma1);

    // RooRealVar mean2("mean2", "mean2", 6.2, 6.15, 6.25);
    // RooRealVar sigma2("sigma2", "sigma2", .05, 0, 0.1);
    // RooGaussian sig2("sig2", "sig2", MASS, mean2, sigma2);

    // RooRealVar sig1frac("sig1frac", "sig1frac", 0.8, 0., 1.);
    // RooAddPdf signal("sig", "sig", RooArgList(sig1, sig2), sig1frac);

    //Background
    RooRealVar a0("a0", "a0", -10,10);
    RooRealVar a1("a1", "a1", -10,10);

    RooChebychev bkg("bkg", "bkg", MASS, RooArgList(a0, a1));

    //Numero de eventos
    RooRealVar nsig("nsig","nsig",1000,0.,13000);
    RooRealVar nbkg("nbkg","nbkg",2000,0.,13000);
    
    //Pdf total
    RooAddPdf model("model","sig+bkg", RooArgList(sig1,bkg), RooArgList(nsig,nbkg));

    //Fit
    RooFitResult* result = model.fitTo(massDS, Extended(kTRUE), Save(kTRUE));

    //Grafico
    TCanvas* c = new TCanvas("c", "c", 1000, 800);

    TPad *pad1 = new TPad("p1", "", 0.01,0.411,0.9903769, 0.99 );
    pad1->SetLeftMargin(0.09);   
    pad1->SetRightMargin(0.019);
    pad1->SetTopMargin(0.09);
    pad1->SetBottomMargin(0.0);  
    pad1->Draw();

    TPad *pad2 = new TPad("p2", "", 0.01,0.01,0.9903769,0.41);
    pad2->SetLeftMargin(0.09);
    pad2->SetRightMargin(0.019);  
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.25);
    pad2->Draw();

    pad1->cd();
    RooPlot* mframe = MASS.frame(6.05, 6.5, 60);

    massDS.plotOn(mframe, Name("data"));

    model.plotOn(mframe, Components(sig1), LineColor(kViolet-6), LineStyle(kDashed), Name("signal"));
    model.plotOn(mframe, Components(bkg), LineColor(kAzure+7), LineStyle(kDashed), Name("bkg"));
    model.plotOn(mframe, LineColor(kPink-4), Name("model"));

    TLegend *legMass = new TLegend(0.7,0.85,0.9,0.5);
    legMass->SetTextSize(0.08);
    legMass->SetFillColor(0);
    legMass->SetBorderSize(0);
    legMass->SetFillStyle(0);
    legMass->AddEntry(mframe->findObject("data"), "Data", "pe");
    legMass->AddEntry(mframe->findObject("model"), "Fit", "l");
    legMass->AddEntry(mframe->findObject("signal"),"Signal","l");
    legMass->AddEntry(mframe->findObject("bkg"),"Bkg","l");

    mframe->SetTitle("Mass for B_{c} meson");
    mframe->Draw();
    mframe->SetYTitle("Events");
    mframe->GetYaxis()->CenterTitle();   
    mframe->SetTitleSize(0.05, "Y");
    mframe->SetLabelSize(0.04,"XY");

    legMass->Draw();

    auto mtext = new TLatex();
    mtext->SetTextSize(0.04);
    mtext->SetTextFont(42);
    mtext->DrawLatex(6.1, 80, Form("B_{c Mass} = %1.4f #pm %1.4f GeV", mean1.getVal(), mean1.getError()));

    RooHist* massPull = mframe->pullHist();
    pad2->cd();
    RooPlot* frame2 = MASS.frame(Title("Pull Distribution")) ;
    frame2->addPlotable(massPull,"P") ;
    frame2->SetYTitle("(Data-Fit)/#sigma");
    frame2->SetTitleSize(0.05, "XY");
    frame2->GetYaxis()->CenterTitle();  
    frame2->SetLabelSize(0.06,"XY"); 
    frame2->Draw();
    
    c->Draw();
    c->SaveAs("./plots/massFit.png");

    //Guardar fit
    ofstream salida("fitData.txt");
    salida << a0.getVal() << a0.getError() << " " << a1.getVal() << " " << a1.getError() << " " << mean1.getVal() << " " << mean1.getError() <<
    " " << sigma1.getVal() << " " << sigma1.getError() << " " << nbkg.getVal() << " " << nbkg.getError() << 
    " " << nsig.getVal() << " " << nsig.getError();
    salida.close();
}