#include <iostream>
using namespace RooFit;

// Class generated with treeBc->MakeClass("dataBcTree");
#include "dataBcTree.C"

void multidimensional()
{


    // -------------------------------------------
    // READ TTREE 
    // -------------------------------------------

    // tree with the data
    TChain *chain = new TChain("treeBc", "");
    chain->Add("dataBc.root");

    TTree *treeBc = (TTree *)chain;

    // instance of the class dataBcTree with the tree
    dataBcTree t(treeBc);
    Long64_t nentries = t.fChain->GetEntries();
    std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
    std::cout << "* This tree has " << nentries << " entries.\n\n";

    //------------------------------

    std::cout << "These are the columns Det, En, and Time:" << std::endl;
    for (auto branch : *treeBc->GetListOfBranches())
    {
        std::cout << "Branch: " << branch->GetName() << std::endl;
    }

    //------------------------------

    // Mass Bc 
    Double_t Mmin = 6.05;
    Double_t Mmax = 6.5;
    RooRealVar M("M", " M [GeV]", Mmin, Mmax); 


    // Lifetime Bc
    Double_t Tmin = 0.3;
    Double_t Tmax = 2.6;
    RooRealVar T("T", " T", Tmin, Tmax); 

    // Lifetime Error
    Double_t Te_min = 0.0001;
    Double_t Te_max = 0.5;
    RooRealVar Terr("Terr", " Terr", Te_min, Te_max); 
    
    // Dataset 
    RooDataSet data("data", "data", RooArgSet(M, T));


    Int_t nTen = nentries/10;
    Int_t k=0;
    Int_t nbytes = 0, nb = 0;

    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = t.LoadTree(jentry);
        if (ientry < 0) break;
        nb = t.fChain->GetEntry(jentry);   nbytes += nb;
        if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
        if(jentry==nentries-1) cout<<endl;

        if((t.Bcmass/t.tau)<5.0)continue;

        M = t.Bcmass;
        T = t.tau;
        Terr = t.tauerror;

        data.add(RooArgSet(T, Terr));
        //data.add(RooArgSet(M, T));
    }

    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    data.Print("v");

    // -------------------------------------------
    // MULTI--DIMENSIONAL MODEL
    // -------------------------------------------

    // Build a gaussian resolution model scaled by the per-event error = gauss(t,bias,sigma*terr)
    RooRealVar bias("bias", "bias", 0, -10, 10);
    RooRealVar sigma("sigma", "per-event error scale factor", 1, 0.1, 10);
    RooGaussModel gm("gm1", "gauss model scaled bt per-event error", T, bias, sigma, Terr);
 

    // Construct gauss1(t|terr)
    RooRealVar tau("tau", "tau", 0.507);

    RooDecay decay_gm("decay_gm", "decay", T, tau, gm, RooDecay::DoubleSided);

    // Specify Terr as conditional observable
    decay_gm.fitTo(data, ConditionalObservables(Terr));

    // Make two-dimensional plot of conditional pdf in (T,Terr)
    TH1 *hh_decay = decay_gm.createHistogram("hh_decay", T, Binning(50), YVar(Terr, Binning(50)));
    hh_decay->SetLineColor(kBlue);

    // Draw all frames on canvas
    TCanvas *c = new TCanvas("canvas_multi", "canvas_multi", 1200, 400);
    c->cd();
    gPad->SetLeftMargin(0.20);
    hh_decay->GetZaxis()->SetTitleOffset(2.5);
    hh_decay->Draw("surf");
    c->Print(Form("plots/multi.png"));

}