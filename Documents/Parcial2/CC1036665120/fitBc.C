#include <iostream>
using namespace RooFit;

// Class generated with treeBc->MakeClass("dataBcTree");
#include "dataBcTree.C"


// Save result of params in file 
void save_resultTfit(ofstream& salida, RooRealVar Ns, RooRealVar Nb, RooRealVar a0, RooRealVar a1 , RooRealVar fs, RooRealVar mean,  RooRealVar width, RooRealVar width2)
{
  salida <<  Ns.getVal() << " " << Ns.getError() << " " << Nb.getVal() << " " << Nb.getError() <<" " <<  a0.getVal() << " " 
    <<  a0.getError() << " " << a1.getVal() << " " << a1.getError() << " "  
	 << fs.getVal() << " " << fs.getError() << " " << mean.getVal() << " " << mean.getError() << " "
	 <<  width.getVal() << " " << width.getError() << " " <<  width2.getVal() << " " << width2.getError()
    ;
  cout << " Se ha creado el archivo ... " << endl; 
  return;
}



// -------------------------------------------
// FIT Bc
// -------------------------------------------
void fitBc()
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
    RooRealVar M("M", " M(J/#psi #pi^{+}) [GeV]", Mmin, Mmax); 
    RooDataSet data("data", "data", RooArgSet(M));

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

        data.add(RooArgSet(M));
    }

    cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << endl;
    data.Print("v");



    // -------------------------------------------
    // FITTING
    // -------------------------------------------
  
    // Create two Gaussian PDFs for the signal 
    RooRealVar meanGauss("#muB"," Mass mean Bc", 6.275, 6.250, 6.300,"GeV");
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

    // Fitting
    Ns.setVal(12747.0); 
    Nb.setVal(12747.0);
    
    RooFitResult* fitting = MassModel.fitTo(data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));

    data.Print("v"); 
    fitting->Print("v");

    // Save parameters in file
    ofstream outFileFit(Form("Param_Fit_Bc.txt"));
    outFileFit.is_open();
    save_resultTfit(outFileFit, Ns, Nb, a0, a1, fraction, meanGauss, sigma1, sigma2);



    // -------------------------------------------
    // PLOT DATA & FIT
    // -------------------------------------------

    int H = 800;
    int W = 1000;
    TCanvas *c1 = new TCanvas();
    c1->cd();

    c1->SetLeftMargin(0.005);
    c1->SetRightMargin(0.01);
    c1->SetTopMargin(0.09);
    c1->SetBottomMargin(0.1);

    // Two Pads = fit + pull 
    TPad *pad1 = new TPad("pad1", "pad1", 0.01, 0.41, 0.9903769, 0.99);
    pad1->SetLeftMargin(0.09);   
    pad1->SetRightMargin(0.019);
    pad1->SetTopMargin(0.09);
    pad1->SetBottomMargin(0.0);  

    TPad *pad2 = new TPad("pad2", "pad2", 0.01, 0.01, 0.9903769, 0.41);
    pad2->SetLeftMargin(0.09);
    pad2->SetRightMargin(0.019);  
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.25);
    pad2->SetFillColor(0);
    pad2->SetGridx(0);
    pad2->SetGridy(0);

    pad1->Draw();
    pad2->Draw();
    
    // --- PAD 1 ---
    pad1->cd();

    Double_t supM = Mmax;
    Double_t infM = Mmin;
    Double_t nbin = ((supM - infM) / 0.010);

    RooPlot* Mframe = M.frame(infM,supM,nbin);

    // Data
    data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0));

    // Model
    MassModel.plotOn(Mframe);
    MassModel.plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));
    MassModel.plotOn(Mframe, Components(sumGauss),LineColor(6),LineWidth(2),LineStyle(kDashed),Name("signal")); 
    
    // Background 
    MassModel.plotOn(Mframe,Components(bkg),LineColor(kGreen),LineWidth(2),Name("bkg"));
    
    data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("data"));
    MassModel.plotOn(Mframe);

    Mframe->SetYTitle("Events / 10 MeV");
    Mframe->SetLabelSize(0.06, "XY");
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
    TLegend *legend1 = new TLegend(0.18,0.60,0.38,0.88); 
    legend1->SetTextSize(0.07);
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
    
    TLegend *legend2 = new TLegend(0.6,0.57,0.8,0.88);
    legend2->SetTextSize(0.07);
    legend2->SetFillColor(0);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->AddEntry("",Form("M(B_{c}^{+}) = %1.2f #pm %1.2f MeV", Mpsi, MpsiE),"");
    legend2->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV", G, GE),"");
    legend2->AddEntry("",Form("N_{B_{c}^{+}} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
    legend2->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
    legend2->Draw();



    // -------------------------------------------
    // PLOT PULL
    // -------------------------------------------
    
    RooHist* pullBc = Mframe->pullHist();

    // --- PAD 2 ---
    pad2->cd();
    
    // Frame with the pull distribution 
    RooPlot* framePull = M.frame(infM,supM,nbin);
    framePull->addPlotable(pullBc,"P");

    framePull->SetYTitle("(Data-Fit)/#sigma");
    framePull->SetLabelSize(0.1,"XY");
    framePull->SetTitleSize(0.13,"X");
    framePull->SetTitleSize(0.11,"Y");
    framePull->GetYaxis()->CenterTitle();
    framePull->GetXaxis()->CenterTitle();
    framePull->GetYaxis()->SetNdivisions(505,1);
    framePull->GetXaxis()->SetNdivisions(510,1);
    framePull->GetXaxis()->SetTickLength(0.07);
    framePull->SetTitleOffset(0.9,"X");
    framePull->SetTitleOffset(0.4,"Y");
    framePull->SetMaximum(3.9);
    framePull->SetMinimum(-3.9);
    
    framePull->Draw();

    // line at 0 
    TLine *line1 = new TLine(infM,0.0,supM,0.0);
    line1->SetLineColor(1);
    line1->SetLineWidth(1);
    line1->Draw();


    c1->Modified();
    gPad->Update();
    gPad->RedrawAxis();


    c1->Print(Form("plots/mass_BcFit.png"));




}