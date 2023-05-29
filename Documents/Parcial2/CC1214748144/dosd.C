using namespace RooFit;

void dosd(){
    Float_t mass, time, error;

    TFile f("bc2.root", "READ");
    TTree *tree = (TTree*)f.Get("tree");

    tree->SetBranchAddress("Bcmass", &mass);
    tree->SetBranchAddress("tau", &time);
    tree->SetBranchAddress("tauerror", &error);

    RooRealVar M("M", "Masa [GeV]", 6.05, 6.5); 
    RooRealVar tau("tau","tiempo [ps]",0.3,2.6);
  

    RooDataSet data("data","data",RooArgSet(M,tau));;


    int nentries = tree->GetEntries();

    for (int evt=0; evt < nentries; evt++) 
    {
        tree->GetEvent (evt);
        if(mass<=6.0 || mass>=6.5) continue;
        if((time/error)<5.0) continue;
        if(time<0.3 || time>2.6) continue;
        if(error<0.0001 || error>0.5) continue;
        M=mass;
        tau=time;   
        data.add(RooArgSet(M,tau)); //
    }



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

    RooFitResult* fitres = pdf2d.fitTo(data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4), ConditionalObservables(decay,tau));

    

    TH1 *hh_decay = pdf2d.createHistogram(" ", M, Binning(50), YVar(tau, Binning(50)));

    //Gráfic
    RooPlot *frame1 = tau.frame(Title("Ajuste de la vida media"));
    data.plotOn(frame1);
    pdf2d.plotOn(frame1);

    RooPlot *frame2 = M.frame(Title("Masa invariante"));
    data.plotOn(frame2);
    pdf2d.plotOn(frame2);



    TCanvas *c = new TCanvas("c", "c", 1200, 400);


    c->Divide(3);

    c->cd(1);
    gPad->SetLeftMargin(0.20);
    hh_decay->SetTitleOffset(2.5,"Z");
    hh_decay->SetTitleOffset(2,"Y");
    hh_decay->SetTitleOffset(2,"X");
    hh_decay->SetTitle("masa invariante y la vida media");
    hh_decay->SetLabelSize(0.03,"X");
    hh_decay->GetZaxis()->SetNdivisions(6,1);
    hh_decay->Draw("surf");

    c->cd(2);
    gPad->SetLeftMargin(0.15);

    frame1->GetYaxis()->SetTitleOffset(1.6);
    frame1->Draw();

    TLegend *leg1 = new TLegend(0.48,0.75,0.8,0.88); 
    leg1->SetFillStyle(0);
    leg1->SetTextSize(0.04);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->AddEntry("",Form("#tau = %1.3f #pm %1.3f ps",-1/tmed.getVal(),1/tmed.getError()),"");
    leg1->Draw();

    c->cd(3);
    gPad->SetLeftMargin(0.15);
    frame2->GetYaxis()->SetTitleOffset(1.6);
    frame2->Draw();

    TLegend *leg2 = new TLegend(0.42,0.75,0.8,0.88); 
    leg2->SetFillStyle(0);
    leg2->SetTextSize(0.03);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry("",Form("M(B_{c})= %1.2f #pm %1.2f MeV",mean.getVal()*1000,mean.getError()*1000),"");
    leg2->Draw();

    c->Draw();
    c->SaveAs("2D.svg");
}