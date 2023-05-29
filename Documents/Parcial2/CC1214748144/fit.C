using namespace RooFit;

void fit(){
    float Bcmass, tau, tauerror;
    TFile f("bc2.root", "READ");
    TTree *tree = (TTree*) f.Get("tree");
    tree->SetBranchAddress("Bcmass", &Bcmass);
    tree->SetBranchAddress("tau", &tau);
    tree->SetBranchAddress("tauerror", &tauerror);

    RooRealVar M("M", "M", 6.05, 6.5); 
    RooDataSet data("dataSet", "dataSet", RooArgSet(M));

    //Se crea dataset para la masa
    for(int event = 0; event < tree->GetEntries(); event++){
        tree->GetEntry(event);
        if(Bcmass<=6.0 || Bcmass>=6.5) continue;
        if((tau/tauerror)<5) continue;
        if(tau<0.3 || tau>2.6) continue;
        if(tauerror<0.0001 || tauerror>0.5) continue;
        M.setVal(Bcmass);
        data.add(RooArgSet(M));
    }   
    

    //gaussiana de la distribución de la masa
    RooRealVar mean("mean", "mean", 6.27, 6.25, 6.3);
    RooRealVar sigma("sigma", "sigma", 0.015, 0.01, 0.05);
    RooGaussian gaus("gaus", "gaus", M, mean, sigma);

    
    //línea del background
    RooRealVar a0("a0", "a0", -10,10);
    RooRealVar a1("a1", "a1", -10,10);
    RooChebychev bkg("bkg", "bkg", M, RooArgList(a0, a1));

    //número de eventos
    RooRealVar nsig("nsig","nsig",0.,13000);
    RooRealVar nbkg("nbkg","nbkg",0.,13000);


    RooAddPdf massmodel("massmodel","massmodel",RooArgList(gaus,bkg),RooArgList(nsig,nbkg));
    RooFitResult* fitres = massmodel.fitTo(data,Extended(kTRUE),Save(kTRUE), NumCPU(4));

    //guardamos los datos resultado del fit en un archivo para luego usarlos en el montecarlo. 
    ofstream myfile;
    myfile.open ("fit.txt");
    myfile << nsig.getVal() << " " << nsig.getError()<<" " << nbkg.getVal() << " " << nbkg.getError()<<" "<<mean.getVal() << " " << mean.getError()<< " "<<sigma.getVal() << " " << sigma.getError();
    myfile <<" " <<a0.getVal() << " " << a0.getError()<< " " <<a1.getVal() << " " << a1.getError();
    myfile.close();


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
    RooPlot* mframe = M.frame(6.05, 6.5, 60);
    data.plotOn(mframe, Name("data"));

    mframe->SetYTitle(" Events");
    massmodel.plotOn(mframe, Components(gaus), LineColor(kViolet-6), LineStyle(kDashed), Name("signal"));
    massmodel.plotOn(mframe, Components(bkg), LineColor(kAzure+7), LineStyle(kDashed), Name("bkg"));
    massmodel.plotOn(mframe, LineColor(kPink-4), Name("model"));

    TLegend *legMass = new TLegend(0.18,0.58,0.38,0.88);
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
    legMass->Draw();
    TLegend *legfit = new TLegend(0.6,0.68,0.8,0.88); 
    legfit->SetTextSize(0.06);
    legfit->SetFillColor(0);
    legfit->SetBorderSize(0);
    legfit->SetFillStyle(0);
    legfit->AddEntry("",Form("B_{c} mass= %1.2f #pm %1.2f MeV",mean.getVal()*1000.0, mean.getError()*1000.0 ),"");
    legfit->Draw();

    //construcción del pull

    pad2->cd();
    RooHist* massPull = mframe->pullHist();
    RooPlot* frame2 = M.frame(Title(" ")) ;
    frame2->addPlotable(massPull,"P") ;
    frame2->SetYTitle(" (Data-Fit)/#sigma");
    frame2->GetXaxis()->CenterTitle();
    frame2->GetYaxis()->CenterTitle();
    frame2->SetXTitle("M(B_{c})[Gev]"); 
    frame2->GetYaxis()->SetNdivisions(505,1);
    frame2->GetXaxis()->SetNdivisions(505,1);
    frame2->GetXaxis()->SetTickLength(0.07);   
    frame2->SetTitleOffset(0.35,"Y");
    frame2->SetTitleSize(0.08,"Y");
    frame2->SetLabelSize(0.06,"XY");
    frame2->SetTitleSize(0.08,"X");
    frame2->Draw();

    c->Draw();
    c->SaveAs("canvas.svg");

}