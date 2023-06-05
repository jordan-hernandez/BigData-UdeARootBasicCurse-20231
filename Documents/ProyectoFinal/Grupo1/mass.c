using namespace RooFit;

void save_resultTfit(ofstream& salida, RooRealVar mean,  RooRealVar width, RooRealVar width2)
{
  salida <<" "  << mean.getVal() << " " << mean.getError() << " "
	 <<  width.getVal() << " " << width.getError() << " " << width2.getVal() << " " << width2.getError()
    ;
  cout << " el archivo se escribio bien" << endl; 
  return;
}
int mass()
{
    // Cargar el archivo
    TFile *f = new TFile("Data.root", "READ");
    f->cd();

    TTree *t = (TTree*)f->Get("Delphes");

    // Crear nuevo archivo y árbol
    TFile *file = TFile::Open("Data_m.root", "RECREATE");
    TTree *tm = new TTree("MuonData", "Muon Data");

    // Variables para lectura de datos
    Float_t Mu_pt[2], Mu_eta[2], Mu_phi[2], Mu_E[2], Mu_iso[2];
    Int_t charge[2], mu_size;

    // Variables para escritura de datos
    Float_t Mu_pt_f[2], Mu_eta_f[2], Mu_phi_f[2], Mu_E_f[2], Mu_iso_f[2];
    Int_t charge_f[2], Mu_size_f,muda;



    // Asociar las variables a las ramas del árbol original
    t->SetBranchAddress("Muon.PT", &Mu_pt);
    t->SetBranchAddress("Muon.Eta", &Mu_eta);
    t->SetBranchAddress("Muon.Phi", &Mu_phi);
    t->SetBranchAddress("Muon.Charge", &charge);
    t->SetBranchAddress("Muon_size", &mu_size);
    t->SetBranchAddress("Muon.IsolationVar", &Mu_iso);

    // Asociar las variables a las ramas del árbol nuevo
    tm->Branch("Muon.PT", &Mu_pt_f,"Muon.PT[2]/F");
    tm->Branch("Muon.Eta", &Mu_eta_f,"Muon.Eta[2]/F");
    tm->Branch("Muon.Phi", &Mu_phi_f,"Muon.Phi[2]/F");
    tm->Branch("Muon.Charge", &charge_f,"Muon.Charge[2]/I");
    tm->Branch("Muon_size", &Mu_size_f,"Muon_size/I");

    

    TLorentzVector muon1;
    TLorentzVector muon2;
    Double_t m_invariant;
    TH1 *h = new TH1F("h", "Invariant mass", 100, 0, 200);


    //--------------------------------------------------------------------------------
    RooRealVar M("M","",40,140);
    RooDataSet data("data","data",RooArgSet(M));


    //--------------------------------------------------------------------------------


    // Loop sobre los muones
    Long64_t nMuons = t->GetEntries();
    for (Long64_t iMuon = 0; iMuon < nMuons; ++iMuon)

    {
        t->GetEntry(iMuon);
        if (mu_size < 2) continue;
        for(int ipart=0; ipart<2; ipart++)
        {
            Mu_pt_f[ipart] = Mu_pt[ipart];
            Mu_eta_f[ipart] = Mu_eta[ipart];
            Mu_phi_f[ipart] = Mu_phi[ipart];
            charge_f[ipart] = charge[ipart];
            Mu_iso_f[ipart] = Mu_iso[ipart];
        }
            
        if (fabs(Mu_pt[0]<30)) continue;
        if (fabs(Mu_pt[1]<30)) continue;
        if (fabs(Mu_eta[0])>2.4) continue;
        if (fabs(Mu_eta[1])>2.4) continue;
        if (Mu_iso[0]>0.1) continue;
        if (Mu_iso[1]>0.1) continue; 
        if (charge[0] == charge[1]) continue;

        tm->Fill();   
        muon1.SetPtEtaPhiE(Mu_pt[0], Mu_eta[0], Mu_phi[0], Mu_E[0]);
        muon2.SetPtEtaPhiE(Mu_pt[1], Mu_eta[1], Mu_phi[1], Mu_E[1]);

        // m_invariant = (muon1 + muon2).M();
        // // cout<<m_invariant<<endl;
        // h->Fill(m_invariant);

        // Obtener las variables de las dos partículas
        Double_t eta1 = muon1.Eta();
        Double_t phi1 = muon1.Phi();
        Double_t pt1 = muon1.Pt();
        Double_t eta2 = muon2.Eta();
        Double_t phi2 = muon2.Phi();
        Double_t pt2 = muon2.Pt();

        // Calcular la masa invariante utilizando la fórmula
        Double_t delta_eta = eta1 - eta2;
        Double_t delta_phi = TMath::Abs(phi1 - phi2);
        if (delta_phi > TMath::Pi())
        {
            delta_phi = 2 * TMath::Pi() - delta_phi;
        }
        Double_t invariant_pt = pt1 * pt2 * (TMath::CosH(delta_eta) - TMath::Cos(delta_phi));
        if (invariant_pt > 0)
        {
        h->Fill(TMath::Sqrt(2 * invariant_pt));
        }
        M = TMath::Sqrt(2 * invariant_pt);    
        data.add(RooArgSet(M));
    }
    // Escribir y cerrar el archivo
    // tm->Scan();
    TCanvas *c = new TCanvas("c", "c", 800, 600);
    h->Fit("gaus");
    h->Draw();
    c->Draw();
    c->Print("invariant_mass.png");
    
    
    
    //// RoFit
    //--------------------

    //FIT
    RooRealVar mean("mean"," Mass mean",91.8,80,95,"GeV");
    RooRealVar width("width"," Mass width",5,0,10,"GeV");
    RooGaussian sig1("Sig1"," Signal PDF",M,mean,width);

    RooRealVar mean2("mean2"," Mass mean",70,80,95,"GeV");
    RooRealVar width2("width2"," Mass width",1.4,0,9,"GeV");
    RooGaussian sig2("Sig2"," Signal PDF2",M,mean,width2);

    RooRealVar fs("fsum","fsum",0,1);
    
    

    RooAddPdf sumgau("sumgau","sumgau",RooArgList(sig1,sig2),RooArgList(fs));

    RooFitResult* fitres = sumgau.fitTo(data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));

    data.Print("v"); 
    fitres->Print("v");


    // TXT

    //parameters output file
    ofstream salida_TotalFit("datos_fit.txt");
    salida_TotalFit.is_open();
    save_resultTfit(salida_TotalFit, mean, width, width2);
    
    
    //PLOT
    Int_t Mmax = 140;
    Int_t Mmin = 40;
    // Double_t nbin = ((Mmax-Mmin)/0.04)+1; 
    Double_t nbin = 50;
    TCanvas* c2 = new TCanvas("c2","c2",800,600);
    TPad *pad1 = new TPad("pad1", "padi",0.01,0.411,0.9903769, 0.99 );
    pad1->SetLeftMargin(0.09);   
    pad1->SetRightMargin(0.019);
    pad1->SetTopMargin(0.09);
    pad1->SetBottomMargin(0.0);  

    TPad *pad2 = new TPad("pad2", "pad2",0.01,0.01,0.9903769,0.41);
    pad2->SetLeftMargin(0.09);
    pad2->SetRightMargin(0.019);  
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.25);
    //pad2->SetTickx(0);
    pad2->SetFillColor(0);
    pad2->SetGridx(0);
    pad2->SetGridy(0);
    
    // RooPlot* frame = M.frame(Title("")) ;
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    RooPlot* frame = M.frame(Mmin,Mmax,nbin) ;
    data.plotOn(frame,DataError(RooAbsData::SumW2),MarkerSize(.75),XErrorSize(0));
    
    sumgau.plotOn(frame,LineColor(kRed),LineWidth(2),Name("Fit")); 
    RooHist* hpull = frame->pullHist() ; 

    sumgau.plotOn(frame,Components(sig1),LineStyle(kDashed)  ,LineColor(kGreen),LineWidth(2),Name("Sig1")); 
    sumgau.plotOn(frame,Components(sig2),LineStyle(kDashed) ,LineColor(kBlue),LineWidth(2),Name("Sig2")); 
    
    
    frame->SetTitle("");
    frame->SetTitleSize(1100);
    frame->SetYTitle("Events");        
    frame->SetLabelSize(0.06,"XY");
    frame->SetTitleSize(0.29,"XY");
    frame->GetXaxis()->SetRangeUser(65,120);
    frame->GetYaxis()->CenterTitle();   
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->SetNdivisions(505,1);
    frame->GetXaxis()->SetNdivisions(505,1);   
    frame->GetXaxis()->SetDecimals(1);
    frame->GetXaxis()->SetTickLength(0.05); 
    frame->SetTitleSize(0.19,"XY");
    frame->SetTitleSize(0.11,"X");
    frame->SetTitleSize(0.1,"Y");  
    frame->SetTitleOffset(0.9,"X");
    frame->SetTitleOffset(0.42,"Y");
    frame->SetMinimum(-40.0);  
    frame->Draw();

    TLegend *leg = new TLegend(0.18,0.58,0.38,0.88); 
    leg->SetTextSize(0.06);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(frame->findObject("Data")," Data","ep"); 
    leg->AddEntry(frame->findObject("Fit")," Fit result","l");
    leg->AddEntry(frame->findObject("Sig1")," Gauss1","l");
    leg->AddEntry(frame->findObject("Sig2"),"Gauss2","l");
    leg->Draw();

    TLatex *   tex1 = new TLatex(0.95,0.93,"L = 61.6 fb^{-1} (#sqrt{s} = 13 TeV)");
    tex1->SetNDC();
    tex1->SetTextAlign(31);
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.053); 
    tex1->SetLineWidth(2);

    tex1->Draw();

    TLatex *   tex2 = new TLatex(0.63,0.93,"Fit to the Z^{0} mass");
    tex2->SetNDC();
    tex2->SetTextAlign(31);
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.073); 
    tex2->SetLineWidth(2);

    tex2->Draw();

    TLatex *   tex3 = new TLatex(0.21,0.93,"CMS");
    tex3->SetNDC();
    tex3->SetTextAlign(31);
    tex3->SetTextFont(61);
    tex3->SetTextSize(0.073);
    tex3->SetLineWidth(2);

    tex3->Draw();

    Double_t Mpsi = mean.getVal();
    Double_t MpsiE = mean.getError();    
    Double_t G = sqrt( fs.getVal()*width.getVal()*width.getVal() + (1-fs.getVal())*width2.getVal()*width2.getVal() );
    Double_t GE = (1/G)*sqrt( (fs.getVal()*fs.getVal())*(width.getVal()*width.getVal())*(width.getError()*width.getError()) + ((1-fs.getVal())*(1-fs.getVal()))*(width2.getVal()*width2.getVal())*(width2.getError()*width2.getError())  );
    
    TLegend *legpar = new TLegend(0.6,0.58,0.8,0.88);
    legpar->SetTextSize(0.06);
    legpar->SetTextFont(42);
    legpar->SetFillColor(0);
    legpar->SetBorderSize(0);
    legpar->SetFillStyle(0);
    legpar->AddEntry("",Form("M(Z^{0}) = %1.2f #pm %1.2f GeV",Mpsi,MpsiE),"");
    legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f GeV",G,GE),"");
    legpar->Draw();




    pad2->cd();
       
    RooPlot* frame2 = M.frame(Mmin,Mmax,nbin) ; 
    frame2->addPlotable(hpull,"P") ;
    frame2->GetYaxis()->SetRangeUser(-4.5,4.5);
    frame2->GetXaxis()->SetRangeUser(65,120);
    frame2->GetYaxis()->SetNdivisions(5);
    // frame2->GetXaxis()->SetNdivisions(5);
    frame2->SetYTitle(" (Data-Fit)/#sigma");
    frame2->SetLabelSize(0.1,"XY");
    frame2->SetTitleSize(0.11,"X");
    frame2->SetTitleSize(0.135,"Y");
    // frame2->GetYaxis()->SetNdivisions(505,1);
    frame2->GetXaxis()->SetNdivisions(505,3);
  
    frame2->GetYaxis()->SetLabelSize(0.08);
    frame2->GetXaxis()->SetLabelSize(0.08);
    frame2->SetTitleOffset(0.31,"Y");
    frame2->GetXaxis()->SetTitle("M_{#mu#mu} (GeV)");
    frame2->GetXaxis()->SetTitleSize(0.13);
    frame2->SetTitle("");
    frame2->GetYaxis()->CenterTitle();   
    frame2->GetXaxis()->CenterTitle();
    frame2->SetTitleOffset(0.9,"X");
    frame2->SetTitleOffset(0.31,"Y");
    frame2->GetXaxis()->SetTickLength(0.06); 
    frame2->Draw();


    TLine *line1 = new TLine(65,0.0,120,0.0);
    line1->SetLineColor(1);
    //line1->SetLineStyle(2);
    line1->SetLineWidth(1);
    line1->Draw();
    c2->Update();
    c2->Draw();
    c2->Print("invariant_mass_fit.png");
    //--------------------
    tm->Print();
    file->Write();
    file->Close();
    f->Close();

    return 0;
}
