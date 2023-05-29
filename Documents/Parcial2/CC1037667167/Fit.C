
{   using namespace RooFit;
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
        
    // Construimos nuestro objeto DataSet para la masa
    RooRealVar Mass("Mass", "Mass", 6.05, 6.5);
    RooDataSet Data("Data", "Data", tree, RooArgSet(Mass));

    Data.Print();
    /* Construccion del modelo de masa---------------------------------------------------------------------------------*/
    RooRealVar Mean("Mean", "Mean", 6, 6.6);
    RooRealVar Sigma("Sigma", "Sigma", 0.00001, 1);
    RooGaussian Signal("Signal", "Signal", Mass, Mean, Sigma);

    RooRealVar C("C", "C", -10, 10);
    RooExponential Backg("Backg", "Backg", Mass, C);

    RooRealVar Nsig("Nsig", "Nsig", 0, 13000);
    RooRealVar Nbkg("Nbkg", "Nbkg", 0, 13000);
    RooAddPdf Model("Model", "Model", RooArgList(Signal, Backg), RooArgList(Nsig, Nbkg));
    RooFitResult* FitResult = Model.fitTo(Data, Extended(kTRUE), Save(kTRUE));

    if (!FitResult) {
        std::cout << "Error: No se pudo fitear el modelo !!" << std::endl;
        file.Close();
        exit(1);
    }

    // Obtener los parametros del ajuste
    double massVal = Mean.getVal();
    double massErr = Mean.getError();

    double nsigVal = Nsig.getVal();
    double nsigErr = Nsig.getError();

    double nbkgVal = Nbkg.getVal();
    double nbkgErr = Nbkg.getError();

    Int_t nbin = 100; //numero de bins para los graficos.
    /*Graficación------------------------------------------------------------------------------------------------------*/
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas->Divide(1, 2 , 0 ,0);

    canvas->cd(1); gPad->SetRightMargin(0.01);
    RooPlot* frame = Mass.frame(Title(" "));
    Data.plotOn(frame);
    Model.plotOn(frame);
    frame->Draw("A");
    frame->GetYaxis()->SetTitle((std::string(frame->GetYaxis()->GetTitle()) + "[GeV]").c_str());
    frame->GetYaxis()->CenterTitle(); 
    frame->GetYaxis()->SetTitleSize(0.07); 
    frame->GetYaxis()->SetTitleOffset(0.5);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->SetXTitle(Form("M(J/#psi #pi^{+})[Gev]"));
    
    TLatex* latex = new TLatex();
    latex->SetNDC();
    //latex->SetTextFont(42);
    latex->SetTextSize(0.06);
    latex->DrawLatex(0.135, 0.9, "B_{c}^{+} #rightarrow  J/#psi #pi^{+}");
    latex->DrawLatex(0.65, 0.9, Form("#chi^{2} = %.2f" ,frame->chiSquare()));
    latex->DrawLatex(0.65, 0.8, Form("#mu = (%.4f #pm %.4f) [GeV]" ,massVal , massErr));
    latex->DrawLatex(0.65, 0.7, Form("Nsig = %.0f #pm %.0f" ,nsigVal , nsigErr));
    latex->DrawLatex(0.65, 0.6, Form("Nbkg = %.0f #pm %.0f" ,nbkgVal , nbkgErr));
    canvas->Update();

    canvas->cd(2);gPad->SetRightMargin(0.01); gPad->SetBottomMargin(0.3);
    RooHist* pullHist = frame->pullHist();
    RooPlot* pullFrame = Mass.frame(Title(" "));
    pullFrame->addPlotable(pullHist, "P");
    pullFrame->GetYaxis()->SetNdivisions(505);
    pullFrame->GetYaxis()->SetTitle("(Data - Fit) /#sigma");
    pullFrame->GetYaxis()->CenterTitle(); 
    pullFrame->GetYaxis()->SetTitleSize(0.07); 
    pullFrame->GetYaxis()->SetTitleOffset(0.5);
    pullFrame->GetYaxis()->SetLabelSize(0.045);
    pullFrame->GetXaxis()->SetTitle(Form("M(J/#psi #pi^{+})[Gev]"));
    pullFrame->GetXaxis()->CenterTitle(); 
    pullFrame->GetXaxis()->SetTitleSize(0.07); 
    pullFrame->GetXaxis()->SetTitleOffset(0.9);
    pullFrame->GetXaxis()->SetLabelSize(0.045);

    TLine* zeroLine = new TLine(6.05, 0, 6.5, 0);
    zeroLine->SetLineStyle(2);
    pullFrame->Draw("A");
    zeroLine->Draw("same");

    canvas->cd(1);
    Model.plotOn(frame,Components(Signal),LineColor(6),LineWidth(2),LineStyle(kDashed),Name("Signal")); 
    Model.plotOn(frame,Components(Backg),LineColor(8),LineWidth(2),LineStyle(kDotted),Name("Bkg"));
    frame->Draw("Same");
    canvas->Update();

    TLegend *legend = new TLegend(0.15, 0.1, 0.35, 0.4);
    legend->AddEntry(frame->getObject(0), "Data", "PES");
    legend->AddEntry(frame->getObject(1), "Fit", "l");
    legend->AddEntry(frame->getObject(2), "B_{c}^{+} Signal", "l");
    legend->AddEntry(frame->getObject(3), "Comb. background", "l");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetMargin(0.1);
    legend->Draw();
    canvas->Update();
    
    canvas->SaveAs("fit_result.png");
    //canvas->SaveAs("fit_result.root");

    // Limpiar memoria
    delete canvas;
    delete legend;
    delete latex;
    delete zeroLine;
    delete FitResult;

    /*Estudio Mc--------------------------------------------------------------------------------------------------------*/
    gStyle->SetOptStat(0);
    TCanvas *MC_canvas = new TCanvas("Estudio MC", "Estudio MC", 900, 900);
    MC_canvas->Divide(3, 2);

    RooMCStudy *MC = new RooMCStudy(Model, Mass, Binned(false), Silence(true), Extended(true), FitOptions(Save(true), PrintEvalErrors(0)));
    MC->generateAndFit(10000); //---- FIT DE 10 HORAS!!----
    //MC->generateAndFit(10); //---Fit "Rapido"---

    RooPlot *MassframeMean = MC->plotParam(Mean, Bins(nbin));
    RooPlot *MassframeMeanE = MC->plotError(Mean, Bins(nbin));
    RooPlot *MassframeMeanPull = MC->plotPull(Mean, Bins(nbin), FitGauss(true));

    RooPlot *NsigframeMean = MC->plotParam(Nsig, Bins(nbin));
    RooPlot *NsigframeMeanE = MC->plotError(Nsig, Bins(nbin));
    RooPlot *NsigframeMeanPull = MC->plotPull(Nsig, Bins(nbin), FitGauss(true));

    RooPlot *NbkgframeMean = MC->plotParam(Nbkg, Bins(nbin));
    RooPlot *NbkgframeMeanE = MC->plotError(Nbkg, Bins(nbin));
    RooPlot *NbkgframeMeanPull = MC->plotPull(Nbkg, Bins(nbin), FitGauss(true));

    MC_canvas->cd(1); gPad->SetLeftMargin(0.1);
    MassframeMean->SetTitle("Toy Mc Mass");
    MassframeMean->GetYaxis()->SetTitleOffset(1.4);
    MassframeMean->Draw();

    MC_canvas->cd(2); gPad->SetLeftMargin(0.1);
    NsigframeMean->SetTitle("Toy Mc Nsig");
    NsigframeMean->GetYaxis()->SetTitleOffset(1.4);
    NsigframeMean->Draw();

    MC_canvas->cd(3); gPad->SetLeftMargin(0.1);
    NbkgframeMean->SetTitle("Toy Mc Nbkg");
    NbkgframeMean->GetYaxis()->SetTitleOffset(1.4);
    NbkgframeMean->Draw();

    MC_canvas->cd(4); gPad->SetLeftMargin(0.1);
    MassframeMeanPull->SetTitle("Toy Mc Mass Pull");
    MassframeMeanPull->GetYaxis()->SetTitleOffset(1.4);
    MassframeMeanPull->Draw();

    MC_canvas->cd(5); gPad->SetLeftMargin(0.1);
    NsigframeMeanPull->SetTitle("Toy Mc Nsig Pull");
    NsigframeMeanPull->GetYaxis()->SetTitleOffset(1.4);
    NsigframeMeanPull->Draw();

    MC_canvas->cd(6); gPad->SetLeftMargin(0.1);
    NbkgframeMeanPull->SetTitle("Toy Mc Nbkg Pull");
    NbkgframeMeanPull->GetYaxis()->SetTitleOffset(1.4);
    NbkgframeMeanPull->Draw();

    MC_canvas->SaveAs("ToyMcPlot.png");
    //MC_canvas->SaveAs("ToyMcPlot.root");

    // Liberar memoria
    delete MC;
    delete MassframeMean;
    delete MassframeMeanE;
    delete MassframeMeanPull;
    delete NsigframeMean;
    delete NsigframeMeanE;
    delete NsigframeMeanPull;
    delete NbkgframeMean;
    delete NbkgframeMeanE;
    delete NbkgframeMeanPull;
    delete MC_canvas;

    file.Close();
}


