
using namespace std;
using namespace RooFit;

void save_resultTfit(ofstream& salida, RooRealVar Ns, RooRealVar Nb, RooRealVar a, RooRealVar mean,  RooRealVar width)
{
  salida <<  Ns.getVal() << " " << Ns.getError() << " " << Nb.getVal() << " " << Nb.getError() <<" " <<  a.getVal() << " " << a.getError() << " "  
	 <<" "<< mean.getVal() << " " << mean.getError() << " "
	 <<  width.getVal() << " " << width.getError()
    ;
  cout << " el archivo se escribio bien" << endl; 
  return;
}


int leerarbol(){
  
  Float_t mass, time, error;

  TFile f("masayT.root", "READ");
  TTree *tree = (TTree*)f.Get("tree");

  tree->SetBranchAddress("mass", &mass);
  tree->SetBranchAddress("time", &time);
  tree->SetBranchAddress("error", &error);


  //Masa del Bc: 6.2751GeV
  Double_t Mmin = 6.05; 
  Double_t Mmax = 6.5; 

  RooRealVar M("M","M(B_{c}^{+})",Mmin,Mmax);
  RooDataSet data("data","data",RooArgSet(M));

  Long64_t nentries = tree->GetEntries();
  cout<<" Entries : "<<nentries<<endl;

  for (int evt=0; evt < nentries; evt++) 

  {
      tree->GetEvent (evt);

      if (error==0) continue;
      if (isnan(error)==1) continue;
      if((time/error)<5.0) continue;

      M=mass;   
      data.add(RooArgSet(M));
  }
  data.Print("v");

  // tree->ResetBranchAddresses(); //

  RooRealVar mean("mean"," Mass mean",6.2751,6.25,6.3,"GeV");
  RooRealVar width("width"," Mass width",0.015,0.01,0.04,"GeV");
  RooGaussian gauss("Sig"," Signal PDF",M,mean,width);


  RooRealVar a("a","a",-10.0,10.0);
  RooChebychev bkg("bkg","Background",M,a);

  RooRealVar Ns("Ns","Ns",0.,2000);
  RooRealVar Nb("Nb","Nb",0.,2000); 

  RooAddPdf MassModel("MassModel","MassModel",RooArgList(gauss,bkg),RooArgList(Ns,Nb));

  RooFitResult* fitres = MassModel.fitTo(data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));

  data.Print("v");
  fitres->Print("v");


  // TCanvas* canv_nominalpull = CreateCanvasNomPull("canv_nominalpull", fitres, data, M, Mmin, Mmax, MassModel, gauss, bkg, Ns, Nb, width, mean);  
  // canv_nominalpull->Print(Form("Hijueputa.png"));


  ofstream salida_TotalFit(Form("output_TotalFit_xd.txt"));
  salida_TotalFit.is_open();
  save_resultTfit(salida_TotalFit, Ns, Nb, a, mean, width);




//************************************************************************************************************
//************************************************************************************************************
//************************************************************************************************************
//************************************************************************************************************
//************************************************************************************************************

//GRÁFICA

   Double_t nbin = ((Mmax-Mmin)/0.008)+1;

  int H = 600;
  int W = 800;

  TCanvas *c1 = new TCanvas("c","c",50,50,W,H);
  //TCanvas *c1 = new TCanvas(cname,cname,W,H);
  //c1->Divide(1,2);
  //c1->cd(1) ;
  c1->cd() ;  
  c1->SetLeftMargin(0.005);
  c1->SetRightMargin(0.01);
  c1->SetTopMargin(0.09);
  c1->SetBottomMargin(0.1);

  TPad *pad1 = new TPad("pad1", "padi",0.01,0.411,0.9903769, 0.99 );
  pad1->SetLeftMargin(0.09);   
  pad1->SetRightMargin(0.019);
  pad1->SetTopMargin(0.09);
  pad1->SetBottomMargin(0.013);  

  TPad *pad2 = new TPad("pad2", "pad2",0.01,0.01,0.9903769,0.41);
  pad2->SetLeftMargin(0.09);
  pad2->SetRightMargin(0.019);  
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.25);
  //pad2->SetTickx(0);
  pad2->SetFillColor(0);
  pad2->SetGridx(0);
  pad2->SetGridy(0);

  pad1->Draw();
  pad2->Draw();
  pad1->cd(); 

  RooPlot* Mframe = M.frame(Mmin,Mmax,nbin);
  data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0));
  //MassModel.plotOn(Mframe);
  MassModel.plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));
  
  RooHist* hpullm2 = Mframe->pullHist() ;
  
    // MassModel.plotOn(Mframe,Components(bkg),DrawOption("F"),FillColor(kCyan),LineColor(kCyan),LineWidth(2),Name("bkg"));

  MassModel.plotOn(Mframe,Components(gauss),LineColor(kRed),LineWidth(2),LineStyle(kDashed),Name("Signal")); 
  MassModel.plotOn(Mframe,Components(bkg),LineColor(kGreen),LineWidth(2), LineStyle(kDashed),Name("bkg")); 
  // MassModel.plotOn(Mframe,DrawOption("F"),FillColor(kBlue-9),LineColor(1),LineWidth(1),Name("fittotal"));
  // MassModel.plotOn(Mframe,Components(gauss),DrawOption("F"),FillColor(kBlue-9),LineColor(kBlue-9),LineWidth(2),Name("Signal")); 
  
  data.plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("Data"));
  MassModel.plotOn(Mframe);
  
  Mframe->SetYTitle("Events / 10 MeV");
  // Mframe->SetXTitle("M(B_{c}^{+})");  
  Mframe->SetTitle("");  
  Mframe->SetLabelSize(0.05,"XY");
  Mframe->SetTitleSize(0.08,"X");
  Mframe->SetTitleSize(0.05,"Y"); 
  Mframe->GetYaxis()->CenterTitle();   
  // Mframe->GetXaxis()->CenterTitle();
  Mframe->GetYaxis()->SetNdivisions(505,1);
  // Mframe->GetXaxis()->SetNdivisions(505,1);
  Mframe->GetXaxis()->SetTickLength(0.06);    
  Mframe->GetXaxis()->SetDecimals(1); 
  Mframe->SetTitleOffset(0.8,"X");
  Mframe->SetTitleOffset(0.65,"Y");
  Mframe->SetMinimum(-1.0); 
  Mframe->Draw();
  
  //TLegend *leg = new TLegend(0.15,0.68,0.35,0.88);
  TLegend *leg = new TLegend(0.18,0.58,0.38,0.88); 
  leg->SetTextSize(0.06);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(Mframe->findObject("Data")," Data","ep"); 
  leg->AddEntry(Mframe->findObject("fittotal")," Fit result","l");
  leg->AddEntry(Mframe->findObject("Signal")," Signal","l");
  leg->AddEntry(Mframe->findObject("bkg"),"Comb. backg.","l");
  //leg->AddEntry(Mframe->findObject("Signal")," Signal","f");
  //leg->AddEntry(Mframe->findObject("bkg"),"Comb. backg.","f");
  leg->Draw();
  
  TLegend *legpar = new TLegend(0.6,0.58,0.8,0.88);
  legpar->SetTextSize(0.06);
  legpar->SetTextFont(42);
  legpar->SetFillColor(0);
  legpar->SetBorderSize(0);
  legpar->SetFillStyle(0);
  legpar->AddEntry("",Form("M(B_{c}^{+})= %1.3f #pm %1.3f MeV",mean.getVal()*1000,mean.getError()*1000),"");
  legpar->AddEntry("",Form("#sigma = %1.3f #pm %1.3f MeV",width.getVal()*1000,width.getError()*1000),"");
  legpar->AddEntry("",Form("N_{B_{c}^{+}} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
  legpar->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
  legpar->Draw();
  
 
  pad2->cd();
  
  // Create a new frame to draw the pull distribution 
  RooPlot* framem2 = M.frame(Mmin,Mmax,nbin) ;
  framem2->addPlotable(hpullm2,"P") ;
  framem2->SetTitle("");
  framem2->SetTitleSize(0.1);
  framem2->SetYTitle(" (Data-Fit)/#sigma");
  framem2->SetXTitle("M(J/#psi#pi^{+}) [GeV]");  
  framem2->SetLabelSize(0.08,"XY");
  framem2->SetTitleSize(0.08,"X");
  framem2->SetTitleSize(0.08,"Y");  
  framem2->GetYaxis()->CenterTitle();   
  framem2->GetXaxis()->CenterTitle();
  framem2->GetYaxis()->SetNdivisions(505,1);
  // framem2->GetXaxis()->SetNdivisions(505,1);
  framem2->GetXaxis()->SetTickLength(0.07);   
  framem2->SetTitleOffset(0.9,"X");
  framem2->SetTitleOffset(0.4,"Y");
  framem2->SetMaximum(4.9);
  framem2->SetMinimum(-4.9);
  
  framem2->Draw();

  TLine *line1 = new TLine(Mmin,0.0,Mmax,0.0);
  line1->SetLineColor(2);
  //line1->SetLineStyle(2);
  line1->SetLineWidth(1);
  line1->Draw();
  
  c1->cd();

  TLatex *   tex1 = new TLatex(0.88,0.95,"L = 19.7 fb^{-1} (#sqrt{s} = 8 TeV)");
  tex1->SetNDC();
  tex1->SetTextAlign(31);
  tex1->SetTextFont(42);
  tex1->SetTextSize(0.04); 
  tex1->SetLineWidth(2);
  
  TLatex *tex2 = new TLatex(0.15,0.95,"CMS");
  tex2->SetNDC();
  tex2->SetTextFont(61);
  tex2->SetTextSize(0.04); 
  tex2->SetLineWidth(2);
  
  TLatex *tex3 = new TLatex(0.22,0.95,"Preliminary");
  tex3->SetNDC();
  tex3->SetTextFont(52);
  tex3->SetTextSize(0.04); 
  tex3->SetLineWidth(2);
  

  tex1->Draw();  
  tex2->Draw();
  tex3->Draw();
  
  c1->Modified();
  gPad->Update();
  gPad->RedrawAxis();
  TLine l;
  l.DrawLine(gPad->GetUxmax(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymin());
  c1->Print(Form("masafit.png"));


  return 0;
}
