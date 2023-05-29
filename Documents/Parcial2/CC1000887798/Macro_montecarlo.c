using namespace RooFit;

TCanvas* CreateCanvasNom(TString cname, RooFitResult* result, 
RooDataSet *data, RooRealVar M,Double_t supM, Double_t infM , 
RooAddPdf MassModel, RooAddPdf sumgau, RooChebychev bkg1, 
RooRealVar Ns, RooRealVar Nb, RooRealVar width, RooRealVar mean)
{
  Double_t nbin = ((supM-infM)/0.010);
  
  int H = 600;
  int W = 800;
  TCanvas *c1 = new TCanvas("c1","",50,50,W,H);
  c1->cd();
  c1->SetLeftMargin(0.13);
  c1->SetRightMargin(0.02);
  c1->SetTopMargin(0.09);
  c1->SetBottomMargin(0.13); 
  
  RooPlot* Mframe = M.frame(infM,supM,nbin);
  data->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0));
  MassModel.plotOn(Mframe,DrawOption("F"),FillColor(0),LineWidth(2),Name("fittotal"));
  MassModel.plotOn(Mframe,Components(sumgau),LineColor(kRed),LineWidth(2),Name("Signal")); 
  MassModel.plotOn(Mframe,Components(bkg1),LineColor(kGreen),LineWidth(2),Name("bkg")); 
  data->plotOn(Mframe,DataError(RooAbsData::SumW2),MarkerSize(1.0),XErrorSize(0),Name("Data"));
  MassModel.plotOn(Mframe);
  
  Mframe->SetYTitle("Events / 10 MeV"); 
  Mframe->SetLabelSize(0.04,"XY");
  Mframe->SetTitleSize(0.05,"XY");
  Mframe->GetYaxis()->CenterTitle();   
  Mframe->GetXaxis()->CenterTitle();
  Mframe->GetYaxis()->SetNdivisions(505,1);
  Mframe->GetXaxis()->SetNdivisions(505,1);   
  Mframe->GetXaxis()->SetDecimals(1); 
  Mframe->SetTitleOffset(0.9,"X");
  Mframe->SetTitleOffset(1.1,"Y");
  Mframe->SetTitleSize(0.06,"XY");
  Mframe->SetMinimum(1.0); 
  Mframe->Draw();
  
  //TLegend *leg = new TLegend(0.15,0.68,0.35,0.88);
  TLegend *leg = new TLegend(0.18,0.68,0.38,0.88); 
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(Mframe->findObject("Data")," Data","ep"); 
  leg->AddEntry(Mframe->findObject("fittotal")," Fit result","l");
  leg->AddEntry(Mframe->findObject("Signal")," Signal","l");
  leg->AddEntry(Mframe->findObject("bkg"),"Comb. backg.","l");
  leg->Draw();
  
  Double_t Mpsi = mean.getVal()*1000.0;
  Double_t MpsiE = mean.getError()*1000.0;
  Double_t G = sqrt(width.getVal()*width.getVal())*1000.0;
  Double_t GE = (1/G)*sqrt( (width.getVal()*width.getVal())*(width.getError()*width.getError())*1000.0);  //OJO
  
  TLegend *legpar = new TLegend(0.6,0.68,0.8,0.88);
  legpar->SetTextSize(0.035);
  legpar->SetFillColor(0);
  legpar->SetBorderSize(0);
  legpar->SetFillStyle(0);
  legpar->AddEntry("",Form("M(B_{s}^{0}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
  legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
  legpar->AddEntry("",Form("N_{B_{s}^{0}} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
  legpar->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
  legpar->Draw();
  
  TLegend *legMass = new TLegend(0.64,0.57,0.83,0.65);
  legMass->SetTextFont(43); 
  legMass->SetTextSize(20);  
  legMass->SetFillColor(0); 
  legMass->SetBorderSize(0);
  legMass->SetFillStyle(0); 
  legMass->Draw(); 
  

  TLatex *   tex1 = new TLatex(0.92,0.926,"ToyMC");
  tex1->SetNDC();
  tex1->SetTextAlign(31);
  tex1->SetTextFont(42);
  tex1->SetTextSize(0.05); 
  tex1->SetLineWidth(2);
  
  TLatex *tex2 = new TLatex(0.2,0.926,"CMS");
  tex2->SetNDC();
  tex2->SetTextFont(61);
  tex2->SetTextSize(0.05); 
  tex2->SetLineWidth(2);
  
  TLatex *tex3 = new TLatex(0.29,0.926,"Preliminary");
  tex3->SetNDC();
  tex3->SetTextFont(52);
  tex3->SetTextSize(0.05); 
  tex3->SetLineWidth(2);
  
  tex1->Draw();  
  tex2->Draw();
  tex3->Draw();
  
  c1->Modified();
  gPad->Update();
  gPad->RedrawAxis();
  TLine l;
  l.DrawLine(gPad->GetUxmax(), gPad->GetUymax(), gPad->GetUxmax(), gPad->GetUymin());
  return c1; 
 
}

void Macro_montecarlo(Int_t ntotal=1, Int_t seed=5)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //gStyle->SetErrorX(0);

  Double_t Mmin = 6.05; 
  Double_t Mmax = 6.5; 

  RooRealVar M("M","",Mmin,Mmax);
  
  
  Double_t supM = Mmax;
  Double_t infM = Mmin;
  
RooRealVar mean("mean"," Mass mean",6.2751,6.25,6.3,"GeV");
RooRealVar width("width"," Mass width",0.015,0.01,0.04,"GeV");
RooGaussian Sig("Sig"," Signal PDF",M,mean,width);

RooRealVar c("c","c",-10.0,10.0);
RooChebychev  bkg1("bkg1","Exp. Background",M,c);

RooRealVar Ns("Ns","Ns",0.,2000);
RooRealVar Nb("Nb","Nb",0.,2000);   

RooAddPdf MassModel("MassModel","MassModel",RooArgList(Sig,bkg1),RooArgList(Ns,Nb));

RooRandom::randomGenerator()->SetSeed(seed);


RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
RooMsgService::instance().setSilentMode(kTRUE);
RooMsgService::instance().setStreamStatus(1,false);
RooMsgService::instance().getStream(1).removeTopic(Integration) ;  
RooMsgService::instance().getStream(1).removeTopic(Minimization) ;  
RooMsgService::instance().getStream(1).removeTopic(Fitting) ;  
RooMsgService::instance().getStream(1).removeTopic(NumIntegration) ;
RooMsgService::instance().getStream(1).removeTopic(Optimization) ;
RooMsgService::instance().getStream(1).removeTopic(ObjectHandling) ;
RooMsgService::instance().getStream(1).removeTopic(Eval) ;
RooMsgService::instance().Print() ;
  
  Double_t Nsig,NsigE,Nbkg,NbkgE,g1,g1E,mu,muE,exp,expE, Nspull, Nsdif, Nbpull, Mupull;
  Int_t status,covQual,badNll;
  
  TFile *file = TFile::Open("ToyMC.root","RECREATE");
  TTree *tree1 = new TTree("tree1","tree1");
  file->cd();

  tree1->Branch("Nspull",&Nspull);
  tree1->Branch("Nsdif",&Nsdif);
  tree1->Branch("Nbpull",&Nbpull);    
  tree1->Branch("Mupull",&Mupull);
  
  tree1->Branch("Nsig",&Nsig);
  tree1->Branch("NsigE",&NsigE);
  
  tree1->Branch("Nbkg",&Nbkg);
  tree1->Branch("NbkgE",&NbkgE);

  tree1->Branch("mu",&mu);
  tree1->Branch("muE",&muE);

  tree1->Branch("g1",&g1);
  tree1->Branch("g1E",&g1E);
  
  tree1->Branch("exp",&exp);
  tree1->Branch("expE",&expE);

  
  tree1->Branch("status",&status);
  tree1->Branch("covQual",&covQual);
  tree1->Branch("badNll",&badNll);

  //****************************************
  //    Input values from DATA fit in pt bins
  //*****************************************
  Double_t nsi, nsie, nbi, nbie, ci, cie;
  Double_t mui, muie, g1i, g1ie ;
  ifstream entrada ("datos_fit.txt");
  //entrada.is_open();
  if ( !entrada ) 
    {
      cout << "No se pudo abrir el archivo entrada fit pr bins" << endl;
      exit( 1 );
    }
  entrada>> nsi >> nsie >> nbi >> nbie >> ci >> cie >> mui >> muie >> g1i >> g1ie  ;
  entrada.close();
  cout << "Ns, Nb, g2 values: " <<  nsi << " " << nbi  << endl;
  //return;
  
  Int_t nruns=0;
  
  for(Int_t n=0;n<ntotal;n++)
  {
    
    Ns.setVal(nsi);
    Ns.setError(nsie);
    
    Nb.setVal(nbi);
    Nb.setError(nbie);

    c.setVal(ci);
    c.setError(cie);
        
    mean.setVal(mui);
    mean.setError(muie);
    
    width.setVal(g1i);
    width.setError(g1ie);

 
    RooDataSet *dataToy = MassModel.generate(RooArgSet(M), Extended(kTRUE));
 
    
    RooFitResult* fitres = MassModel.fitTo(*dataToy,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));
   

    Nspull = (Ns.getVal()-nsi)/Ns.getError();
    Nsdif = (Ns.getVal()-nsi);
    Nbpull = (Nb.getVal()-nbi)/Nb.getError();   
    Mupull = (mean.getVal()-mui)/mean.getError();

    Nsig = Ns.getVal();
    NsigE = Ns.getError();
    
    Nbkg = Nb.getVal();
    NbkgE = Nb.getError();

    g1 = width.getVal();
    g1E = width.getError();

    mu = mean.getVal();
    muE = mean.getError();


    exp = c.getVal();
    expE = c.getError();
    
    status = fitres->status();
    covQual = fitres->covQual();
    badNll = fitres->numInvalidNLL();


    tree1->Fill();

    nruns++;
    /*
    TCanvas* canv_nominal = CreateCanvasNom("canv_nominal", fitres, dataToy, M, supM, infM, MassModel, sumgau, bkg1, Ns, Nb, width, width2, fs, mean,ptl,pth); 
    canv_nominal->Print(Form("plots/mass_BsFit_ptbins_ToyMC_%1.0f_%1.0f.png",ptl,pth));
    canv_nominal->Print(Form("plots/mass_BsFit_ptbins_ToyMC_%1.0f_%1.0f.pdf",ptl,pth));
    */
    
    delete dataToy;
    delete fitres;

  }

 tree1->Write();
 //file->Write("",TObject::kOverwrite);
 
 cout<<" for end: "<<nruns<<endl;

}//End analysis
