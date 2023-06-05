using namespace RooFit;

TCanvas* CreateCanvasNom(TString cname, RooFitResult* result, 
RooDataSet *data, RooRealVar M,Double_t supM, Double_t infM , RooRealVar fs,
RooAddPdf MassModel,RooRealVar width,RooRealVar width2, RooRealVar mean)
{
//   Double_t nbin = ((supM-infM)/0.010);
    Double_t nbin = 100;
  
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
//   leg->AddEntry(Mframe->findObject("Sig1")," Gaus1","l");
//   leg->AddEntry(Mframe->findObject("Sig2"),"Gaus2","l");
  leg->Draw();
  
  Double_t Mpsi = mean.getVal();
  Double_t MpsiE = mean.getError();
  

  Double_t G = sqrt( fs.getVal()*width.getVal()*width.getVal() + (1-fs.getVal())*width2.getVal()*width2.getVal() );
  Double_t GE = (1/G)*sqrt( (fs.getVal()*fs.getVal())*(width.getVal()*width.getVal())*(width.getError()*width.getError()) + ((1-fs.getVal())*(1-fs.getVal()))*(width2.getVal()*width2.getVal())*(width2.getError()*width2.getError())  );

  
  TLegend *legpar = new TLegend(0.6,0.68,0.8,0.88);
  legpar->SetTextSize(0.035);
  legpar->SetFillColor(0);
  legpar->SetBorderSize(0);
  legpar->SetFillStyle(0);
  legpar->AddEntry("",Form("M(B_{s}^{0}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
  legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
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

void mcz(Int_t ntotal=1, Int_t seed=4)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //gStyle->SetErrorX(0);

  Double_t Mmin = 40; 
  Double_t Mmax = 140; 

  RooRealVar M("M","",40,140);
  
  
  Double_t supM = Mmax;
  Double_t infM = Mmin;
  
    //FIT
    
    RooDataSet data("data","data",RooArgSet(M));
    RooRealVar mean("mean"," Mass mean",91.8,80,95,"GeV");
    RooRealVar width("width"," Mass width",5,0,10,"GeV");
    RooGaussian sig1("Sig1"," Signal PDF",M,mean,width);

    // RooRealVar mean2("mean2"," Mass mean",70,80,95,"GeV");
    RooRealVar width2("width2"," Mass width",1.4,0,9,"GeV");
    RooGaussian sig2("Sig2"," Signal PDF2",M,mean,width2);

    RooRealVar fs("fsum","fsum",0,1);



    RooAddPdf sumgau("sumgau","sumgau",RooArgList(sig1,sig2),RooArgList(fs));

    RooFitResult* fitres = sumgau.fitTo(data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));

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
  
  Double_t g1,g1E,mu,muE,g2,g2E,Mupull;
  Int_t status,covQual,badNll;
  
  TFile *file = TFile::Open("ToyMC.root","RECREATE");
  TTree *tree1 = new TTree("tree1","tree1");
  file->cd();

 
  tree1->Branch("Mupull",&Mupull);
  


  tree1->Branch("mu",&mu);
  tree1->Branch("muE",&muE);

  tree1->Branch("g1",&g1);
  tree1->Branch("g1E",&g1E);
  
    tree1->Branch("g2",&g2);
    tree1->Branch("g2E",&g2E);



  
  tree1->Branch("status",&status);
  tree1->Branch("covQual",&covQual);
  tree1->Branch("badNll",&badNll);

  //****************************************
  //    Input values from DATA fit in pt bins
  //*****************************************
//   Double_t nsi, nsie, nbi, nbie, ci, cie;
  Double_t mui, muie, g1i, g1ie , g2i, g2ie;
  ifstream entrada ("datos_fit.txt");
  //entrada.is_open();
  if ( !entrada ) 
    {
      cout << "No se pudo abrir el archivo entrada fit pr bins" << endl;
      exit( 1 );
    }
  entrada>> mui >> muie >> g1i >> g1ie >> g2i >> g2ie;
  entrada.close();


  
  Int_t nruns=0;
  
  for(Int_t n=0;n<ntotal;n++)
  {
        
    mean.setVal(mui);
    mean.setError(muie);
    
    width.setVal(g1i);
    width.setError(g1ie);

    width2.setVal(g2i);
    width2.setError(g2ie);

 
    RooDataSet *dataToy = sumgau.generate(RooArgSet(M), Extended(kTRUE));
 
    
    RooFitResult* fitres = sumgau.fitTo(*dataToy,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));
   

    Mupull = (mean.getVal()-mui)/mean.getError();

    g1 = width.getVal();
    g1E = width.getError();

    g2 = width2.getVal();
    g2E = width2.getError();

    mu = mean.getVal();
    muE = mean.getError();
    
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