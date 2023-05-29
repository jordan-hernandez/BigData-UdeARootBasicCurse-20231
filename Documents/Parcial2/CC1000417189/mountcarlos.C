using namespace RooFit;
using namespace std;

void mountcarlos(Int_t ntotal=1000, Int_t seed=5)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //gStyle->SetErrorX(0);

Double_t Mmin = 6.05; 
Double_t Mmax = 6.5; 

RooRealVar M("M"," M",Mmin,Mmax);
  
RooRealVar mean("mean"," Mass mean",6.2751,6.25,6.3,"GeV");
RooRealVar width("width"," Mass width",0.015,0.01,0.04,"GeV");
RooGaussian gauss("Sig"," Signal PDF",M,mean,width);

// RooRealVar width2("width2"," Mass width2 ",0.020,0.015,0.05,"GeV");
// RooGaussian Sig2("Sig2"," Signal PDF B",M,mean,width2);


RooRealVar a0("a0","a0",-10.0,10.0);
RooChebychev bkg("bkg","Background",M,a0);

RooRealVar Ns("Ns","Ns",0.,2000);
RooRealVar Nb("Nb","Nb",0.,2000); 

RooAddPdf MassModel("MassModel","MassModel",RooArgList(gauss,bkg),RooArgList(Ns,Nb));

  //------------ Fit procedure -------------------
  
  //RooRandom::randomGenerator()->SetSeed(3);
  RooRandom::randomGenerator()->SetSeed(seed);
  // esto que sigue solo es para que no imprima en pantalla todo lo que hace del ajuste
  
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
  
  Double_t llh,edm,Nsig,NsigE,Nbkg, NbkgE, a, aE, mu, muE, sigma, sigmaE;
  Double_t  Nspull, Nsdif, Nbpull, Mupull;
  Int_t status,covQual,badNll;
  //TFile *file = new TFile("ToyMC_Bs.root","RECREATE");
  //TFile *file = new TFile(Form("ToyMC_Bs_ptbins_%1i_%1.0f_%1.0f.root",seed,ptl,pth),"RECREATE");
  TFile *file = TFile::Open(Form("ToyMC_Bc_bins.root"),"RECREATE");
  TTree *tree = new TTree("tree","tree");
  file->cd();

  tree->Branch("Nspull",&Nspull);
  tree->Branch("Nsdif",&Nsdif);
  tree->Branch("Nbpull",&Nbpull);    
  tree->Branch("Mupull",&Mupull);
  
  tree->Branch("Nsig",&Nsig);
  tree->Branch("NsigE",&NsigE);
  
  tree->Branch("Nbkg",&Nbkg);
  tree->Branch("NbkgE",&NbkgE);

  tree->Branch("mu",&mu);
  tree->Branch("muE",&muE);

  tree->Branch("sigma",&sigma);
  tree->Branch("sigmaE",&sigmaE);

  tree->Branch("a",&a);
  tree->Branch("aE",&aE);

  tree->Branch("edm",&edm);
  tree->Branch("llh",&llh);
  tree->Branch("status",&status);
  tree->Branch("covQual",&covQual);
  tree->Branch("badNll",&badNll);

  //****************************************
  //    Input values from DATA fit in pt bins
  //*****************************************
  Double_t nsi, nsie, nbi, nbie, ai, aie, devi,devie;
  Double_t mui, muie;
 
  ifstream entrada (Form("output_TotalFit_xd.txt"));
  //entrada.is_open();
  if ( !entrada ) 
    {
      cout << "No se pudo abrir el archivo entrada fit pr bins" << endl;
      exit( 1 );
    }
    // Ns, Nb, a, mean, width
  entrada>> nsi >> nsie >> nbi >> nbie >> ai >> aie >> mui >> muie >> devi >> devie;
  entrada.close();
  cout << "Ns, Nb values: " <<  nsi << " " << nbi << endl;
  //return;
  
  Int_t nruns=0;
  
  for(Int_t n=0;n<ntotal;n++)
  {
    //Ns.setConstant(kTRUE);
    Ns.setVal(nsi);
    Ns.setError(nsie);
    
    Nb.setVal(nbi);
    Nb.setError(nbie);

    a0.setVal(ai);
    a0.setError(aie);
        
    mean.setVal(mui);
    mean.setError(muie);
    
    width.setVal(devi);
    width.setError(devie);

    //RooRandom::randomGenerator()->SetSeed(n);
    RooDataSet *dataToy = MassModel.generate(RooArgSet(M), Extended(kTRUE));
    //RooDataSet *dataToy = MassModel.generate(RooArgSet(M),Ns.getVal() +  Nb.getVal() );
    //RooDataSet dataToy = MassModel.generate(RooArgSet(M), Ns.getVal() + Nb.getVal() );
    
    RooFitResult* fitres = MassModel.fitTo(*dataToy,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(6));
    //RooFitResult* fitres = MassModel.fitTo(data,Extended(),Minos(kTRUE),Save(kTRUE), NumCPU(4));
    //dataToy->Print("v");
    //fitres->Print("v");

    Nspull = (Ns.getVal()-nsi)/Ns.getError();
    Nsdif = (Ns.getVal()-nsi);
    Nbpull = (Nb.getVal()-nbi)/Nb.getError();   
    Mupull = (mean.getVal()-mui)/mean.getError();

    Nsig = Ns.getVal();
    NsigE = Ns.getError();
    
    Nbkg = Nb.getVal();
    NbkgE = Nb.getError();

    sigma = width.getVal();
    sigmaE = width.getError();

    mu = mean.getVal();
    muE = mean.getError();

    a = a0.getVal();
    aE = a0.getError();
    
    llh = fitres->minNll();
    status = fitres->status();
    covQual = fitres->covQual();
    badNll = fitres->numInvalidNLL();
    edm = fitres->edm();

    tree->Fill();

    nruns++;
    /*
    TCanvas* canv_nominal = CreateCanvasNom("canv_nominal", fitres, dataToy, M, supM, infM, MassModel, sumgau, bkg1, Ns, Nb, width, width2, fs, mean,ptl,pth); 
    canv_nominal->Print(Form("plots/mass_BsFit_ptbins_ToyMC_%1.0f_%1.0f.png",ptl,pth));
    canv_nominal->Print(Form("plots/mass_BsFit_ptbins_ToyMC_%1.0f_%1.0f.pdf",ptl,pth));
    */
    if (nruns%10==0) cout<<nruns<<endl;
    delete dataToy;
    delete fitres;

  }

 tree->Write();
 //file->Write("",TObject::kOverwrite);
 
 cout<<" for end: "<<nruns<<endl;

}