using namespace RooFit;

void save_resultTfit(ofstream& salida, RooRealVar Ns, RooRealVar Nb, RooRealVar c, RooRealVar mean,  RooRealVar width)
{
  salida <<  Ns.getVal() << " " << Ns.getError() << " " << Nb.getVal() << " " << Nb.getError() <<" " <<  c.getVal() << " " << c.getError() << " "  
	 << " " << mean.getVal() << " " << mean.getError() << " "
	 <<  width.getVal() << " " << width.getError() 
    ;
  cout << " el archivo se escribio bien" << endl; 
  return;
}

void save_result(ofstream& salida, RooRealVar Ns, RooRealVar Nb)
{
  salida <<  Ns.getVal() << " " << Ns.getError() << " " << Nb.getVal() << " " << Nb.getError();
  cout << " el archivo se escribio bien" << endl; 
  return;
}
int Macro2()
{
Float_t masa, tau , errtau;


TFile f ("parcial2.root", "READ");
TTree *tree= (TTree*)f.Get("tree");
tree->SetBranchAddress ("masa", &masa);
tree->SetBranchAddress ("tau", &tau);
tree->SetBranchAddress ("errtau", &errtau);

Double_t Mmin = 6.05; 
Double_t Mmax = 6.5; 
 

Long64_t nentries = tree->GetEntries();
cout<<" Entries : "<<nentries<<endl;
  
//------------------------------
RooRealVar M("M","",Mmin,Mmax);
RooDataSet data("data","data",RooArgSet(M));


for(int evt = 0; evt < nentries; evt++)
    {
      tree->GetEvent(evt);
      if(isnan(errtau)) continue;
      if((tau/errtau)<5.0) continue;
      M=masa;    
      data.add(RooArgSet(M));
    }
tree->ResetBranchAddresses();

 //FIT
RooRealVar mean("mean"," Mass mean",6.2751,6.25,6.3,"GeV");
RooRealVar width("width"," Mass width",0.015,0.01,0.04,"GeV");
RooGaussian gaus("Sig"," Signal PDF",M,mean,width);

RooRealVar c("c","c",-10.0,10.0);
RooChebychev  bkg1("bkg1","Exp. Background",M,c);

RooRealVar Ns("Ns","Ns",0.,2000);
RooRealVar Nb("Nb","Nb",0.,2000);   

RooAddPdf MassModel("MassModel","MassModel",RooArgList(gaus,bkg1),RooArgList(Ns,Nb));
RooFitResult* fitres = MassModel.fitTo(data,Extended(),Minos(kFALSE),Save(kTRUE), NumCPU(4));

data.Print("v"); 
fitres->Print("v");
Double_t nbin = ((Mmax-Mmin)/0.007)+1;  // 0.005 is the bin size

TCanvas *c1 = new TCanvas("basic", "basic",50,50, 800, 600);
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
pad2->SetFillColor(0);
pad2->SetGridx(0);
pad2->SetGridy(0);

pad1->Draw();
pad2->Draw();
pad1->cd(); 

RooPlot* framemm = M.frame(Title("")) ;
RooPlot* frame = M.frame(Mmin,Mmax,nbin) ;
data.plotOn(frame,DataError(RooAbsData::SumW2),MarkerSize(.75),XErrorSize(0));
MassModel.plotOn(frame,LineColor(kRed),LineWidth(2),Name("Fit")); 

RooHist* hpullm2 = frame->pullHist() ;


MassModel.plotOn(framemm,Components(gaus),LineColor(kBlue),LineWidth(2),Name("Signal")); 
MassModel.plotOn(frame,Components(bkg1),LineStyle(kDashed)  ,LineColor(kGreen),LineWidth(2),Name("bkg")); 
MassModel.plotOn(frame,Components(gaus),LineStyle(kDashed) ,LineColor(kBlue),LineWidth(2),Name("Signal")); 
framemm->Draw();

frame->SetTitle("Fit to the B_{c}^{+} mass");
frame->SetTitleSize(0.8,"");
frame->SetYTitle("Events");        
frame->SetLabelSize(0.06,"XY");
frame->SetTitleSize(0.19,"XY");
frame->GetYaxis()->CenterTitle();   
frame->GetXaxis()->CenterTitle();
frame->GetYaxis()->SetNdivisions(505,1);
frame->GetXaxis()->SetNdivisions(505,1);   
frame->GetXaxis()->SetDecimals(1);
frame->GetXaxis()->SetTickLength(0.06); 
frame->SetTitleSize(0.11,"X");
frame->SetTitleSize(0.1,"Y");  
frame->SetTitleOffset(0.9,"X");
frame->SetTitleOffset(0.42,"Y");
frame->SetMinimum(-4.0);   
frame->Draw();

//TLegend *leg = new TLegend(0.15,0.68,0.35,0.88);
TLegend *leg = new TLegend(0.18,0.58,0.38,0.88); 
leg->SetTextSize(0.06);
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->AddEntry(frame->findObject("Data")," Data","ep"); 
leg->AddEntry(frame->findObject("Fit")," Fit result","l");
leg->AddEntry(frame->findObject("Signal")," Signal","l");
leg->AddEntry(frame->findObject("bkg"),"Comb. backg.","l");
leg->Draw();
Double_t Mpsi = mean.getVal()*1000.0;
Double_t MpsiE = mean.getError()*1000.0;
Double_t G = sqrt( width.getVal()*width.getVal() + (width.getVal()*width.getVal() ))*1000.0;
Double_t GE = (1/G)*sqrt( (width.getVal()*width.getVal())*(width.getError()*width.getError()) + (width.getVal()*width.getVal())*(width.getError()*width.getError()) )*1000.0*1000.0;

TLegend *legpar = new TLegend(0.6,0.58,0.8,0.88);
legpar->SetTextSize(0.06);
legpar->SetTextFont(42);
legpar->SetFillColor(0);
legpar->SetBorderSize(0);
legpar->SetFillStyle(0);
legpar->AddEntry("",Form("M(B_{c}^{+}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
legpar->AddEntry("",Form("N_{B_{c}^{+}} = %1.0f #pm %1.0f",Ns.getVal(),Ns.getError()),"");
legpar->AddEntry("",Form("N_{bkg} = %1.0f #pm %1.0f",Nb.getVal(),Nb.getError()),"");
legpar->Draw();


pad2->cd();


RooPlot* framem2 = M.frame() ;
framem2->addPlotable(hpullm2,"P") ;
framem2->SetTitle("");
framem2->SetYTitle(" (Data-Fit)/#sigma");
framem2->SetXTitle("Mass B_{c}^{+}  Gev");
framem2->SetLabelSize(0.1,"XY");
framem2->SetTitleSize(0.11,"X");
framem2->SetTitleSize(0.135,"Y");  
framem2->GetYaxis()->CenterTitle();   
framem2->GetXaxis()->CenterTitle();
framem2->GetYaxis()->SetNdivisions(505,1);
framem2->GetXaxis()->SetNdivisions(505,1);
framem2->GetXaxis()->SetTickLength(0.07);   
framem2->SetTitleOffset(0.9,"X");
framem2->SetTitleOffset(0.31,"Y");
framem2->SetMaximum(4.9);
framem2->SetMinimum(-4.9);

framem2->Draw();

TLine *line1 = new TLine(Mmin,0.0,Mmax,0.0);
line1->SetLineColor(1);
//line1->SetLineStyle(2);
line1->SetLineWidth(1);
line1->Draw();


c1->Modified();
c1->Print("fit.png");

//parameters output file
ofstream salida_TotalFit("datos_fit.txt");
salida_TotalFit.is_open();
save_resultTfit(salida_TotalFit, Ns, Nb, c, mean, width);




return 0;
}