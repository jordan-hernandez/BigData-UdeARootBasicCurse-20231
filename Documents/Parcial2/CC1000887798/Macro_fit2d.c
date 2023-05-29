using namespace RooFit;
void Macro_fit2d()
{
Float_t masa, tau, errtau;

TFile f("parcial2.root", "READ");
TTree *tree = (TTree*) f.Get("tree");

tree->SetBranchAddress("masa", &masa);
tree->SetBranchAddress("tau", &tau);
tree->SetBranchAddress("errtau", &errtau);

Double_t Mmin = 6.05; 
Double_t Mmax = 6.5; 

RooRealVar M("M","Mass B_{c}^{+}",Mmin,Mmax);
RooRealVar Tau("Tau","#tau (ps)",0.3,2.6);
  
RooRealVar mean("mean"," Mass mean",6.2751,6.25,6.3,"GeV");
RooRealVar width("width"," Mass width",0.015,0.01,0.04,"GeV");
RooGaussian Sig("Sig"," Signal PDF",M,mean,width);

RooRealVar c("c","c",-10.0,10.0);
RooChebychev  bkg1("bkg1","Exp. Background",M,c);

RooRealVar Ns("Ns","Ns",0.,2000);
RooRealVar Nb("Nb","Nb",0.,2000);  

RooDataSet data("data","data",RooArgSet(M,Tau ));

RooAddPdf MassModel("MassModel","MassModel",RooArgList(Sig,bkg1),RooArgList(Ns,Nb));

RooRealVar pa("pa","pa",-10,10);
RooExponential tauexp("tauexp","tauexp",Tau,pa);

// RooProdPdf model("model","model",RooArgList(MassModel,tauexp));
RooProdPdf model("model", "model", MassModel, Conditional(tauexp,Tau ));





for(int event = 0; event < tree->GetEntries(); event++)
{
   tree->GetEntry(event);
   if(masa<=6.0 || masa>=6.5) continue;
   if(tau<=0.3 || tau>=2.6) continue;
   if(errtau<=.0001 || errtau>=.5) continue;
   if((tau/errtau)<5.0) continue;

   M = masa;
   Tau = tau;
   errtau = errtau;
   data.add(RooArgSet(M,Tau));
}

data.Print();

RooFitResult* r = model.fitTo(data,Save(),Minos(kTRUE),Extended(kTRUE));


RooPlot* framem = M.frame();
framem->SetTitle("Mass for B_{c}^{+} meson");
framem->SetXTitle("Mass (GeV)");

data.plotOn(framem,MarkerSize(0.5),Name("Data"));
model.plotOn(framem, LineColor(kBlue),Name("Fit"));
framem->GetXaxis()->SetNdivisions(505);
framem->Draw();

RooPlot* framet = Tau.frame();
framet->SetTitle("Life time for B_{c}^{+} meson");
framet->SetXTitle("#tau (ps)");
data.plotOn(framet,MarkerSize(0.5));
model.plotOn(framet, LineColor(kBlue),Name("fit"));
framet->Draw();


TH1 *hh_decay = model.createHistogram("hh_decay", M, Binning(50), YVar(Tau, Binning(50)));
hh_decay->SetLineColor(kBlue);


TCanvas *c2 = new TCanvas("c2", "", 1200, 400);
c2->Divide(3);
c2->cd(1);
gPad->SetLeftMargin(0.20);
hh_decay->SetTitle("Decay time vs Mass");
hh_decay->GetZaxis()->SetTitleOffset(2.5);
hh_decay->GetZaxis()->SetTitle("Events");
hh_decay->GetXaxis()->SetTitleOffset(2);
hh_decay->GetYaxis()->SetTitleOffset(2.3);
hh_decay->GetXaxis()->SetLabelSize(0.03);
hh_decay->GetZaxis()->SetNdivisions(505);
hh_decay->GetYaxis()->SetNdivisions(505);
hh_decay->GetXaxis()->SetNdivisions(505);
hh_decay->GetYaxis()->CenterTitle();
hh_decay->GetXaxis()->CenterTitle();
hh_decay->Draw("surf");


c2->cd(2);
gPad->SetLeftMargin(0.15);
framet->GetYaxis()->SetTitleOffset(1.6);
framet->Draw();

TLegend *leg2 = new TLegend(0.42,0.550,0.8,0.88); 
leg2->SetTextSize(0.04);
leg2->SetFillColor(0);
leg2->SetBorderSize(0);
leg2->SetFillStyle(0);
leg2->AddEntry(framet->findObject("data"),"       data","ep"); 
leg2->AddEntry(framet->findObject("fit"),"        Fit result","l");
leg2->AddEntry("","        ","");
leg2->AddEntry("","        ","");
leg2->AddEntry("","        ","");
leg2->AddEntry("",Form("#tau = %1.3f #pm %1.3f ps",pa.getVal(),pa.getError()),"");
leg2->Draw();


c2->cd(3);
gPad->SetLeftMargin(0.15);
framem->GetYaxis()->SetTitleOffset(1.6);
framem->Draw();

TLegend *leg3 = new TLegend(0.15,0.70,0.8,0.88); 
leg3->SetTextSize(0.04);
leg3->SetFillColor(0);
leg3->SetBorderSize(0);
leg3->SetFillStyle(0);
leg3->AddEntry(framem->findObject("Data"),"data","ep"); 
leg3->AddEntry(framem->findObject("Fit"),"Fit result","l");
leg3->Draw();



c2->SaveAs("fit2d.png");

}