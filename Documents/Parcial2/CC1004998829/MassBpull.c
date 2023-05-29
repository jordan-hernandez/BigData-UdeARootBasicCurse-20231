#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TFile.h"

using namespace RooFit; 

void save_resultTfit(ofstream& salida, RooRealVar nsig, RooRealVar nbkg, RooRealVar a0, RooRealVar mean,  RooRealVar sigma1)
{
  salida <<  nsig.getVal() << " " << nsig.getError() << " " << nbkg.getVal() << " " << nbkg.getError() <<" " <<  a0.getVal() << " " << a0.getError() << " " << mean.getVal() << " " << mean.getError() << " "
	<<  sigma1.getVal() << " " << sigma1.getError();
  cout << " el archivo 'resulTfit' se escribio bien" << endl; 
  return;
}

void MassBpull(){
  // I m p o r t   T T r e e   i n t o   a   R o o D a t a S e t
  // -----------------------------------------------------------
  TFile f("partc2.root");
  TTree* t = (TTree*) f.Get ("mytree");
  Float_t bcm, tau, etau;

  t->SetBranchAddress("bcm",&bcm);
  t->SetBranchAddress("tau",&tau);
  t->SetBranchAddress("etau",&etau);

  std::cout << "These are the columns bcm, tau, and etau:" << std::endl;
  for (auto branch : *t->GetListOfBranches()) {
    std::cout << "Branch: " << branch->GetName() << std::endl;
  }

  RooRealVar Bcm("bcm","bcm",6.05,6.5);
  
  RooDataSet *datatre = new RooDataSet("datatre", "datatre", RooArgSet(Bcm));

  for (int i=0; i<t->GetEntries();i++){
    t->GetEvent(i);
    if((tau/etau)<5.0) continue;
    Bcm=bcm;
    datatre->add(RooArgSet(Bcm));
  }

  datatre->Print("***\n V \n ***");
  ///
  RooRealVar mean("mean","mean of gaussians",6.05,6.5) ;
  RooRealVar sigma1("sigma1","sigma1 of gaussians",0.01,0.15) ;
  RooGaussian sig1("sig1","Signal component 1",Bcm,mean,sigma1) ;
  
  // Build Chebychev polynomial pdf
  RooRealVar a0("a0","a0",-10,10) ;
  RooChebychev bkg("bkg","Background",Bcm,a0) ;
    
  // Associated nsig/nbkg as expected number of events with sig/bkg _in_the_range_ "signalRange"
  RooRealVar nsig("nsig","number of signal events",0.,1500) ;
  RooRealVar nbkg("nbkg","number of background events",0,1500) ;
  //model
  RooAddPdf  model("model","sig+bkg", RooArgList(sig1,bkg), RooArgList(nsig,nbkg)) ;
  
  
  // S a m p l e   d a t a ,   f i t   m o d e l
  // -------------------------------------------

  RooFitResult* result = model.fitTo(*datatre,Extended(kTRUE),Save(kTRUE)); 

  // Plot datatre and PDF overlaid
  TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1000);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.42, 1, 1);
  //pad1->SetBottomMargin(0); // Sin margen inferior en el primer pad
  pad1->Draw();

  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.43);
  pad2->SetTopMargin(0); // Sin margen superior en el segundo pad
  pad2->SetBottomMargin(0.25); // Margen inferior más grande en el segundo pad
  pad2->Draw();

  pad1->cd();

  RooPlot* xframe = Bcm.frame();

  datatre->plotOn(xframe,Name("Datos"));
  model.plotOn(xframe,Name("Fit"));
  RooHist* hpullm2 = xframe->pullHist();
  model.plotOn(xframe,Components(sig1),LineColor(kRed),LineWidth(2),Name("Signal")); 
  model.plotOn(xframe, Components(bkg), LineColor(kGreen), LineStyle(kDashed),Name("Background")); 
  xframe->GetXaxis()->CenterTitle();
  xframe->GetYaxis()->SetNdivisions(505,1);
  xframe->GetXaxis()->SetNdivisions(505,1);   
  xframe->Draw();

  TLegend *legend = new TLegend(0.17, 0.67, 0.37, 0.87);
  legend->AddEntry(xframe->findObject("Signal"), "Signal");
  legend->AddEntry(xframe->findObject("Datos"), "Data");
  legend->AddEntry(xframe->findObject("Fit"), "Fit");
  legend->AddEntry(xframe->findObject("Background"), "Background");
  legend->SetBorderSize(0);
  legend->Draw();

  Double_t Mpsi = mean.getVal()*1000.0;
  Double_t MpsiE = mean.getError()*1000.0;

  Double_t G = sigma1.getVal()*1000.0;
  Double_t GE = sigma1.getError()*1000.0;

  TLegend *legpar = new TLegend(0.55,0.68,0.75,0.88);
  legpar->SetTextSize(0.035);
  legpar->SetFillColor(0);
  legpar->SetBorderSize(0);
  legpar->SetFillStyle(0);
  legpar->AddEntry("",Form("Mass(B_{c}^{+}) = %1.2f #pm %1.2f MeV",Mpsi,MpsiE),"");
  legpar->AddEntry("",Form("#sigma = %1.2f #pm %1.2f MeV",G,GE),"");
  legpar->AddEntry("",Form("N_{B_{c}^{+}} = %1.2f #pm %1.2f",nsig.getVal(),nsig.getError()),"");
  legpar->AddEntry("",Form("N_{bkg} = %1.2f #pm %1.2f",nbkg.getVal(),nbkg.getError()),"");
  legpar->Draw();
  xframe->SetXTitle("M(J/\\psi\\pi^{+}) [GeV]");
  xframe->SetYTitle("Events/10 MeV");

  datatre->Print("v");
  result->Print("v");

  pad2->cd();
  RooPlot* framem2 = Bcm.frame(Title("Pull Distribution"));
  framem2->addPlotable(hpullm2, "P");
  framem2->SetYTitle("(Data-Fit)/#sigma");
  framem2->SetXTitle("M(J/\\psi\\pi^{+}) [GeV]");
  framem2->GetXaxis()->CenterTitle();
  framem2->GetXaxis()->SetNdivisions(505,1);   
  framem2->Draw();
  TLine *line = new TLine(6.05, 0, 6.5, 0);
  line->SetLineColor(kRed);
  line->SetLineStyle(kDashed);
  line->Draw();


  canvas->SaveAs("mass&pull.pdf");

  //parameters output file
  ofstream salida_TotalFit("output_resulTFit_1_2.txt");
  salida_TotalFit.is_open();
  save_resultTfit(salida_TotalFit, nsig, nbkg, a0, mean, sigma1);
}