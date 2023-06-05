#include <TROOT.h>
#include "TMath.h"
#include <iostream>
#include <fstream>
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooAbsPdf.h"
#include "RooArgusBG.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TLatex.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooNumIntConfig.h"
#include "TObject.h"

#ifndef SIMULATOR_H 
#define SIMULATOR_H

//using namespace RooFit;


  /*Esta clase recibe como input un modelo pdf con sus respectivos argumentos y genera a partir de ello
    1.) datos sinteticos
    2.) un fit a esos datos sinteticos
    3.) un pull y montecarlo de dichos datos
    
    todo esto solo con findes educativos de formacion*/

class Simulator : public TObject {
  
  private:
  int nbin;
  RooRealVar   *Obs; /*todo esto hay que definirlo como puntero porque los objetos de root borraron el operador asignacion =*/
  RooAbsPdf    *Model;
  RooDataSet   *DataSet;
  RooFitResult *FitResult;
    
  public:
  Simulator(); // Default constructor
  Simulator(RooAbsPdf &  , RooRealVar &  , int const&); 
  //TCanvas* MainPlot(RooAbsPdf *Model ,RooDataSet *DataSet);            
  TCanvas* MainPlot(); 
  TCanvas* McPlot(const RooArgSet & ); 

  ~Simulator();

  
     

ClassDef(Simulator, 1);
};

#endif





