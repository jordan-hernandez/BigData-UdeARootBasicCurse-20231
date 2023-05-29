#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooGaussModel.h"
#include "RooGenericPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooAbsCollection.h"

using namespace RooFit;

void Fit2D() {
    Float_t Masa, Tau, TauE;

    // Abrir archivo ROOT y obtener el árbol de datos:
    TFile f("DatosCSV.root", "READ");
    TTree* tree = (TTree*)f.Get("MiArbol");

    // Asignar las ramas del árbol a las variables correspondientes:
    tree -> SetBranchAddress("Masa", &Masa);
    tree -> SetBranchAddress("Tau", &Tau);
    tree -> SetBranchAddress("TauE", &TauE);
    
    // Variables:
    RooRealVar M("M", "Masa", 6.05, 6.5);
    RooRealVar T("T", "Tau", 0.3, 2.6);

    // Conjunto de datos para almacenar los valores:
    RooDataSet data("data", "data", RooArgSet(M, T));
    
    Int_t entradas = tree -> GetEntries();
    cout << "Entradas: " << entradas << endl;
    
    // Iterar sobre las entradas del árbol:
    for (int evt = 0; evt < tree -> GetEntries(); evt++) {
        tree -> GetEvent(evt);
        
        // Aplicar condiciones de selección:
        if (isnan(TauE)) continue;                  // Omitir eventos con entrada vacía
        if (Masa < 6.05 || Masa > 6.5) continue;     // Condicion sobre Masa
        if (Tau < 0.3 || Tau > 2.6) continue;        // Condicion sobre Tau
        if (TauE < 0.0001 || TauE > 0.5) continue;   // Condicion sobre TauE
        if (Tau / TauE < 5.0) continue;              // Condicion sobre Tau y TauE

		// Agregar los valores al conjunto de datos:
        M = Masa;
        T = Tau;
        data.add(RooArgSet(M, T));
	    }
	    
	// Imprimir el conjunto de datos:
    data.Print();

    // Modelo para el ajuste (suma de señal y fondo):
    RooRealVar mean("mean", "mean", 6.27, 6.25, 6.3);
    RooRealVar sigma("sigma", "sigma", 0.015, 0.01, 0.04);
    RooGaussian gauss("gauss", "gauss", M, mean, sigma);

    RooRealVar c("c", "c", -10.0, 10.0);
    RooChebychev bg("bg", "Background", M, c);

    RooRealVar Ns("Ns", "Ns", 0., 1000);
    RooRealVar Nb("Nb", "Nb", 0., 1000);
    RooAddPdf ModeloMasa("ModeloMasa", "ModeloMasa", RooArgList(gauss, bg), RooArgList(Ns, Nb));

    // Modelo de tiempo de vida (exponencial):
    RooRealVar c2("c2", "c2", -10, 10);
    RooExponential Lifetime("Lifetime", "Lifetime", T, c2);

    // Modelo de masa y modelo de tiempo de vida:
    RooProdPdf MLifetime("MLifetime", "MLifetime", ModeloMasa, Conditional(Lifetime, T));

    // Fit (modelo de masa * modelo de tiempo de vida):
    RooFitResult* FitML = MLifetime.fitTo(data, Extended(), Minos(kFALSE), Save(kTRUE), NumCPU(4));
    FitML -> Print("v");

    RooPlot* frame = T.frame(Title("Tau"));

    // Crear un lienzo y dos pads para las gráficas:
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    TPad *pad1 = new TPad("p1", "", 0.01, 0.05, 0.49, 0.99);
	TPad *pad2 = new TPad("p2", "", 0.51, 0.05, 0.99, 0.99);
	
    // Establecer márgenes y dibujar los pads:
    pad1 -> SetLeftMargin(0.09);
    pad1 -> SetRightMargin(0.019);
    pad1 -> SetTopMargin(0.09);
    pad1 -> Draw();

    pad2 -> SetLeftMargin(0.09);
    pad2 -> SetRightMargin(0.019);
    pad2 -> SetTopMargin(0.09);
    pad2 -> Draw();
    
    pad1 -> cd();
    
    RooPlot* frame1 = M.frame(Title("Tau"));
    
    frame1 -> SetTitle("Masa B_{c} ");
    frame1 -> SetYTitle("Events");
    frame1 -> SetXTitle("M(B_{c})[GeV]");
    frame1 -> GetXaxis() -> CenterTitle();
    frame1 -> GetYaxis() -> CenterTitle();
    frame1 -> GetXaxis() -> SetNdivisions(505, 1);
    frame1 -> GetYaxis() -> SetNdivisions(505, 1);
    frame1 -> GetXaxis() -> SetTickLength(0.07);
    
    data.plotOn(frame1, Name("Data"), MarkerSize(0.8));
    MLifetime.plotOn(frame1, LineColor(kRed+2), LineWidth(2), Name("Fit"), XErrorSize(0));
    frame1 -> Draw();

    TLegend *l1 = new TLegend(0.18, 0.68, 0.58, 0.88);
    l1 -> SetTextSize(0.04);
    l1 -> SetFillColor(0);
    l1 -> SetBorderSize(0);
    l1 -> SetFillStyle(0);
    l1 -> AddEntry(frame1 -> findObject("Data"), " Data", "ep");
    l1 -> AddEntry(frame1 -> findObject("Fit"), " Fit", "l");
    l1 -> Draw();

    pad2 -> cd();

    RooPlot* frame2 = T.frame();
    frame2 -> SetTitle("Lifetime B_{c}");
    frame2 -> SetXTitle("Tau [ps]");
    frame2 -> SetYTitle("Events");
    frame2 -> GetYaxis() -> CenterTitle();
    frame2 -> GetXaxis() -> CenterTitle();
    frame2 -> GetYaxis() -> SetNdivisions(505, 1);
    frame2 -> GetXaxis() -> SetNdivisions(505, 1);
    frame2 -> GetXaxis() -> SetTickLength(0.07);
    
    data.plotOn(frame2, MarkerSize(0.8), XErrorSize(0), Name("Data"));
    MLifetime.plotOn(frame2, LineWidth(2), Name("Fit"), LineColor(kRed+2), MarkerSize(0.8), XErrorSize(0));
    frame2 -> Draw();

    TLegend *leg = new TLegend(0.18, 0.68, 0.58, 0.88);
    leg -> SetTextSize(0.04);
    leg -> SetFillColor(0);
    leg -> SetBorderSize(0);
    leg -> SetFillStyle(0);
    leg -> AddEntry(frame2 -> findObject("Data"), " Data", "ep");
    leg -> AddEntry(frame2 -> findObject("Fit"), " Fit", "l");
    leg -> Draw();

    canvas -> Draw();
    canvas -> Print("fit_tau.png");

    // Histograma 2D:
    TH1* Histo = MLifetime.createHistogram(" ", M, Binning(50), YVar(T, Binning(50)));
    Histo -> SetLineColor(kViolet-1);

    RooPlot* frame3 = M.frame(Title("2D Fit"));
    data.plotOn(frame3);
    MLifetime.plotOn(frame3, ProjWData(data, true));

    TCanvas* canvas1 = new TCanvas("canvas", "canvas", 1200, 400);
    gPad -> SetLeftMargin(0.20);
    Histo -> GetXaxis() -> SetTitleOffset(2.5);
    Histo -> GetYaxis() -> SetTitleOffset(4.5);
    Histo -> GetZaxis() -> SetTitleOffset(1.5);
    Histo -> SetTitle("Lifetime B_{c}");
    Histo -> GetYaxis() -> CenterTitle();
    Histo -> GetXaxis() -> CenterTitle();
    Histo -> GetZaxis() -> CenterTitle();
    Histo -> GetYaxis() -> SetNdivisions(505, 1);
    Histo -> GetXaxis() -> SetNdivisions(505, 1);
    Histo -> GetZaxis() -> SetNdivisions(505, 1);
    Histo -> GetXaxis() -> SetTickLength(0.07);
    Histo -> SetTitleOffset(0.5, "Y");
    Histo -> Draw("surf");

    canvas1 -> Draw();
    canvas1 -> Print("fit_2d.png");
	}
