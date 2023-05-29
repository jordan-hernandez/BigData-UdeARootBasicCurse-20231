#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TMathBase.h"

using namespace RooFit;
using namespace TMath;

int FitMasaBc() {
    Float_t Masa, Tau, TauE;

    // Abrir archivo ROOT y obtener el árbol de datos:
    TFile f("DatosCSV.root", "READ");
    TTree* tree = (TTree*)f.Get("MiArbol");

    // Asignar las ramas del árbol a las variables correspondientes:
    tree -> SetBranchAddress("Masa", &Masa);
    tree -> SetBranchAddress("Tau", &Tau);
    tree -> SetBranchAddress("TauE", &TauE);

    // Variable para la masa:
    RooRealVar M("M", "Masa", 6.05, 6.5);
    // Conjunto de datos para almacenar los valores de masa:
    RooDataSet data("data", "data", RooArgSet(M));

    Int_t entradas = tree -> GetEntries();
    cout << "Entradas: " << entradas << endl;

    // Iterar sobre las entradas del árbol:
    for (int evt = 0; evt < tree -> GetEntries(); evt++) {
        tree -> GetEvent(evt);
        
        // Aplicar condiciones de selección:
        if (Tau / TauE < 5.0) continue; // Condición para el error
        if (isnan(TauE)) continue; // Omitir eventos con entrada vacía

        // Agregar el valor de masa al conjunto de datos:
        M = Masa;
        data.add(RooArgSet(M));
	    }

    // Imprimir el conjunto de datos:
    data.Print();

    // Restaurar las ramas del árbol:
    tree -> ResetBranchAddresses();

    // Modelo de la señal (Gaussiana):
    RooRealVar mean("mean", "mean", 6.27, 6.25, 6.3);
    RooRealVar sigma("sigma", "sigma", 0.015, 0.01, 0.04);
    RooGaussian gauss("gauss", "Gaussian", M, mean, sigma);

    // Modelo del background (Chebychev):
    RooRealVar c("c", "c", -10.0, 10.0);
    RooChebychev bg("bg", "Background", M, c);

    // Número de eventos de señal y fondo:
    RooRealVar Ns("Ns", "Ns", 0., 2000);
    RooRealVar Nb("Nb", "Nb", 0., 2000);

    // Modelo para el ajuste (suma de señal y fondo):
    RooAddPdf ModeloMasa("ModeloMasa", "ModeloMasa", RooArgList(gauss, bg), RooArgList(Ns, Nb));

    // Realizar el ajuste al conjunto de datos:
    RooFitResult* Fit = ModeloMasa.fitTo(data, Extended(), Minos(kFALSE), Save(kTRUE), NumCPU(4));

    // Imprimir resultados del ajuste
    data.Print("v");
    Fit -> Print("v");

    // Guardar los parámetros del ajuste en un archivo de texto:
    ofstream SalidaFit(Form("fit_param.txt"));
    SalidaFit.is_open();
    SalidaFit << Ns.getVal() << " " << Ns.getError() << " " << Nb.getVal() << " " << Nb.getError() << " " << c.getVal() << " " 
			  << c.getError() << " " << mean.getVal() << " " << mean.getError() << " " << sigma.getVal() << " " << sigma.getError();

    // Crear un marco para la gráfica:
    RooPlot* frame = M.frame(Title("Fit para la masa del meson Bc"));

    // Graficar los datos, el ajuste, la señal y el fondo en el marco:
    data.plotOn(frame, MarkerSize(0.8), XErrorSize(0), Name("Data"));
    ModeloMasa.plotOn(frame, LineColor(kBlue-2), LineWidth(2), Name("Fit"));
    ModeloMasa.plotOn(frame, Components(gauss), LineColor(kMagenta+2), LineWidth(2), Name("Signal"));
    ModeloMasa.plotOn(frame, Components(bg), LineColor(kCyan+1), LineWidth(2), LineStyle(kDashed), Name("Background"));

    // Crear un lienzo y dos pads para las gráficas:
    TCanvas* canvas = new TCanvas("canvas", "canvas", 50, 50, 1200, 800);
    TPad* pad1 = new TPad("pad1", "pad1", 0.01, 0.411, 0.9903769, 0.99);
    TPad* pad2 = new TPad("pad2", "pad2", 0.01, 0.01, 0.9903769, 0.41);

    // Establecer márgenes y dibujar los pads:
    pad1 -> SetLeftMargin(0.09);
    pad1 -> SetRightMargin(0.019);
    pad1 -> SetTopMargin(0.09);
    pad1 -> SetBottomMargin(0.0);
    pad1 -> Draw();
        
    pad2 -> SetLeftMargin(0.09);
    pad2 -> SetRightMargin(0.019);
    pad2 -> SetTopMargin(0.0);
    pad2 -> SetBottomMargin(0.25);
    pad2 -> Draw();

    pad1 -> cd();

    // Configurar el marco de la gráfica:
    frame -> SetYTitle("Events");
    frame -> SetTitleSize(0.06, "Y");
    frame -> SetTitleOffset(0.5, "Y");
	frame -> GetXaxis() -> CenterTitle();
    frame -> GetYaxis() -> CenterTitle();
    frame -> GetXaxis() -> SetNdivisions(505, 1);
    frame -> GetYaxis() -> SetNdivisions(505, 1);
    frame -> GetXaxis() -> SetTickLength(0.07);
    frame -> SetLabelSize(0.04, "XY");
    frame -> SetMinimum(-5);
    frame -> Draw("A");

    // Crear una leyenda para los objetos de la figura:
    TLegend* lfig = new TLegend(0.18, 0.58, 0.38, 0.88);
    lfig -> SetTextSize(0.06);
    lfig -> SetFillColor(0);
    lfig -> SetBorderSize(0);
    lfig -> SetFillStyle(0);
    lfig -> AddEntry(frame -> findObject("Data"), " Data", "ep");
    lfig -> AddEntry(frame -> findObject("Fit"), " Fit", "l");
    lfig -> AddEntry(frame -> findObject("Signal"), " Signal", "l");
    lfig -> AddEntry(frame -> findObject("Background"), "Background", "l");
    lfig -> Draw();

    // Crear una leyenda para el resultado del ajuste:
    TLegend* lfit = new TLegend(0.65, 0.70, 0.9, 0.88);
    lfit -> SetTextSize(0.06);
    lfit -> SetFillColor(0);
    lfit -> SetBorderSize(0);
    lfit -> SetFillStyle(0);
    lfit -> AddEntry("", Form("M(B_{c}) = %1.2f #pm %1.2f MeV", mean.getVal() * 1000.0, mean.getError() * 1000.0), "");
    lfit -> AddEntry("", Form("N_{B_{c}} = %1.0f #pm %1.0f", Ns.getVal(), Ns.getError()), "");
    lfit -> AddEntry("", Form("N_{bkg} = %1.0f #pm %1.0f", Nb.getVal(), Nb.getError()), "");
    lfit -> Draw();

    // Crear un marco para el pull:
    RooPlot* Mframe = M.frame(6.05, 6.5, ((6.5 - 6.05) / 0.008) + 1);
    data.plotOn(Mframe, DataError(RooAbsData::SumW2), MarkerSize(1.0), XErrorSize(0));
    ModeloMasa.plotOn(Mframe, DrawOption("F"), FillColor(0), LineWidth(2), Name("FitPull"));
    RooHist* hpullm2 = Mframe -> pullHist();

    pad2 -> cd();

    // Configurar el marco para el pull:
    RooPlot* framem2 = M.frame(Title(" "));
    framem2 -> SetYTitle("(Data-Fit)/#sigma");
    framem2 -> addPlotable(hpullm2, "P");
    framem2 -> GetXaxis() -> CenterTitle();
    framem2 -> GetYaxis() -> CenterTitle();
    framem2 -> SetXTitle("M(B_{c})[GeV]");
    framem2 -> GetYaxis() -> SetNdivisions(505, 1);
    framem2 -> GetXaxis() -> SetNdivisions(505, 1);
    framem2 -> GetXaxis() -> SetTickLength(0.07);
    framem2 -> SetTitleOffset(0.35, "Y");
    framem2 -> SetTitleSize(0.08, "Y");
    framem2 -> SetLabelSize(0.06, "XY");
    framem2 -> SetTitleSize(0.08, "X");
    framem2 -> Draw();

    // Guardar el lienzo:
    canvas -> Draw();
    canvas -> Print("fit_meson.png");

    return 0;
}
