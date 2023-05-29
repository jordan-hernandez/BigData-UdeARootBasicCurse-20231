#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"

using namespace std;

int Leyendo_Data()
{     
      // Verificamos la existencia y la correcta aperture del archivo
      ifstream in;

      // El archivo a trabajar será en formato .txt porque el original no coincide con formato .csv
      in.open("bc2dFit.txt");
      if ( !in ) 
      {     
            // Si no exite, verificar que exista el archivo en formato .csv
            in.open("../bc2dFit.csv");

            if ( !in){ 
            cout << "No se pudo abrir el archivo" << endl;
            exit( 1 );
            }
            else {
            // Cambiamos al formato al archivo 
                  filesystem::copy_file("../bc2dFit.csv", "bc2dFit.txt");
                  cout << "Archivo renombrado y abierto correctamente" << endl;}
      }
      else {
            cout << "Archivo abierto correctamente" << endl;
      }
      
      // Creamos el archivo .root
      TFile File1( "RootFile1.root", "RECREATE" );

      // Creamos el TTree
      auto *Tree1 = new TTree( "Tree1" , "Tree of data for my analysis" );

      // Llenamos el TTree con los datos del archivo .txt
      // Las columnas del archivo son: Masa Bc (GeV), Tau (ps), Error de Tau (ps)
      Tree1->ReadFile("bc2dFit.txt","mass/F:Tau/F:Tau_err/F");                                                             

      // Verificando que el TTree se haya llenado correctamente
      Tree1->Print();

      // Para ver a mayor profundidad:
      // Tree1->Scan();

      // Cerrando el archivo utilizado
      in.close();

      // Escribir cambios en el archivo .root
      File1.Write();
      
      return 0;
}   
