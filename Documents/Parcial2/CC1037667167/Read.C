//Este macro lee el archivo de datos bc2dFit.csv, guardamos sus datos en un arbol e imprime el arbol en la terminal!!.

#include<iostream>

{
//Crear TFile que contendrá el TTree
TFile* file = new TFile("Data_b_Hadron.root" , "RECREATE");
TTree* tree = new TTree("Data_tree" , "Data_tree");

//Abrimos archivo lectura y verificamos que este bien abierto.
fstream ArchivoEntrada("bc2dFit.txt" , ios::in);
if(!ArchivoEntrada){cout<<"No se pudo abrir el archivo: "<<endl;exit(1);}  

//leemos los datos del archivo:
tree->ReadFile("bc2dFit.txt" , "Mass/D:Tau/D:eTau/D");

//Imprimimos cosas
tree->Print();
tree->Scan();

//Guardar y cerrar!
tree->Write();
ArchivoEntrada.close();
file->Close();
}