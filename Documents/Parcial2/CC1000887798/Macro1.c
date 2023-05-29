int Macro1() {

ifstream in;
in.open("bc2dFit.csv");
if ( !in ) 
{
  cout << "No se pudo abrir el archivo" << endl;
  exit( 1 );
}

TFile file("parcial2.root", "RECREATE");


TTree tree("tree", "tree");


Float_t masa, tau, errtau;
Int_t entradas = 0;

tree.Branch("masa", &masa, "masa/F");
tree.Branch("tau", &tau, "tau/F");
tree.Branch("errtau", &errtau, "errtau/F");

string line;
getline(in, line);

while (getline(in, line)) {
    istringstream iss(line);
    iss >> masa >> tau >> errtau;
    if (!in.good()) break;  
    tree.Fill();
    entradas++;
}
cout << Form(" found %d points",entradas) << endl;

tree.Write();
file.Close();

return 0;
}
