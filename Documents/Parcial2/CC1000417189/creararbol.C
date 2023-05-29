void creararbol(){
ifstream in;
in.open("bc2dFit.csv");
if ( !in ) 
{
    cout << "No se pudo abrir el archivo" << endl;
    exit( 1 );
}

Float_t mass, time, error;
Int_t entradas = 0;

// auto f = TFile::Open("masayT.root","RECREATE");
TFile f ("masayT.root", "RECREATE");
TTree tree ("tree", "tree");

tree.Branch ("mass" , &mass, "mass/F");
tree.Branch ("time" , &time, "time/F");
tree.Branch ("error" , &error, "error/F");

string line;
getline(in,line);


while (getline(in, line)) {
    istringstream iss(line);
    iss >> mass >> time >> error;   
    if (!in.good()) break;  
    tree.Fill();
    entradas++;
}
cout << Form(" found %d points",entradas) << endl;

in.close();
f.Write();
}