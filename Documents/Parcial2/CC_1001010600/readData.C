void readData(){
  Float_t mass, tau, tErr;

  ifstream in;
  in.open("../bc2dFit.csv");
  if (!in) {
    cout << "No se pudo abrir el archivo" << endl;
    exit(1);
  }

  //Crea el root file y el tree
	TFile f("BcData.root", "RECREATE");
  TTree tree("tree", "tree");

  tree.Branch("mass", &mass, "mass/F");
  tree.Branch("tau", &tau, "tau/F");
  tree.Branch("tErr", &tErr, "tErr/F");

  //Lee la primera linea que contiene los nombres de las variables
  string uslessLine;
  getline(in, uslessLine);

  while(1){
    in >> mass >> tau >> tErr;
    if(!in.good()) break;
    tree.Fill();
  }

  in.close();
  f.Write();
}
