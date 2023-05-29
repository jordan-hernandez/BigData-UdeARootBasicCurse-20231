void read() {
  float Bcmass, tau, tauerror; 
  
  ifstream in;
  in.open("../bc2dFit.csv");

  if (!in) {
    cout << "No se pudo abrir el archivo" << endl;
    exit(1);
  }

	TFile f("bc2.root", "RECREATE");
  TTree tree("tree", "tree");

  tree.Branch("Bcmass", &Bcmass, "Bcmass/F");
  tree.Branch("tau", &tau, "tau/F");
  tree.Branch("tauerror", &tauerror, "tauerror/F");

  string uslessLine;
  getline(in, uslessLine);

  while(1){
    in >> Bcmass >> tau >> tauerror;
    if(!in.good()) break;
    tree.Fill();
  }

  in.close();
  f.Write();
  f.Close();
}