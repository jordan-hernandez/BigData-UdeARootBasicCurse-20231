void read2csv() {
  float Bcmass, tau, tauerror;

  ifstream in("bc2dFit.csv");
  if (!in) {
    cout << "No se pudo abrir el archivo" << endl;
    exit(1);
  }

  TFile f("bc2_data.root", "RECREATE");
  TTree tree("tree", "tree");
  tree.Branch("Bcmass", &Bcmass, "Bcmass/F");
  tree.Branch("tau", &tau, "tau/F");
  tree.Branch("tauerror", &tauerror, "tauerror/F");

  string uslessLine;
  getline(in, uslessLine);

  while (in >> Bcmass >> tau >> tauerror) {
    tree.Fill();
  }

  in.close();
  f.Write();
  f.Close();
}
