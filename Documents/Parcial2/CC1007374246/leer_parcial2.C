void leer_parcial2() {
// Este código permite leer el archivo bc2Fit.csv y almacena sus datos en Datos_organizados.root
   ifstream in;
   in.open("bc2dFit.csv");
    if ( !in ) 
    {
      cout << "No se pudo abrir el archivo" << endl;
      exit( 1 );
    }
 
   Float_t mass, tau, taue;
   Int_t nlines = 0;
   auto f = TFile::Open("Datos_organizados.root","RECREATE");
   TTree mytree( "mytree" , "data from csv file" );
   mytree.Branch ("mass" , &mass, "mass/F");
   mytree.Branch ("tau" , &tau, "tau/F");
   mytree.Branch ("taue" , &taue, "taue/F");

   string line;
   getline(in,line);

   while (getline(in, line)) {
      istringstream iss(line);
      iss >> mass >> tau >> taue;
      if (!in.good()) break;
      mytree.Fill();
      nlines++;
   }
   //printf(" found %d points\n",nlines);
   cout << Form(" found %d points",nlines) << endl;

   in.close();
   f->Write();
}
