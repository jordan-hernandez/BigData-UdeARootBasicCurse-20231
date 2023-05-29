{  // read file 
   ifstream in;
   in.open("bc2dFit.csv");
    if ( !in ) 
    {
      cout << "No se pudo abrir el archivo" << endl;
      exit( 1 );
    }
   Float_t bcm, tau,etau;
   auto f = TFile::Open("partc2.root","RECREATE");
   TTree mytree( "mytree" , "data from csv file" );
   mytree.Branch ("bcm" , &bcm, "bcm/F");
   mytree.Branch ("tau" , &tau, "tau/F");
   mytree.Branch ("etau" , &etau, "etau/F");  
   // Saltar la primera línea del archivo
   std::string line;
   std::getline(in, line);

   while (getline(in,line)) {
      istringstream iss(line);
      iss >> bcm >> tau >> etau;
      if (!in.good()) break;    
      mytree.Fill();
    }

   in.close();
   f->Write();
}
