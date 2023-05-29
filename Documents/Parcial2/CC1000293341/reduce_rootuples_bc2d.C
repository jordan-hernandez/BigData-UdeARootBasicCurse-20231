using namespace std;

void reduce_rootuples_bc2d() {
// read file 
   ifstream in;
   in.open("bc2dFit.csv");
    if ( !in ) 
    {
      cout << "No se pudo abrir el archivo" << endl;
      exit( 1 );
    }
 
   Float_t Bcmass,tau, tauerror;
   Int_t nlines = 0;

   //auto f = TFile::Open("reducetree_Bc2d_AOD_best.root","RECREATE");
   TString outfileName(Form("reducetree_Bc2d_AOD_best.root")); 
    
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TTree mytree( "mytree" , "data from csv file" );
   mytree.Branch ("Bcmass" , &Bcmass, "Bcmass/F");
   mytree.Branch ("tau" , &tau, "tau/F");
   mytree.Branch ("tauerror" , &tauerror, "tauerror/F");

    // Skip the header line.
   std::string header;
   std::getline(in, header);

   while (std::getline(in, header)) {
      istringstream iss(header);
      in >> Bcmass >> tau >> tauerror;
      if (!in.good()) break;
      mytree.Fill();
      nlines++;
   }

   outputFile->Write("",TObject::kOverwrite);

   cout << Form(" found %d points",nlines) << endl;

   in.close();
 
   //f->Write();
}