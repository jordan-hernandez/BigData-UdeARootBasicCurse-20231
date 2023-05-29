void readData2File()
{

    std::ifstream inputFile("bc2dFit.csv");
    if (!inputFile)
    {
        std::cerr << "Error: No se pudo leer el archivo csv." << std::endl;
        return;
    }
    
    TFile* outputFile = TFile::Open("Bc_data.root", "RECREATE");
    TTree* tree = new TTree("tree", "Bc mass data");
    
    Double_t bcmass, tau, tauerror;
    tree->Branch("Bcmass", &bcmass);
    tree->Branch("Tau", &tau);
    tree->Branch("TauError", &tauerror);
    
    std::string line;
    Int_t nlines = 0;
    getline(inputFile, line);
    while (std::getline(inputFile, line))
    {
        
        std::istringstream iss(line);
        iss >> bcmass >> tau >> tauerror;
        if (nlines < 5) cout<< Form("x=%1f, y=%1.3f, z=%1.3f",bcmass,tau,tauerror) << endl;

        tree->Fill();
        nlines++;
    }
    
    cout << Form(" found %d points",nlines) << endl;
    inputFile.close();
    
    tree->Write();
    tree->MakeClass("myTreeClass");

    outputFile->Close();
    //delete outputFile;
}