void createTTree()
{
    TFile dataBc( "dataBc.root", "RECREATE" );
    TTree treeBc("treeBc", "treeBc");

    // ReadFile fills the tree    from txt file
    treeBc.ReadFile("bc2dFit.txt", "Bcmass/F:tau/F:tauerror/F");

    // Scan prints out data
    treeBc.Scan("Bcmass:tau:tauerror");

    // Write saves tree in file
    dataBc.Write();
}