void LeerCSV() {
	// Abrir el archivo CSV:
	ifstream file;
	file.open("bc2dFit.csv");
	if (!file) {
		cout << "Error al leer el archivo" << endl;
		exit(1);
		}
	
	// Declarar variables para almacenar los datos del CSV:
	Float_t Masa, Tau, TauE;
	Int_t lineas = 0;
	
	// Crear un archivo ROOT para almacenar los datos organizados:
	auto f = TFile::Open("DatosCSV.root", "RECREATE");
	
	// Crear un árbol en el archivo ROOT:
	TTree mytree("MiArbol", "Datos del archivo csv");
	mytree.Branch("Masa", &Masa, "Masa/F"); 
	mytree.Branch("Tau", &Tau, "Tau/F");    
	mytree.Branch("TauE", &TauE, "TauE/F"); 
	
	// Leer la primera línea del archivo CSV (encabezados) y descartarla:
	string line;
	getline(file, line);
	
	// Leer el resto de líneas del archivo CSV
	while (getline(file, line)) {
		istringstream iss(line);
		iss >> Masa >> Tau >> TauE;
		
		if (!file.good()) break;
		mytree.Fill(); // Llenar el árbol con los datos de la línea actual
		lineas++;
		}

	// Cerrar el archivo CSV y guardar el archivo ROOT:
	file.close();
	f -> Write();
	}
