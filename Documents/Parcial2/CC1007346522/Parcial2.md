# Parcial 2
Big Data - UdeA
#### Santiago Galvis Duque
santiago.galvis@udea.edu.co

### Orden de los archivos  
- **Leyendo_Data.C:** Lectura del archivo .txt y conversión a .root
- **Data._Fit.C:** Fitting de la Masa del mesón $B_c$ con respectivo Pull; el plot se guarda como *MasaBc_FiteoyPull.png* en la carpeta ***plots***
- **ToyMonteCarlo.C:** Generación de datos aleatorios para probar el fit; crea archivo *ToyMC.root*
- **Analisis_ToyMonteCarlo.C:** Analisis de los datos obtenidos con el método de MonteCarlo, crea tres plots que muestran como varía los parámetros del fit, la señal y el background en los archivos *Pull_Mu_ToyMC.png*, *Pull_Yi_ToyMC.png* y *Pull_Yib_ToyMC.png* en la carpeta ***plots*** respectivamente.
- **Data._Fit_2D.C:** Fitting 2D de la Masa y Tiempo de vida $\tau$ de manera que la distribución 2D es condicional; el plot se guarda como *MasayTau_Fiteo.png* en la carpeta ***plots***