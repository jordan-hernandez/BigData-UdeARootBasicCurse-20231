#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>

using namespace std;

int main(){

  // grupos del 01 al 05
  
  srand( time( NULL ) );
  
  for ( int tiro = 1; tiro <= 20; tiro++ ){
    cout << 1 + rand() % 5 << endl; 
  }
  

  
return 0; 
}
