#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::sin;

const double error = 1e-7;    //precisión requerida

double f(double x){           //función a la que le voy a hallar los ceros
  return sin(x)/x;
}

int main(){
  double a = 2, b = 4, fa;    //intervalo y función evaluada en a
  double m, fm;               //declaración del promedio y la función evaluada en él

  fa = f(a);

  while(b-a > error){         //mientras b-a sea mayor que la precisión requerida
    m = (b+a)/2;              //se hace el promedio entre los dos valores
    fm = f(m);                //y se evalúa la función ahí
    
    if(fa*fm > 0){            //si la función evaluada en a por la función evaluada en m es mayor que cero, el cero no está en el intervalo,
      fa = fm;                //tenemos que correr a y fa
      a = m;
    }

    else
      b = m;                  //si es menor que cero, el cero está en el intervalo, así que corro b a m y vuelvo a iterar.
  }
  
  cout << "El cero es = " << (b+a)/2 << endl;
  return 0;
}
