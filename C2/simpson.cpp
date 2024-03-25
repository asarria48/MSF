#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::cos;

double f(double x){          //función a integrar
  return cos(x);
}

double simpson(double a, double b, double n){    //método de integración
  n = 2*n;                   //el número de intervalos debe ser par
  double h = (b-a)/n;        //tamaño del intervalo 
  double sum = 0;

  for(int i = 0; i<=n; i++){
    double x = a + i*h;      //puntos en los que se va a evaluar la función
    if(i==0 or i==n){        //para el primer y último valor
	sum += f(x);         
      }
    else if(i%2 == 0){       //para los valores pares 
      sum += 2*f(x);
    }
    else{                    //para los que sobran (impares)
      sum += 4*f(x);
    }
  }
  return sum*h/3;            //cojo toda la suma y la multiplico por h/3
}

int main(){

  double a = 0;
  double b = M_PI/2;
  double n = 50;

  cout << "La integral por Simpson da " << simpson(a,b,n) << endl;
    
  return 0;
}





