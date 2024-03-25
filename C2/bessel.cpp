#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::cos;
using std::sin;

double f(double t, double x, double n){  //función a integrar
  return cos(n*t - x*sin(t));             //las funciones de Bessel pueden escribirse de esta forma, pero falta integrar
}

double simpson(double a, double b, double nsteps, double x, double n){    //método de integración
  nsteps = 2*nsteps;                   //el número de intervalos debe ser par
  double h = (b-a)/nsteps;        //tamaño del intervalo 
  double sum = 0;

  for(int i = 0; i<=nsteps; i++){
    double t = a + i*h;      //puntos en los que se va a evaluar la función
    if(i==0 or i==nsteps){        //para el primer y último valor
      sum += f(t,x,n);         
      }
    else if(i%2 == 0){       //para los valores pares 
      sum += 2*f(t,x,n);
    }
    else{                    //para los que sobran (impares)
      sum += 4*f(t,x,n);
    }
  }
  return sum*h/3;            //cojo toda la suma y la multiplico por h/3
}

double bessel(double n, double x){    //función de Bessel completa, integrada por método de Simpson
  return (1/M_PI*simpson(0,M_PI,100,x,n));
}

int main(){

  double n = 0;              //orden de la función de Bessel
  double a = 0;
  double b = M_PI;
  double nsteps = 9999;
  double h = (b-a)/nsteps;

  for(int i = 0; i < nsteps; i++){
    double x = a + i*h;      //losx puntos a evaluar en la función de Bessel

    cout << x << " " << bessel(n,x) << endl;
  }
    
  return 0;
}





