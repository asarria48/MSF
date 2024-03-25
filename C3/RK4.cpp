#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

double f(double t, double x){    //función a integrar
  return x/2;
}

void RK4(double & t, double & x, double dt){    //paso por referencia para no hacer copias de cada paso que involucre a t y a x

  //definición de RK4
  double dx1 = dt*f(t,x);
  double dx2 = dt*f(t + (dt/2), x + (dx1/2));
  double dx3 = dt*f(t + (dt/2), x + (dx2/2));
  double dx4 = dt*f(t+dt, x+dx3);

  x += (dx1 + 2*dx2 + 2*dx3 + dx4)/6;
  t += t*dt;

  //la función es void, por lo que hace cosas pero no retorna nada
}

int main(void){
  double dt = 0.1;     //paso de tiempo
  double t, x; 

  for(t = 0, x = 1; t < 10;){
    RK4(t,x,dt);
  }

  cout << "x("<<t<<")_approx="<<x<<endl;        //resultado numérico
  cout << "x("<<t<<")_exacto="<<exp(t/2)<<endl; //resultado analítico

  return 0;
}

