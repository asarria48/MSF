#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

const double GM=1.0;

//Deaclaración de la clase
class Cuerpo;

//Deaclaración de la interfase
class Cuerpo{
private:
  vector3D r,V,F; double m;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
	      double m0);
  void CalculeFuerza(void);
  void Muevase(double dt);
  double Getx(void){return r.x();}; // Inline
  double Gety(void){return r.y();}; // Inline
};
//Implementación de las funciones
void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
	      double m0){
  r.load(x0,y0,z0); V.load(Vx0,Vy0,Vz0); m=m0;
}
void Cuerpo::CalculeFuerza(void){
  double norm = sqrt(r.x()*r.x() + r.y()*r.y() + r.z()*r.z());
  double aux = -GM*m/pow(norm,3);
  F.load(aux*r.x(),aux*r.y(),aux*r.z());
}
void Cuerpo::Muevase(double dt){
  //Algoritmo de Euler
  r+=V*dt;  V+=F*(dt/m);
}

//----------- Funciones Globales -----------


int main(){
  double t,dt=0.001,ttotal=200;
  int Ncuadros=200; double tdibujo,tcuadro=ttotal/Ncuadros;
  int r0 = 10;
  double omega = sqrt(GM/(pow(r0,3)));
  double T = 2*M_PI/omega;
  double V0 = omega*r0;
  Cuerpo Balon;
  
  //----------(x0,y0,z0,Vx0,Vy0,Vz0,m0   ,R0)
  Balon.Inicie( r0, 0, 0, 0, 0.5*V0,  0,0.453);
  for(t=tdibujo=0;t<ttotal;t+=dt,tdibujo+=dt){
    if(tdibujo>tcuadro){

      tdibujo=0;
    }
    cout<<Balon.Getx()<<" "<<Balon.Gety()<<endl;
    Balon.CalculeFuerza();
    Balon.Muevase(dt);

  }
  return 0;
}
