#include <iostream>
#include <cmath>
#include "vector.h"
using namespace std;

//Constantes del problema físico
const double GM=1.0;
//Constantes del algoritmo de integración
const double xi=0.1786178958448091;
const double lambda=-0.2123418310626054;
const double chi=-0.06626458266981849;
const double Um2lambdau2=(1-2*lambda)/2;
const double Um2chiplusxi=1-2*(chi+xi);
//Deaclaración de la clase
class Cuerpo;

//Deaclaración de la interfase
class Cuerpo{
private:
  vector3D r,V,F; double m,R;
public:
  void Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,double m0,double R0);
  void CalculeFuerza(void);
  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; // Inline
  double Gety(void){return r.y();}; // Inline
};
//Implementación de las funciones
void Cuerpo::Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,double m0,double R0){
  r.load(x0,y0,z0);  V.load(Vx0,Vy0,Vz0); m=m0; R=R0;
}
void Cuerpo::CalculeFuerza(void){
  double aux=-GM*m*pow(r.norm2(),-1.5);
  F=r*aux;
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(coeficiente*dt/m);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}

//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'UnBalon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-11:11]"<<endl;
  cout<<"set yrange[-11:11]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<endl;
}


int main(){
  double r0=10,m0=1;
  double omega=sqrt(GM/(r0*r0*r0)); double T=2*M_PI/omega, V0=omega*r0;
  double t,dt=0.1,ttotal=1.1*T;
  int Ncuadros=2000; double tdibujo,tcuadro=ttotal/Ncuadros;
  Cuerpo Balon;

  //  InicieAnimacion();
  
  //----------(x0,y0,z0,Vx0,   Vy0,Vz0,m0,R0)
  Balon.Inicie(r0, 0, 0,  0,0.5*V0,  0, 1,0.5);
  for(t=tdibujo=0;t<ttotal;t+=dt,tdibujo+=dt){

    if(tdibujo>tcuadro){
      /*
      InicieCuadro();
      Balon.Dibujese();
      TermineCuadro();
      */
      tdibujo=0;
    }
    cout<<Balon.Getx()<<" "<<Balon.Gety()<<endl;
    
    Balon.Mueva_r(dt,xi);
    Balon.CalculeFuerza(); Balon.Mueva_V(dt,Um2lambdau2);
    Balon.Mueva_r(dt,chi);
    Balon.CalculeFuerza(); Balon.Mueva_V(dt,lambda);
    Balon.Mueva_r(dt,Um2chiplusxi);
    Balon.CalculeFuerza(); Balon.Mueva_V(dt,lambda);
    Balon.Mueva_r(dt,chi);
    Balon.CalculeFuerza(); Balon.Mueva_V(dt,Um2lambdau2);
    Balon.Mueva_r(dt,xi);
    
  }
  return 0;
}
