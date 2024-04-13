#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "vector.h"
#include "Random64.h"
using namespace std;

//Constantes del problema físico
double Lx=200, Ly=200; 
const int Nx=6, Ny=1;
const int N=21, Ntot=N+4;
const double g=9.8, KHertz=1.0e4, Gamma=0.5, Kcundall=500, mu=0.4;
const double Kcentral = 1;
const double Krepulsion = 1000;
const double R_anillo = 40;
const double V = 1; //magnitud de la velocidad inicial para el punto 4


//Constantes del algoritmo de integración
const double xi=0.1786178958448091;
const double lambda=-0.2123418310626054;
const double chi=-0.06626458266981849;
const double Um2lambdau2=(1-2*lambda)/2;
const double Um2chiplusxi=1-2*(chi+xi);

//--------------- Declarar las clases-----------
class Cuerpo;
class Colisionador;

//--------- Declarar las interfases de las clases---------
class Cuerpo{
private:
  vector3D r,V,F; double m,R,Q; double theta,omega,tau,I;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,
	      double theta0,double omega0,double m0,double R0,double Q0);
  void BorreFuerza(void){F.load(0,0,0); tau=0;};// Inline
  void SumeFuerza(vector3D dF,double dtau){F+=dF; tau+=dtau;};// Inline
  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; // Inline
  double Gety(void){return r.y();}; // Inline
  friend class Colisionador;
};
class Colisionador{
private:
  double xCundall[Ntot][Ntot],sold[Ntot][Ntot];
public:
  void Inicie(void);
  void CalculeTodasLasFuerzas(Cuerpo * Grano,double dt);
  void CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2,
			  double & xCundall,double & sold,double dt);
};

//-------Implementar las funciones de las clases------
//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,
		    double theta0,double omega0,double m0,double R0,double Q0){
  r.load(x0,y0,0);  V.load(Vx0,Vy0,0); m=m0; R=R0; Q=Q0;
  theta=theta0; omega=omega0; I=2.0/5.0*m*R*R;
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(coeficiente*dt); theta+=omega*(coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(coeficiente*dt/m); omega+=tau*(coeficiente*dt/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
      <<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";
}
//------- Funciones de la clase Colisionador --------
/*void Colisionador::Inicie(void){
  int i,j; //j>i
  for(i=0;i<Ntot;i++)
    for(j=0;j<Ntot;j++)
      xCundall[i][j]=sold[i][j]=0;
      }*/
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Grano,double dt){
  int i,j;
  //Borro las fuerzas de todos los granos
  for(i=0;i<N+4;i++)
    Grano[i].BorreFuerza();
  //sumo fuerza central
  // vector3D Fc;
  for(i=0;i<N;i++){
    vector3D Fcentral = -Kcentral*Grano[i].r;
    //Fc.load(Fcentral,Fcentral,0);
    Grano[i].SumeFuerza(Fcentral,0);
  }
  for(i=0;i<N;i++){
    vector3D Fdisipa = -Gamma*Grano[i].V;
    //Fc.load(Fcentral,Fcentral,0);
    Grano[i].SumeFuerza(Fdisipa,0);
  }
  
  //Recorro por parejas, calculo la fuerza de cada pareja y se la sumo a los dos
  for(i=0;i<N;i++)
    for(j=i+1;j<N+4;j++)
      CalculeFuerzaEntre(Grano[i],Grano[j],xCundall[i][j],sold[i][j],dt);
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2,
				      double & xCundall,double & sold,double dt){
  //Determinar si hay colision
  vector3D r21=Grano2.r-Grano1.r; double d=r21.norm();
  double R1=Grano1.R,R2=Grano2.R;
  double s=(R1+R2)-d;
  if(d>1e-10){ //Si hay colisión

    //Vectores unitarios
    vector3D n=r21*(1.0/d),t,k; t.load(n.y(),-n.x(),0); k.load(0,0,1);

    double Frepul = Krepulsion*Grano1.Q*Grano2.Q/(d*d);
    vector3D F_repulsiva = n*Frepul;

    Grano1.SumeFuerza(F_repulsiva*(-1),0);
    Grano2.SumeFuerza(F_repulsiva,0);
}
}
//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'UnBalon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-110:"<<110<<"]"<<endl;
  cout<<"set yrange[-110:"<<110<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
    cout<<endl;
}

int main(){
  Cuerpo Grano[N];
  Colisionador Hertz;
  Crandom ran64(1);
  int i,ix,iy;
  //Parametros de la simulación
  double m0=1.0; double R0=3.0; double Q0 = 1.0;
  double T = 2*M_PI*sqrt(m0/Kcentral);
  double kT=10; 
  //Variables auxiliares para la condición inicial
  double dx=Lx/(Nx+1),dy=Ly/(Ny+1);
  double theta; double V0=0, omega0=3*M_PI;//sqrt(kT/m0);
  double x0=40,y0=0,Vx0,Vy0;
  //Variables auxiliares para las paredes
  double Rpared=100*Lx, Mpared=100*m0;
  //Variables auxiliares para correr la simulacion
  int Ncuadros=10000; double t,tdibujo,dt=T/1000,tmax=5*T,tcuadro=tmax/Ncuadros;

  vector<double> r_maxorigen;
  vector<int> Nparticulas = {6,9,14,21,32,48};

  for(int N : Nparticulas) {
  
    //     InicieAnimacion();
  for(i = 0; i<N; i++){
    double angulo = (2*M_PI/N)*i;
    double angulovel = (2*M_PI*ran64.r());
    x0 = R_anillo*cos(angulo);
    y0 = R_anillo*sin(angulo);
    Vx0 = V*cos(angulovel);
    Vy0 = V*sin(angulovel);
    Grano[i].Inicie(x0,y0,Vx0,Vy0,0,0,m0,R0,Q0);
  }
      
  //CORRO
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){

    if(tdibujo>tcuadro){

      //   cout << t << " " << Grano[0].Getx() << endl;
      
      /*     InicieCuadro();
      for(i=0;i<N;i++) Grano[i].Dibujese();
      TermineCuadro();*/
      
      tdibujo=0;
    }
    //cout<<Grano[1].Getx()<<" "<<Grano[1].Gety()<<endl;
    
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,xi);    
    Hertz.CalculeTodasLasFuerzas(Grano,dt); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeTodasLasFuerzas(Grano,dt); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,Um2chiplusxi);
    Hertz.CalculeTodasLasFuerzas(Grano,dt); for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeTodasLasFuerzas(Grano,dt); for(i=0;i<N;i++)Grano[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,xi);
       
  }

  double r_max = 0;
  for(int i = 0; i<N; i++){
    double r = sqrt(Grano[i].Getx()*Grano[i].Getx() + Grano[i].Gety()*Grano[i].Gety());
    if(r>r_max)
      r_max = r;
  }
  r_maxorigen.push_back(r_max);
  }

ofstream data_file("punto5grafica.dat");
for(size_t i = 0; i < Nparticulas.size(); i++){
    data_file << Nparticulas[i] << " " << r_maxorigen[i] << endl;
  }
  data_file.close();
  
  return 0;
}
