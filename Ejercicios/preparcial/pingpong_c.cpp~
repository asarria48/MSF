#include <iostream>
#include <cmath>
#include "../vector.h"
#include "../Random64.h"
using namespace std;


//Constantes del problema físico
const int Nx=1; 
const int Ny=1; 
const double Lx=20, Ly=60; 
const int N = Nx*Ny;
const int T = 5;
const int A = 1.0;
const double omega = 2*M_PI/T;
const double g = 9.80;
const double K=1.0e4;
const double Gamma = 10; 

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
  vector3D r,V,F; double m,R; double theta, omega;
public:
  void Inicie(double x0,double y0, double Vx0,double Vy0, double m0,double R0);
  void BorreFuerza(void){F.load(0,0,0);};
  void SumeFuerza(vector3D dF){F+=dF;};
  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; 
  double Gety(void){return r.y();}; 
  friend class Colisionador;
};

class Colisionador{
private:
public:
  void Inicie(void);
  void CalculeTodasLasFuerzas(Cuerpo * Granos , double dt, double t);
  void CalculeFuerzaEntre(Cuerpo & Granos1, Cuerpo & Granos2, double dt);
};

//-------Implementar las funciones de las clases------
//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0,double y0, double Vx0,double Vy0, double m0,double R0){
  r.load(x0,y0,0);  V.load(Vx0,Vy0,0); m=m0; R=R0;
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(coeficiente*dt/m); 
}
void Cuerpo::Dibujese(void){
  //Dibuja un circulo de radio R con centro en r y un segmento que parte de r 
  //con longitud R en la dirección theta
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
      <<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";
}

//------- Funciones de la clase Colisionador --------
/*void Colisionador::Inicie(void){
  for(int i=0;i<N;i++)
    for(int j=0;j<;j++)
      xCundall[i][j]=sold[i][j]=0;
      }*/
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Granos , double dt, double t){
  int i;
  //Borro las fuerzas de todos los Granoss
  for(i=0;i<N;i++)
    Granos[i]. BorreFuerza();

  // Fuerza de gravedad
  vector3D Fg;
  for(i=0;i<N;i++) {
    Fg.load(0,-Granos[i].m*g,0);  
    Granos[i].SumeFuerza(Fg);
  }
  

  //Recorro por parejas, calculo la fuerza de cada pareja y se la sumo a los dos
  for(i=0;i<N;i++){
    double altura = Granos[i].r.y()-A*sin(omega*t);
  
  //Determinar si hay colisión
  double R1 = Granos[i].R;
  double s = R1 - altura;

  if(s>0){
    
    double Fy = K*pow(s,1.5)-Gamma*Granos[i].m*Granos[i].V.y()*sqrt(s);

    vector3D F;
    if(Fy < 0){
      F.load(0,0,0);
      Granos[i].SumeFuerza(F);
    }

    else{
      F.load(0,Fy,0);
      Granos[i].SumeFuerza(F);
    }
    
  }

}
}

/*
void Colisionador::CalculeFuerzaEntre(Cuerpo & Granos1, Cuerpo & Granos2 , double & xCundall
                                        , double & sold, double dt){
  
  //Calcular el vector normal
  vector3D r21=Granos2.r-Granos1.r; 
  double d=r21.norm();
  double R1 = Granos1.R, R2 = Granos2.R;

  for(i = 0; i<N; i++){
    double altura = Granos[i].r.y()-A*sin(omega*t);
  
  //Determinar si hay colisión
  double s= R1 - altura;
  
  if(s>0){// Si hay colision

    //Vectores unitarios
    vector3D n = r21*(1.0/d),t; t.load(n.y(),-n.x(),0);

  

    // Calculo la velocidad de contacto
    vector3D Rw; Rw.load(0,0,R2*Granos2.omega + R1*Granos1.omega);
    vector3D Vc = Granos2.V - Granos1.V - (Rw^n);
    double Vcn = (Vc*n), Vct = Vc*t;

    // Calculo la fuerza

    double m0 = (Granos1.m * Granos2.m) / (Granos1.m + Granos2.m);
    double Fn = K*pow(s,1.5) - (Gamma*sqrt(s)*m0)*(Vc*n);
    vector3D F1 = n*Fn;

    //sumo la fuerza
    
     Granos2.SumeFuerza(F1,0);  Granos1.SumeFuerza(F1,0);
    
    /*
    double m0;

    //Calculo  la fuerza de deformacion plastica normal (Hertz-Kuratomo-Kano)
    double Fn =K*pow(s,1.5) - (Gamma*sqrt(s)*m0)*(Vc*n); 

    //Calculo la fuerza tangencial (Cundall)
    xCundall += Vct*dt; double Ft = 0; double Ftmax = mu*fabs(Fn);
    if(fabs(Ft)>Ftmax) Ft=Ft/fabs(Ft)*Ftmax;

    //Variables a calcular
    vector3D F1,F2, tau1,tau2;
    F2 = n*Fn + t*Ft; tau2 = ((n*(-R2))^F2); F1 = F2*(-1); tau1 = ((n*R1)^F1);
    //Calcula y Cargue las fuerzas

    Granos2.SumeFuerza(F2,tau2*k);  Granos1.SumeFuerza(F1,tau1*k);

  }

  if(sold>=0 && s<0) xCundall=0;
  sold=s;
}
}

*/


//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
  // cout<<"set terminal gif animate"<<endl; 
  // cout<<"set output 'Granoss.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:10]"<<endl;
  cout<<"set yrange[-10:60]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    //cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    //  cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    //cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    //cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
    
}
void TermineCuadro(void){
    cout<<endl;
}

int main(){
  
  Cuerpo Granos[N];
  Colisionador Hertz;
  Crandom ran64(1);
 //Parametros de la simulación
  int i;
  double m0=1.0; double R0=2.0; 
  //Variables auxiliares para la condición inicial
  double x0 = 0, y0 = 30, Vx0 = 0, Vy0 = 0;
  //Variables auxiliares para las paredes
  //double Rpared=1000*Lx, Mpared=1000*m0;
  //Variables auxiliares para correr la simulacion
  int Ncuadros=10000; double t,dt=1e-2,tmax=200,tcuadro=tmax/Ncuadros;
  double tdibujo;


   InicieAnimacion();
  
  //INICIO
  //Inicializar las paredes
  /*Granos[N].Inicie(Lx/2,Ly+Rpared, 0,  0, 0,  0,Mpared,Rpared); //Pared arriba
  Granos[N + 1].Inicie(Lx/2,-Rpared, 0,  0, 0,  0,Mpared,Rpared); //Pared abajo
  Granos[N + 2].Inicie(Lx+Rpared,Ly/2, 0,  0, 0,  0,Mpared,Rpared); //P<ared derecha
  Granos[N + 3].Inicie(-Rpared,Ly/2, 0,  0, 0,  0,Mpared,Rpared); //Pared izquierda*/
  


  /*  for(ix=0;ix<Nx;ix++)
    for(iy=0;iy<Ny;iy++){
      theta=2*M_PI*ran64.r();
      x0=(ix+1)*dx; y0=(iy+1)*dy; Vx0=V0*cos(theta); Vy0=V0*sin(theta);
      //----------------(x0,y0,z0,Vx0,Vy0,Vz0,m0,R0)
      Granos[iy*Nx+ix].Inicie( x0, y0, Vx0, Vy0, 0, 0,m0,R0);	
      }*/
  

  //---------------(x0,y0,Vx0,Vy0,m0,R0)
   Granos[0].Inicie(x0,y0,Vx0,Vy0,m0,R0);
  
  //CORRO
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){

    if(tdibujo>=tcuadro){

      //        cout << t << " " << Granos[0].Gety() << endl;
      
      InicieCuadro();
      for(i=0;i<N;i++) Granos[i].Dibujese();
      cout << " , " << -Lx/2 <<" + "<<Lx/7<<"*t,"<<sin(omega*t);
      TermineCuadro();
      
      tdibujo=0;
    }
    // cout<<Granos[1].Getx()<<" "<<Granos[1].Gety()<<endl;
    
    for(i=0;i<N;i++) Granos[i].Mueva_r(dt,xi);    
    Hertz.CalculeTodasLasFuerzas(Granos,dt,t); 
    for(i=0;i<N;i++) Granos[i].Mueva_V(dt,Um2lambdau2);

    for(i=0;i<N;i++) Granos[i].Mueva_r(dt,chi);
    Hertz.CalculeTodasLasFuerzas(Granos,dt,t); 
    for(i=0;i<N;i++) Granos[i].Mueva_V(dt,lambda);
    
    for(i=0;i<N;i++) Granos[i].Mueva_r(dt,Um2chiplusxi);
    Hertz.CalculeTodasLasFuerzas(Granos,dt,t); 
    for(i=0;i<N;i++)Granos[i].Mueva_V(dt,lambda);

    for(i=0;i<N;i++) Granos[i].Mueva_r(dt,chi);
    Hertz.CalculeTodasLasFuerzas(Granos,dt,t); 
    for(i=0;i<N;i++)Granos[i].Mueva_V(dt,Um2lambdau2);

    for(i=0;i<N;i++) Granos[i].Mueva_r(dt,xi);
    
  }
  return 0;

}
