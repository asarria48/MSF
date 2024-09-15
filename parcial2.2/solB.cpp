#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

//---------Constantes globales---------

const int Lx = 256;  
const int Ly = 64;

const int Q = 5;  

const double W0 = 1.0/3.0;    
const double C0 = 0.24;         
const double C2 = C0*C0;       
const double AUX0 = 1-3*C2*(1-W0);
const double D = 0.95;

//Constantes de la evolución

const double tau = 0.5;       
const double Utau = 1.0/tau;  
const double UmUtau = 1.0 - Utau;  

//---------Clase Lattice Boltzmann---------

class LatticeBoltzman{
private:
  double w[Q];                
  int Vx[Q], Vy[Q];           
  double *f, *fnew;          
  
public:
  LatticeBoltzman(void);      
  ~LatticeBoltzman(void);     
  int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;}; //Función índice, me dice en qué índice está la componente fi de la celda ix iy, la formula hace la linealización pertinente
  double rho(int ix, int iy, bool UseNew);    
  double Jx(int ix, int iy, bool UseNew);     
  double Jy(int ix, int iy, bool UseNew);     
  double feq(double rho, double Jx0, double Jy0, int i);
  //Funciones de evolución temporal:
  void Start(double rho0, double Jx0, double Jy0);
  void Collision(void);
  void ImposeFields(int t);
  void Advection(void);
  void Print(const char * NameFile);
};


LatticeBoltzman::LatticeBoltzman(void){ 
  
  w[0] = W0; w[1] = w[2] = w[3] = w[4] = (1.0 - W0)/4;
 
  Vx[0] = 0; Vx[1] = 1; Vx[2] = 0; Vx[3] = -1; Vx[4] = 0;
  Vy[0] = 0; Vy[1] = 0; Vy[2] = 1; Vy[3] = 0;  Vy[4] = -1;

  int ArraySize = Lx*Ly*Q;    
  f = new double [ArraySize]; fnew = new double [ArraySize];  
}

LatticeBoltzman::~LatticeBoltzman(void){
  delete[] f; delete[] fnew;
}

//------Campos macroscópicos------

double LatticeBoltzman::rho(int ix, int iy, bool UseNew){  
  
  double sum; int i, n0;
  for(sum = 0, i = 0; i < Q; i++){    
    n0 = n(ix,iy,i);                  
    if(UseNew)                        
      sum += fnew[n0];              
    else
      sum += f[n0];
  }

  return sum;
}

double LatticeBoltzman::Jx(int ix, int iy, bool UseNew){  
  double sum; int i, n0;
  for(sum = 0, i = 0; i < Q; i++){
    n0 = n(ix,iy,i);
    if(UseNew)
      sum += Vx[i]*fnew[n0];         
    else
      sum += Vx[i]*f[n0];
  }

  return sum;
}

double LatticeBoltzman::Jy(int ix, int iy, bool UseNew){ 
  double sum; int i, n0;
  for(sum = 0, i = 0; i < Q; i++){
    n0 = n(ix,iy,i);
    if(UseNew)
      sum += Vy[i]*fnew[n0];
    else
      sum += Vy[i]*f[n0];
  }

  return sum;
}


double LatticeBoltzman::feq(double rho0, double Jx0, double Jy0, int i){
  if(i>0)                                         
    return 3*w[i]*(C2*rho0+Vx[i]*Jx0+Vy[i]*Jy0);  
  else
    return rho0*AUX0;                             
}

//------Funciones de evolución temporal------

void LatticeBoltzman::Start(double rho0, double Jx0, double Jy0){
  int ix, iy, i, n0;
  //Para cada celda
  for(ix = 0; ix < Lx; ix++)
    for(iy = 0; iy < Ly; iy++)
      //En cada dirección
      for(i = 0; i < Q; i++){
	n0 = n(ix,iy,i);
	f[n0] = feq(rho0,Jx0,Jy0,i); 
      }
}


void LatticeBoltzman::Collision(void){
//en este paso, f y fnew se actualizan para cada celda y cada dirección en función de las distribuciones actuales y la función de equilibrio
  int ix, iy, i, n0; double rho0, Jx0, Jy0;
  //Para cada celda
  for(ix = 0; ix < Lx; ix++){
    for(iy = 0; iy < Ly; iy++){
      //Computar los campos macroscópicos
      rho0 = rho(ix,iy,false); Jx0 = Jx(ix,iy,false); Jy0 = Jy(ix,iy,false);  //con los valores actuales de f y no los nuevos de fnew
      //Para cada dirección
      for(i = 0; i < Q; i++){
	n0 = n(ix,iy,i);  //para obtener el índice correspondiente en los arreglos f y fneW
	fnew[n0] = UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,i);
    }  
   }
  }

}

void LatticeBoltzman::ImposeFields(int t){  
  int ix, iy, i, n0; double rho0, Jx0, Jy0, lambda, omega, C, dx, dy,r2;
  double ixc = 64;
  double iyc = Ly/2;
  double R = 12;
  double c0 = 0.24;
  double c1 = 0.25;
  lambda = 256;
  omega = 2*M_PI/lambda*c0; 
  rho0 = sin(omega*t);
  
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      if(ix == 1){
	Jx0 = Jx(ix,iy,false); Jy0 = Jy(ix,iy,false);
	for(i = 0; i < Q; i++){
	n0 = n(ix,iy,i);
	fnew[n0] = feq(rho0,Jx0,Jy0,i);
      }
      }
    }
  }
 
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
     dx = (ix-ixc)*(ix-ixc);
     dy = (iy-iyc)*(iy-iyc);
     r2 = dx + dy;
     if(r2 < R*R){
       C = c0 - (c0-c1)*0.5*(1-tanh((r2-R*R)/4.0));
	Jx0 = Jx(ix,iy,false); Jy0 = Jy(ix,iy,false);
	for(i = 0; i < Q; i++){
	n0 = n(ix,iy,i);
	fnew[n0] = feq(0,Jx0,Jy0,i);
      }
      }
  }
    }
   
  }
    

void LatticeBoltzman::Advection(void) {
  int ix, iy, i, n0, ixnext, iynext, n0next;
  // Para cada celda
  for(ix = 0; ix < Lx; ix++) {
    for(iy = 0; iy < Ly; iy++) {
      // En cada dirección
      for(i = 0; i < Q; i++) {
        n0 = n(ix, iy, i);
	
        if (ix == 0 || ix == Lx-1) {
       
          ixnext = (ix == 0) ? ix + 1 : ix - 1;  
          iynext = (iy + Vy[i] + Ly) % Ly;      
          n0next = n(ixnext, iynext, (i+2)%4);
          fnew[n0] = D*f[n0next];  
        } else {
          // Advección normal en el resto del dominio
          ixnext = (ix + Vx[i] + Lx) % Lx; 
          iynext = (iy + Vy[i] + Ly) % Ly; 
          n0next = n(ixnext, iynext, i);
          f[n0next] = fnew[n0];
        }
      }
    }
  }
}
void LatticeBoltzman::Print(const char * NameFile){
  ofstream MyFile(NameFile);  
  double rho0;
  int ix, iy;
  for(ix = 0; ix < Lx; ix++){
    for(iy = 0; iy < Ly; iy++){
      rho0 = rho(ix,iy,true);
      MyFile<<ix<<" "<<iy<<" "<<rho0<<endl; //imprime la densidad en cada cuadro de la malla
    }
    MyFile<<endl;
  }
  MyFile.close();
}

//---------Funciones globales---------

int main(void){
  LatticeBoltzman Ondas;  //le doy un nombre a la clase
  double lambda = 256;
  double T = lambda/0.24;
  int t, tmax = int(10*T);  //defino un tiempo para la simulación
  int Niter = 2132+tmax;    //un periodo más 
  double rho0 = 1.0, Jx0 = 0, Jy0 = 0;  //inicializo los campos macroscópicos 
  double rhovalues[Lx]; //guardar los valores
  //Start
  Ondas.Start(rho0,Jx0,Jy0);  //inicialización independiente del tiempo
  //Run
  for(t=0;t<tmax;t++){  //comienzo el bucle temporal, los pasos que vienen dependen del tiempo porque deben hacerse en cada paso 
    Ondas.Collision();  
    Ondas.ImposeFields(t);
    Ondas.Advection();
  }
  for(t=tmax;t<Niter;t++){
    Ondas.Collision();  
    Ondas.ImposeFields(t);
    Ondas.Advection();
    for(int ix = 0; ix<Lx; ix++){
    rhovalues[ix] = Ondas.rho(ix,Ly,true);
   }
  }

  //Show
  Ondas.Print("Ondas.dat");

  ofstream rhodata("rhovsix.txt");
  for(int ix=0;ix<Lx;ix++){
    rhodata << ix << " " << rhovalues[ix] << endl;
  }

  rhodata.close();

  return 0;
}
