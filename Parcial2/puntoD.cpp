#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
using namespace std;

//---------Constantes globales---------

//Defino un dominio de simulación

const int Lx = 256;    //tamaño de la grilla en 2D
const int Ly = 64;
const int N = 24; //divisiones cilindro
const double R = 12.0; //radio del cilindro
const double ixc = 64;
const double iyc = Ly/2;

const int Q = 5;              //Número de velocidades, 5 para 2D: arriba, abajo, izquierda, derecha, centro

const double W0 = 1.0/3.0;    //Peso de la velocidad central
//const double C = 0.24;         //Velocidad del sonido, C < 1 por criterio de Courant
//const double C2 = C*C;        //velocidad del sonido al cuadrado 
//const double AUX0 = 1-3*C2*(1-W0);    //Variable auxiliar

//Constantes de la evolución

const double tau = 0.5;       //Tiempo de relajación, para alcanzar f_equilibrio 
const double Utau = 1.0/tau;  //para integración
const double UmUtau = 1.0 - Utau;  //para integración 

//---------Clase Lattice Boltzmann---------

class LatticeBoltzman{
private:
  double w[Q];                //Pesos de las velocidades
  int Vx[Q], Vy[Q];           //Vectores de las velocidades
  double *f, *fnew;           //Datos para todas las celdas, funciones de distribución
  
public:
  LatticeBoltzman(void);      //Constructor para definir los pesos
  ~LatticeBoltzman(void);     //Destructor
  int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;}; //Función índice, me dice en qué índice está la componente fi de la celda ix iy, la formula hace la linealización pertinente
  double rho(int ix, int iy, bool UseNew);    //Campo macroscópico rho
  double Jx(int ix, int iy, bool UseNew);     //Campo macroscópico Jx
  double Jy(int ix, int iy, bool UseNew);     //Campo macroscópico Jy
  double feq(double rho, double Jx0, double Jy0, int i, int ix, int iy); //Función de equilibrio a partir de los campos y de la dirección i
  double C(int ix, int iy);
  //Funciones de evolución temporal:
  void Start(double rho0, double Jx0, double Jy0);
  void Collision(void);
  void ImposeFields(int t);
  void Advection(void);
  void BounceBack(void);
  void Print(const char * NameFile);

  //interpolaciones
  double rhointerpolado(double x, double y);
  double jxinterpolado(double x, double y);
  double jyinterpolado(double x, double y);
};

  void calculardS(double &dSx, double &dSy, int n);
  double calcularFx(LatticeBoltzman &ondas, int n);

//---------Constructor---------
//El constructor se llama cada vez que creo un ejemplar de la clase a la que pertenece

//EL constructor inicializa los pesos y las velocidades, además reserva memoria para guardar las funciones de distribución 

LatticeBoltzman::LatticeBoltzman(void){ 
  //Definir los pesos: los pesos definen la influencia de cada dirección de velocidad en la dinámica del sistema
  w[0] = W0; w[1] = w[2] = w[3] = w[4] = (1.0 - W0)/4;
  //Definir los vectores de velocidad
  Vx[0] = 0; Vx[1] = 1; Vx[2] = 0; Vx[3] = -1; Vx[4] = 0;
  Vy[0] = 0; Vy[1] = 0; Vy[2] = 1; Vy[3] = 0;  Vy[4] = -1;
  //Crear los arreglos dinámicos
  //Son dinámicos para crearse en el tiempo de ejecución, adaptándose al tamaño adecuado que depende de la malla y el número de direcciones 
  int ArraySize = Lx*Ly*Q;     //*Q porque tengo Q funciones en cada celda
  f = new double [ArraySize]; fnew = new double [ArraySize];  //el operador new llos hace dinámicos 
}

//---------Destructor---------
//El destructor se llama cada vez que quiero borrar cuando no estoy usando espacio en la memoria
//libera la memoria reservada para que no se acabe el almacenamiento y haya un óptimo uso de la memoria
LatticeBoltzman::~LatticeBoltzman(void){
  delete[] f; delete[] fnew;
}

//------Campos macroscópicos------

//magnitudes físicas que describen el comportamiento del fluido a escala macroscópica
//estos campos se calculan a partir de las funciones de distribución

double LatticeBoltzman::rho(int ix, int iy, bool UseNew){  //calcula la densidad en la celda bidimensional

  //la densidad se calcula sumando las funciones de distribución en cada dirección i  
  
  double sum; int i, n0;
  for(sum = 0, i = 0; i < Q; i++){    //se inicializa la suma en cero, se itera sobre cada dirección
    n0 = n(ix,iy,i);                  //usa la función n para encocntrar los índices en cada celda ix iy para cada dirección i
    if(UseNew)                        //para usar fnew o f
      sum += fnew[n0];                //usa n0 para acceder a los valores en fnew y en f
    else
      sum += f[n0];
  }

  return sum;
}

double LatticeBoltzman::Jx(int ix, int iy, bool UseNew){  //densidad de flujo en dirección x
  double sum; int i, n0;
  for(sum = 0, i = 0; i < Q; i++){
    n0 = n(ix,iy,i);
    if(UseNew)
      sum += Vx[i]*fnew[n0];          //se suma la componente de la velocidad en x por la función de distribución correspondiente
    else
      sum += Vx[i]*f[n0];
  }

  return sum;
}

double LatticeBoltzman::Jy(int ix, int iy, bool UseNew){  //densidad de flujo en dirección y
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
double LatticeBoltzman::C(int ix, int iy){
  int x0 = 64;
  int y0 = Ly/2;
  double R = 12.0;
  double c0 = 0.24;
  double c1 = 0.25;
  double d = 4.0;

  double r2 = (ix - x0)*(ix - x0) + (iy - y0)*(iy - y0);
  return c0 - ((c0-c1)/2)*(1-tanh((r2-R*R)/d));
}

double LatticeBoltzman::feq(double rho0, double Jx0, double Jy0, int i, int ix, int iy){  //dada la densidad y las corrientes, calcula la función de equilibrio dada una dirección i
  double C2 = C(ix,iy)*C(ix,iy);
  double AUX0 = 1-3*C2*(1-W0);
  if(i>0)                                         //para una dirección que no sea la central
    return 3*w[i]*(C2*rho0+Vx[i]*Jx0+Vy[i]*Jy0);  //asocio la función de equilibrio en esa dirección usando los campos macroscópicos
  else
    return rho0*AUX0;                             //si es la dirección central, aplico esta otra forma 
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
	f[n0] = feq(rho0,Jx0,Jy0,i,ix,iy); //Obligo que el valor de la función sea el valor de equilibrio que se obtiene de los campos macroscópicos, inicializo no a cero sino a los valores de equilibrio predeterminados 
      }
}

void LatticeBoltzman:: BounceBack(void){
  int n0, n0_next;
  double D = 0.95; //factor de amortiguamiento

  //A la izquierda ix = 0:
  for(int iy = 0; iy<Ly;iy++){
    for(int i = 0; i<4; i++){
      n0 = n(0,iy,i);
      n0_next = n(0,iy,(i+2)%4);
      fnew[n0] = D*f[n0_next];
    }
  }

  //A la derecha ix = 127:
    for(int iy = 0; iy<Ly;iy++){
    for(int i = 0; i<4; i++){
      n0 = n(Lx-1,iy,i);
      n0_next = n(Lx-1,iy,(i+2)%4);
      fnew[n0] = D*f[n0_next];
    }
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
	n0 = n(ix,iy,i);  //para obtener el índice correspondiente en los arreglos f y fnew
	fnew[n0] = UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,i,ix,iy); //Ecuación de Boltzmann discretizada (fórmula que nos dan)
      }
   }
 }
  
  BounceBack();
}

void LatticeBoltzman::ImposeFields(int t){  //impone un campo que afecte el fluido, por ejemplo: una oscilación 
  int ix, iy, i, n0; double rho0, Jx0, Jy0, lambda, omega, c0;

  //defino los parámetros asociados al campo que voy a imponer
  lambda = 256; //como manda el ejercicio
  omega = 2*M_PI/lambda;  //velocidad angular de la onda 
  for(iy = 0; iy < Ly; iy++){
    ix = 1;
    rho0 = sin(omega*t); Jx0 = Jx(ix,iy,false); Jy0 = Jy(ix,iy,false);
    for(i = 0; i < Q; i++){
	n0 = n(ix,iy,i);
	fnew[n0] = feq(rho0,Jx0,Jy0,i,ix,iy);
    }
  }
}

void LatticeBoltzman::Advection(void){
  int ix, iy, i, n0, ixnext, iynext, n0next;
  //Para cada celda
   for(ix = 0; ix < Lx; ix++)
    for(iy = 0; iy < Ly; iy++)
      //En cada dirección
      for(i = 0; i < Q; i++){
        //Calculo la posición actual y la siguiente
	//Con condiciones de frontera periodicas
	ixnext = (ix+Vx[i]+Lx)%Lx; iynext = (iy+Vy[i]+Ly)%Ly; //por ejemplo, si la malla fuese de 10x10, y para i = 0, ixnext = ix porque no hay velocidad ni en x ni en y
	n0 = n(ix,iy,i); n0next = n(ixnext,iynext,i);
	//Hago la advección: paso de los fnew a los f
	f[n0next] = fnew[n0];  //para luego usar los actuales (f) para hacer los nuevos (fnew)
      }
   }

double LatticeBoltzman::rhointerpolado(double x, double y){
  int ix = static_cast<int>(x);
  int iy = static_cast<int>(y);
  double u = x - ix;
  double v = y - iy;

  //Para (ix,iy)
  double rho_ix_iy = rho(ix,iy,false);

  //Para (ix+1,iy)
  double rho_ixP1_iy = rho(ix+1,iy,false);

  //Para (ix,iy+1)
  double rho_ix_iyP1 = rho(ix,iy+1,false);

  //Para (ix+1,iy+1)
  double rho_ixP1_iyP1 = rho(ix+1,iy+1,false);

  //Total para rho:

  double interpolatedrho = rho_ix_iy*(1-u)*(1-v)+rho_ixP1_iy*(u)*(1-v)+rho_ix_iyP1*(1-u)*(v)+rho_ixP1_iyP1*u*v;

  return interpolatedrho;
}

double LatticeBoltzman::jxinterpolado(double x, double y){
  int ix = static_cast<int>(x);
  int iy = static_cast<int>(y);
  double u = x - ix;
  double v = y - iy;

  //Para (ix,iy)
  double jx_ix_iy = Jx(ix,iy,false);

  //Para (ix+1,iy)
  double jx_ixP1_iy = Jx(ix+1,iy,false);

  //Para (ix,iy+1)
  double jx_ix_iyP1 = Jx(ix,iy+1,false);

  //Para (ix+1,iy+1)
  double jx_ixP1_iyP1 = Jx(ix+1,iy+1,false);

  //Total para Jx:

  double interpolatedjx = jx_ix_iy*(1-u)*(1-v)+jx_ixP1_iy*(u)*(1-v)+jx_ix_iyP1*(1-u)*(v)+jx_ixP1_iyP1*u*v;

  return interpolatedjx;
}


double LatticeBoltzman::jyinterpolado(double x, double y){
  int ix = static_cast<int>(x);
  int iy = static_cast<int>(y);
  double u = x - ix;
  double v = y - iy;

  //Para (ix,iy)
  double jy_ix_iy = Jy(ix,iy,false);

  //Para (ix+1,iy)
  double jy_ixP1_iy = Jy(ix+1,iy,false);

  //Para (ix,iy+1)
  double jy_ix_iyP1 = Jy(ix,iy+1,false);

  //Para (ix+1,iy+1)
  double jy_ixP1_iyP1 = Jy(ix+1,iy+1,false);

  //Total para Jy:

  double interpolatedjy = jy_ix_iy*(1-u)*(1-v)+jy_ixP1_iy*(u)*(1-v)+jy_ix_iyP1*(1-u)*(v)+jy_ixP1_iyP1*u*v;

  return interpolatedjy;
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

void calculardS(double &dSx, double &dSy, int n){
  double theta = 2*M_PI*n/N; //ángulo del elemento
  dSx = R*cos(theta);
  dSy = R*sin(theta);
}

double calcularFx(LatticeBoltzman &ondas, int n){
  double theta = 2*M_PI*n/N; //ángulo del elemento
  double x = ixc + R*cos(theta);
  double y = iyc + R*sin(theta);
  double c0 = 0.24;

  double rhoP = ondas.rhointerpolado(x,y);
  double JxP = ondas.jxinterpolado(x,y);
  double JyP = ondas.jyinterpolado(x,y);

  double dSx, dSy;
  calculardS(dSx,dSy,n);

  double dFx = (-0.5*(JxP*JxP + JyP*JyP)+ 0.5 * rhoP * rhoP * c0 * c0)*dSx + JxP * (JxP * dSx + JyP * dSy);

  return dFx;
}

//---------Funciones globales---------

int main(void){
  LatticeBoltzman Ondas;  //le doy un nombre a la clase 
  int t, tmax = 20000;  //defino un tiempo para la simulación
  double rho0 = 1.0, Jx0 = 0, Jy0 = 0;  //inicializo los campos macroscópicos
  //  double rho_val[Lx];
  double T = 1066.6; //lambda/c0
  int N_iter = int(10.5*T);

  ofstream Fxvst("fxvst.txt");
  //Start
  Ondas.Start(rho0,Jx0,Jy0);  //inicialización independiente del tiempo
  //Run
  //Oscilar por 10.5 T sin tomar datos
  for(t=0;t<N_iter;t++){  //comienzo el bucle temporal, los pasos que vienen dependen del tiempo porque deben hacerse en cada paso 
    Ondas.Collision();  
    Ondas.ImposeFields(t);
    Ondas.Advection();
  }
  //Empiezo a tomar datos después de 10.5 T
  for(t=N_iter;t<tmax;t++){
    Ondas.Collision();  
    Ondas.ImposeFields(t);
    Ondas.Advection();

    double Fxtotal = 0.0;
     
    for(int n = 0; n<N; n++){
      double dFx = calcularFx(Ondas,n);
      Fxtotal += dFx;
    }

    Fxvst << t << " " << Fxtotal << endl;
    
  }
 
  Fxvst.close();

  return 0;
}

//IMPRIMIR EN CONSOLA
//g++
//gnuplot
//splot 'Ondas.dat' w l
