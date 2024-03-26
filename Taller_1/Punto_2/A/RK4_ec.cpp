#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

//Tengo un sistema de dos ecuaciones acopladas
double f1(double r, double R, double y){
  return y;
}

double f2(double r, double R, double y){
  return -(y/r)-R;
}


void RK4acop(double & r, double & R, double & y, double dr){  

double dR1, dR2, dR3, dR4, dy1, dy2, dy3, dy4;

  dR1 = dr*f1(r,R,y);                    dy1 = dr*f2(r,R,y);
  dR2 = dr*f1(r+dr/2,R+dR1/2,y+dy1/2);   dy2 = dr*f2(r+dr/2,R+dR1/2,y+dy1/2);
  dR3 = dr*f1(r+dr/2,R+dR2/2,y+dy2/2);   dy3 = dr*f2(r+dr/2,R+dR2/2,y+dy2/2);
  dR4 = dr*f1(r+dr,R+dR3,y+dy3);         dy4 = dr*f2(r+dr,R+dR3,y+dy3);

  R += (dR1+2*dR2+2*dR3+dR4)/6;
  y += (dy1+2*dy2+2*dy3+dy4)/6;

  r += dr;
  
}


int main(void){

  double r, R, y;
  double dr = 0.01;

  for(r = 0.01, R = 1.0, y = 0.0; r <= 10;){

    cout << r << " " << R << endl;
    RK4acop(r,R,y,dr);
 
  }
  
  return 0;
}

