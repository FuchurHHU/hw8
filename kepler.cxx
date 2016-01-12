#include<cmath>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>

using namespace std;
void timestep(const double dt, double* p, double*q);
double Hamilton(double* p, double* q);
int main(){
  const int dim = 2;
  const double e = 0.6;
  double p[dim];
  double q[dim];
  const double tmin = 0;
  const double tmax = 20 * M_PI ;
  //const double dt = 0.05;
  const double dt = 0.0005;
  const int Nt = (tmax-tmin)/dt;
 //initial values:
  p[0] = 0; 
  p[1] = sqrt((1 + e)/(1 - e ));
  q[0] = 1 - e;
  q[1] = 0;
  double t = tmin;
  double H = Hamilton(p,q);
ofstream out("kepler.dat");
    out << t << "\t" << q[0] << "t" << q[1] <<"\t" << H << endl;
  
  
  for (int i = 0; i < Nt; i++ ){
    timestep(dt, p, q); 
    t+=dt;
    H = Hamilton(p,q);
    out << t << "\t" << q[0] << "\t" << q[1] <<"\t" << H << endl; 
 
  }
  out.close();
  return 0;
}

void timestep(const double dt, double* p, double*q){
  p[0] -= dt * pow(q[0]*q[0] + q[1]*q[1],-1.5) * q[0]; 
  p[1] -= dt * pow(q[0]*q[0] + q[1]*q[1],-1.5) * q[1];
  
  q[0] += dt*p[0];
  q[1] += dt*p[1]; //hier werden schon die überschriebenen ps genutzt (p_(n+1) aus Aufgabe)
  
  // neue p und q ausgerechnet aus Hamiltonian vgl Unterlagen  
}
double Hamilton(double* p, double* q){
   
  double H = 0.5 * (p[0]*p[0] + p[1]*p[1]) - pow(q[0]*q[0] + q[1]*q[1],-0.5);
  return H;
  
}
