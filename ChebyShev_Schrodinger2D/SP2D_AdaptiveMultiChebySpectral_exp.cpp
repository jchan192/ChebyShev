#include <complex>
#include <vector>
#include <iostream>
#include <cstring>
#include <fstream>
#include <algorithm>
#include "headers_Cheby/SP1D.h"
#include <random>
#include <functional>
#include <chrono>
#include <fftw3.h>
#include <iomanip>

using complex = std::complex<double>; 
using vector_real = std::vector<double>;
using vector_complex = std::vector<std::complex<double>>;
using vector_patch = std::vector<Patches>;

complex linear_interp(double x, double x0, double x1,
		      complex psi1, complex psi2)
{

  double xd = (x-x0)/(x1-x0);
  return psi1*(1-xd) + psi2 * xd;
}

vector_real cheby_exgrid(double N, double a, double b)           
{
  vector_real x(N+1);
  for (int i=0;i<N+1;i++)
    x[i] = (b-a)/2.*(-cos(i*M_PI/N))+(b+a)/2.0;

  return x;
}
double cheby_poly(double x, double n)
{
  double T0 = 1.0;
  double T1 = x;
  double T;

  if (n == 0) return  T0;
  else if (n == 1) return  T1;
  else if (n > 1){
    for (int i=0; i<n-1; i++){
      T = 2*x*T1-T0;
      T0 = T1;
      T1 = T;
    }
    return T;
  }
  return 0;
}

vector_real cheby_exgrid_oni(double N, double a, double b, int ith)
{
  // cheby grid but locating ith point at a & N-ith at b

  vector_real x(N+1);

  double scale = (a-b) / (cos((N-ith)*M_PI/N) - cos(ith*M_PI/N)) ; 
  double mid   = b - scale * (-cos((N-ith)*M_PI/N));      

  for (int i=0;i<N+1;i++)
    x[i] = scale *(-cos(i*M_PI/N)) + mid ;

  return x;
}

//void cheby_coeff(vector_complex y, double N,vector_complex &Coeff)
void cheby_coeff()
{
  
  int N =32;
  vector_real coeff((N+1)*(N+1));

  vector_real x = cheby_exgrid(N,0,1);

  for (int i=0;i<N+1 ;i++){
  for (int j=0;j<N+1 ;j++){
    coeff[i*(N+1)+j] = x[i]*x[i] + x[j]*x[j];
//    std::cout << y[i] <<"\n";
  }}
  fftw_plan plan = fftw_plan_r2r_2d(N+1,N+1, &coeff[0], &coeff[0], FFTW_REDFT00, FFTW_REDFT00,FFTW_ESTIMATE);  
  fftw_execute(plan);

  // re-scale coefficients
  for (int j=0;j<(N+1)*(N+1) ;j++)
    coeff[j] /= N*N;

  for (int j=0;j<(N+1);j++){
    coeff[0*(N+1)+j] *= 0.5;
    coeff[j*(N+1)+0] *= 0.5;
  }
  for (int j=1;j<(N+1);j++){
    coeff[N*(N+1)+j] *= 0.5;
    coeff[j*(N+1)+N] *= 0.5;
  }
  std::cout << "\n";
  for (int i=0;i<N+1 ;i++){
  for (int j=0;j<N+1 ;j++){
    std::cout << coeff[i*(N+1)+j] <<"\n";
  }}

  double interp;
  double x_chebgrid,y_chebgrid;
  double a = 0.;
  double b = 1.;
  
  x_chebgrid = -( 2.0 * x[1] - a  - b ) / ( b-a);
  y_chebgrid = -( 2.0 * x[1] - a  - b ) / ( b-a);

    // interp  = coeff[0]*cheby_poly(x_chebgrid,0) * 0.5;

    // for (int j=1;j<N;j++)
    //   interp += coeff[j]*cheby_poly(x_chebgrid,j);

    // interp += coeff[N]*cheby_poly(x_chebgrid,N) * 0.5;

  interp = 0.0 ;
  for (int i=0;i<N+1;i++)
  for (int j=0;j<N+1;j++)
      interp += coeff[i*(N+1)+j]*cheby_poly(x_chebgrid,j)*cheby_poly(y_chebgrid,i);

 std::cout << interp <<"\n";
}
void cheby_coeff(vector_complex &y, double N, vector_complex &Coeff,Patches &P)
{
  
  vector_real coeff_R((N+1)*(N+1));
  vector_real coeff_I((N+1)*(N+1));

  vector_real x = cheby_exgrid(N,0,1);

  for (int j=0;j<(N+1)*(N+1);j++){
    P.coeff_R[j] = y[j].real();
    P.coeff_I[j] = y[j].imag();
  }  
  // fftw_plan plan_R = fftw_plan_r2r_2d(N+1,N+1, &coeff_R[0], &coeff_R[0], FFTW_REDFT00, FFTW_REDFT00,FFTW_ESTIMATE);  
  // fftw_plan plan_I = fftw_plan_r2r_2d(N+1,N+1, &coeff_I[0], &coeff_I[0], FFTW_REDFT00, FFTW_REDFT00,FFTW_ESTIMATE);  
  fftw_execute(P.plan_R);
  fftw_execute(P.plan_I);
  // fftw_destroy_plan(plan_R);
  // fftw_destroy_plan(plan_I);

  for (int j=0;j<(N+1)*(N+1);j++){
    Coeff[j].real(P.coeff_R[j]) ;
    Coeff[j].imag(P.coeff_I[j]) ;
  }
  
  // re-scale coefficients
 for (int j=0;j<(N+1)*(N+1) ;j++)
   Coeff[j] /= N*N;

  for (int j=0;j<(N+1);j++){
    Coeff[0*(N+1)+j] *= 0.5;
    Coeff[j*(N+1)+0] *= 0.5;
  }
  for (int j=1;j<N+1;j++){
    Coeff[N*(N+1)+j] *= 0.5;
    Coeff[j*(N+1)+N] *= 0.5;
  }

//  Coeff[N*(N+1)+N] *= 0.25;

 //  double x_chebgrid = -( 2.0 * x[1] - 0  -1 ) / ( 1);
 //  double y_chebgrid = -( 2.0 * x[1] - 0  -1 ) / ( 1);

 //  double interp = 0.0 ;
 //  for (int i=0;i<N+1;i++)
 //  for (int j=0;j<N+1;j++)
 //      interp += Coeff[i*(N+1)+j].real()*cheby_poly(x_chebgrid,j)*cheby_poly(y_chebgrid,i);

 // std::cout << interp <<"\n";
 // std::exit(0);
}

vector_complex cheby_coeff(vector_real x,vector_complex y, double N)
{

  vector_complex a;
  complex sum;
   vector_real x_chebgrid = cheby_exgrid(N,-1.0,1.0);
   for (int j=0;j<N+1;j++){
     sum = y[0]*cheby_poly(x_chebgrid[0],j) * 0.5; 

     for (int j_sum=1;j_sum<N;j_sum++)
       sum += y[j_sum]*cheby_poly(x_chebgrid[j_sum],j);    

     sum += y[N]*cheby_poly(x_chebgrid[N],j) * 0.5;
     sum *= 2.0/N;
     a.push_back(sum);
   }
   return a;
}

void cheby_expDn_matrix(vector_real &D, 
			vector_complex &Dc,
			vector_real &x, 
			unsigned long int N, 
			double dt)
{

  complex Coeff = complex(0.0,1.0) * 0.5 * dt;

  vector_complex D2((N+1)*(N+1),0.0),D4((N+1)*(N+1),0.0),D6((N+1)*(N+1),0.0),D8((N+1)*(N+1),0.0); 
  
  D.resize((N+1)*(N+1));
  Dc.resize((N+1)*(N+1));
  vector_real dX((N+1)*(N+1));

  for(int i=0; i<(N+1)*(N+1);i++) {
    D[i] = 0.0;
    Dc[i] = complex(0.0,0.0);
  }

  for(int i=0; i<N+1;i++){
    for(int j=0; j<N+1;j++){
      dX[i*(N+1)+j] = x[j] - x[i];
      if (i==j) dX[i*(N+1)+j] += 1.0;
      if      ((i==0 || i==N) && (j>0  && j<N))   D[i*(N+1)+j] = pow(-1,j+i) * 2.0; 
      else if ((i==0 || i==N) && (j==0 || j==N)) D[i*(N+1)+j] = pow(-1,j+i) * 1.0; 
      else if ((j==0 || j==N) && (i>0  && i<N)) D[i*(N+1)+j] = pow(-1,j+i) * 0.5;
      else D[i*(N+1)+j] = pow(-1,j+i) * 1.0;
      D[i*(N+1)+j] /= dX[i*(N+1)+j]; 
      //      std::cout<< D[i*(N+1)+j] << " "; 
    }
    //    std::cout <<"\n";
  }

  // correction diagonal
  double sum[N+1];
  for(int i=0; i<N+1;i++){
    sum[i] = 0;
    for(int j=0; j<N+1;j++) sum[i] += D[i*(N+1)+j];        
  }
  for(int i=0; i<N+1;i++) D[i*(N+1)+i] -= sum[i];

  // make D^2

  for(int i=0; i<(N+1)*(N+1);i++) Dc[i] = complex(0.0,0.0);

  for(int i=0; i<N+1;i++) {
    for(int j=0; j<N+1; j++){
      dX[i*(N+1)+j] = 0.0;
      for(int k=0; k<N+1; k++){
	dX[i*(N+1)+j] += D[i*(N+1)+k]* D[k*(N+1)+j] ;
      }
    }
  }

  for(int i=0; i<(N+1)*(N+1);i++) {
    //    D[i] = dX[i];
    Dc[i] += dX[i] * Coeff;
    D2[i] = dX[i];
    //    std::cout << dX[i] << "\n";
  }

  // // make D^4
  for(int i=0; i<N+1;i++) {
    for(int j=0; j<N+1; j++){
      dX[i*(N+1)+j] = 0.0;
      for(int k=0; k<N+1; k++){
  	dX[i*(N+1)+j] += D2[i*(N+1)+k].real()* D2[k*(N+1)+j].real() ;
      }
    }
  }
  for(int i=0; i<(N+1)*(N+1);i++) {
    // D[i] = dX[i];
    Dc[i] += dX[i] * Coeff*Coeff * 0.5;
    D4[i] = dX[i];
    //    std::cout << dX[i] << "\n";
  }

  // // make D^6
  for(int i=0; i<N+1;i++) {
    for(int j=0; j<N+1; j++){
      dX[i*(N+1)+j] = 0.0;
      for(int k=0; k<N+1; k++){
  	dX[i*(N+1)+j] += D4[i*(N+1)+k].real()* D2[k*(N+1)+j].real() ;
      }
    }
  }
  for(int i=0; i<(N+1)*(N+1);i++) {
    Dc[i] += dX[i] * Coeff*Coeff*Coeff /6.0;
    D6[i] = dX[i];
  }

  // // make D^8
  for(int i=0; i<N+1;i++) {
    for(int j=0; j<N+1; j++){
      dX[i*(N+1)+j] = 0.0;
      for(int k=0; k<N+1; k++){
  	dX[i*(N+1)+j] += D6[i*(N+1)+k].real()* D2[k*(N+1)+j].real() ;
      }
    }
  }

  for(int i=0; i<(N+1)*(N+1);i++) {
    Dc[i] += dX[i] * Coeff*Coeff*Coeff*Coeff / 24.0;
    // D8[i] = dX[i] * Coeff*Coeff*Coeff*Coeff / 24.0;
    // D6[i] *= Coeff*Coeff*Coeff / 6.0;
    // D4[i] *= 0.5 * Coeff*Coeff;
    // D2[i] *= Coeff;
  }

}

vector_complex cheby_interp(vector_real &x, vector_complex &coeff,double N,double a,double b)
{
  int Nx = x.size();
  vector_complex interp(Nx);
  double x_chebgrid;

  for (int ix=0; ix<Nx;ix++){
    x_chebgrid = ( 2.0 * x[ix] - a  - b ) / ( b-a);

    interp[ix]  = coeff[0]*cheby_poly(x_chebgrid,0) * 0.5;

    for (int j=1;j<N;j++)
      interp[ix] += coeff[j]*cheby_poly(x_chebgrid,j);

    interp[ix] += coeff[N]*cheby_poly(x_chebgrid,N) * 0.5;
  }

  return interp ;	  
}
complex cheby_interp(double &x, double &y, vector_complex &coeff,double N,
		     double ax,double bx, double ay, double by)
{
  //  int Nx = x.size();
  complex interp;
  double x_chebgrid,y_chebgrid;
  double Tx0, Ty0, Tx1, Ty1;
  vector_real polyx(N+1),polyy(N+1);

  x_chebgrid = -( 2.0 * x - ax  - bx ) / ( bx-ax);
  y_chebgrid = -( 2.0 * y - ay  - by ) / ( by-ay);

  for (int i=0;i<N+1;i++){
    polyx[i] = cheby_poly(x_chebgrid,i);
    polyy[i] = cheby_poly(y_chebgrid,i);
  }
    

  interp = complex(0.0,0.0) ;
  for (int i=0;i<N+1;i++)
  for (int j=0;j<N+1;j++)
      interp += coeff[i*(N+1)+j]*polyy[j]*polyx[i];
//      interp += coeff[i*(N+1)+j]*cheby_poly(y_chebgrid,j)*cheby_poly(x_chebgrid,i);
//      interp += coeff[i*(N+1)+j]*polyy[j]*polyx[i];

  return interp ;	  
}

void cheby_DD_matrix(vector_real &D, vector_real &x, unsigned long int N)
{

  D.resize((N+1)*(N+1));
  vector_real dX((N+1)*(N+1));
  for(int i=0; i<N+1;i++){
    for(int j=0; j<N+1;j++){
      dX[i*(N+1)+j] = x[j] - x[i];
      if (i==j) dX[i*(N+1)+j] += 1.0;
      if      ((i==0 || i==N) && (j>0  && j<N))   D[i*(N+1)+j] = pow(-1,j+i) * 2.0; 
      else if ((i==0 || i==N) && (j==0 || j==N)) D[i*(N+1)+j] = pow(-1,j+i) * 1.0; 
      else if ((j==0 || j==N) && (i>0  && i<N)) D[i*(N+1)+j] = pow(-1,j+i) * 0.5;
      else D[i*(N+1)+j] = pow(-1,j+i) * 1.0;
      D[i*(N+1)+j] /= dX[i*(N+1)+j]; 
      //      std::cout<< D[i*(N+1)+j] << " "; 
    }
    //    std::cout <<"\n";
  }

  // correction diagonal
  double sum[N+1];
  for(int i=0; i<N+1;i++){
    sum[i] = 0;
    for(int j=0; j<N+1;j++) sum[i] += D[i*(N+1)+j];        
  }
  for(int i=0; i<N+1;i++) D[i*(N+1)+i] -= sum[i];

  // make D^2
  for(int i=0; i<N+1;i++) {
    for(int j=0; j<N+1; j++){
      dX[i*(N+1)+j] = 0.0;
      for(int k=0; k<N+1; k++){
	dX[i*(N+1)+j] += D[i*(N+1)+k]* D[k*(N+1)+j] ;
      }
    }
  }
  for(int i=0; i<(N+1)*(N+1);i++) D[i] = dX[i];
}
complex GWP_vals(double x,double y, double Time)
{
  
  double hbar = 1.;
  double m =1.;
  double  sigma_px = 8, sigma_py=8;
  //  double  vx0 =  v0_init ; //2e1;
  double  vx0 =  -5e0,  vy0 = -5e0; //2e1;
  double  x0 = 0.5;
  double  y0 = 0.5;

  double sigma_x, phase_x, M_x;
  double sigma_y, phase_y, M_y;

  sigma_x = hbar/(2*sigma_px) * pow(1+4*pow(sigma_px,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
  phase_x = 1./hbar * (m*vx0+sigma_px*sigma_px/(sigma_x*sigma_x)*Time/(2*m)*(x-x0-vx0*Time))*
    (x-x0-vx0*Time) + vx0*m/(2*hbar) * vx0*Time - atan(2*sigma_px*sigma_px*Time/(hbar*m))/ 2 ;
  M_x = 1./(pow(2*M_PI,1./4)*pow(sigma_x,1./2)) * exp(-(x-x0-vx0*Time)*(x-x0-vx0*Time)/(4*sigma_x*sigma_x)) ;

  sigma_y = hbar/(2*sigma_py) * pow(1+4*pow(sigma_py,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
  phase_y = 1./hbar * (m*vy0+sigma_py*sigma_py/(sigma_y*sigma_y)*Time/(2*m)*(y-y0-vy0*Time))*
    (y-y0-vy0*Time) + vy0*m/(2*hbar) * vy0*Time - atan(2*sigma_py*sigma_py*Time/(hbar*m))/ 2 ;
  M_y = 1./(pow(2*M_PI,1./4)*pow(sigma_y,1./2)) * exp(-(y-y0-vy0*Time)*(y-y0-vy0*Time)/(4*sigma_y*sigma_y));

  return complex(M_x*M_y*cos(phase_x+phase_y),M_x*M_y*sin(phase_x+phase_y));
}

void GWP(Params &par, Operators &opr,double Time)
{
  
  unsigned long int N = opr.size;
  double hbar = 1.;
  double m =1.;
  double  sigma_px = 8;
  double  vx0 =  v0_init ; //2e1;
  double  x0 = x0_init;

  //  double Time = 0;
  double sigma_x, phase_x, M_x;


  for (size_t i = 0; i < N+1; ++i){
    sigma_x = hbar/(2*sigma_px) * pow(1+4*pow(sigma_px,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
    phase_x = 1./hbar * (m*vx0+sigma_px*sigma_px/(sigma_x*sigma_x)*Time/(2*m)*(par.x[i]-x0-vx0*Time))*
      (par.x[i]-x0-vx0*Time) + vx0*m/(2*hbar) * vx0*Time - atan(2*sigma_px*sigma_px*Time/(hbar*m))/ 2 ;
    M_x = 1./(pow(2*M_PI,1./4)*pow(sigma_x,1./2)) * exp(-(par.x[i]-x0-vx0*Time)*(par.x[i]-x0-vx0*Time)/(4*sigma_x*sigma_x)) ;

    opr.phi[i] = M_x*M_x; //imag
	 //	 std::cout << opr.wfc[k+N*(j+N*i)] <<  "\n"; 
    if (i==0) {
      // sigma_x = hbar/(2*sigma_px) * pow(1+4*pow(sigma_px,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
      // phase_x = 1./hbar * (m*vx0+sigma_px*sigma_px/(sigma_x*sigma_x)*Time/(2*m)*(-x0-vx0*Time))*
      // 	(-x0-vx0*Time) + vx0*m/(2*hbar) * vx0*Time - atan(2*sigma_px*sigma_px*Time/(hbar*m))/ 2 ;
      // M_x = 1./(pow(2*M_PI,1./4)*pow(sigma_x,1./2)) * exp(-(-x0-vx0*Time)*(-x0-vx0*Time)/(4*sigma_x*sigma_x)) ;

      sigma_x = hbar/(2*sigma_px) * pow(1+4*pow(sigma_px,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
      phase_x = 1./hbar * (m*vx0+sigma_px*sigma_px/(sigma_x*sigma_x)*Time/(2*m)*(par.x[i]-x0-vx0*Time))*
      	(par.x[i]-x0-vx0*Time) + vx0*m/(2*hbar) * vx0*Time - atan(2*sigma_px*sigma_px*Time/(hbar*m))/ 2 ;
      M_x = 1./(pow(2*M_PI,1./4)*pow(sigma_x,1./2)) * exp(-(par.x[i]-x0-vx0*Time)*(par.x[i]-x0-vx0*Time)/(4*sigma_x*sigma_x)) ;
      //opr.wfc[i].real( M_x*cos(phase_x)); 
      //opr.wfc[i].imag( M_x*sin(phase_x)); 
      par.a_R = M_x*cos(phase_x); 
      par.a_I = M_x*sin(phase_x); 
    }
    if (i==N) {
      // sigma_x = hbar/(2*sigma_px) * pow(1+4*pow(sigma_px,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
      // phase_x = 1./hbar * (m*vx0+sigma_px*sigma_px/(sigma_x*sigma_x)*Time/(2*m)*(par.xmax-x0-vx0*Time))*
      // 	(par.xmax-x0-vx0*Time) + vx0*m/(2*hbar) * vx0*Time - atan(2*sigma_px*sigma_px*Time/(hbar*m))/ 2 ;
      // M_x = 1./(pow(2*M_PI,1./4)*pow(sigma_x,1./2)) * exp(-(par.xmax-x0-vx0*Time)*(par.xmax-x0-vx0*Time)/(4*sigma_x*sigma_x)) ;

      sigma_x = hbar/(2*sigma_px) * pow(1+4*pow(sigma_px,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
      phase_x = 1./hbar * (m*vx0+sigma_px*sigma_px/(sigma_x*sigma_x)*Time/(2*m)*(par.x[i]-x0-vx0*Time))*
      	(par.x[i]-x0-vx0*Time) + vx0*m/(2*hbar) * vx0*Time - atan(2*sigma_px*sigma_px*Time/(hbar*m))/ 2 ;
      M_x = 1./(pow(2*M_PI,1./4)*pow(sigma_x,1./2)) * exp(-(par.x[i]-x0-vx0*Time)*(par.x[i]-x0-vx0*Time)/(4*sigma_x*sigma_x)) ;
      // opr.wfc[i].real( M_x*cos(phase_x)); 
      //opr.wfc[i].imag( M_x*sin(phase_x)); 
      par.b_R = M_x*cos(phase_x); 
      par.b_I = M_x*sin(phase_x); 
    }

  }
}
void initial_GWP(vector_real &x, 
		 vector_real &y,
		 vector_complex &psi)
{
  
  unsigned long int N = x.size();
  double Time = 0;

  for (size_t i = 0; i < N; ++i){
  for (size_t j = 0; j < N; ++j){
    psi[i*N+j] = GWP_vals(x[i], y[j], Time); //imag
    //    std::cout << x[i]<< " " << std::abs(psi[i]) <<  "\n"; 
  }}
  //  write_data(par, opr, 0);

}
void Print(Params &par,
	   vector_real &x,
	   vector_real &y,
	   vector_complex  &psi,
	   int tn){

  double size = x.size()-1;

  for (int i=0; i< size+1; i++){
  for (int j=0; j< size+1; j++){
    std::cout << x[i] << " " 
  	      << std::abs(psi[i]) << " " 
  	      << std::abs(GWP_vals(x[i],y[j], tn*par.dt))<< " "
  	      << psi[i*(size+1)+j].real() << " " 
  	      << GWP_vals(x[i],y[j], tn*par.dt).real()<< " " 
  	      << psi[i*(size+1)+j].imag() << " " 
  	      << GWP_vals(x[i],y[j], tn*par.dt).imag()<< "\n"; 
  }}
}
void Print_file(Params &par,
		vector_patch &P,
		unsigned long int NsubD,
		int tn){

  int size,D;
  int Ng = 8;

  std::stringstream filename_fdm;
  filename_fdm << "output/" << tn << ".dat";
  std::ofstream ffdm = std::ofstream(filename_fdm.str());

  if (ffdm)
  {
    std::stringstream data_fdm;

    for (int Di=0;Di<NsubD;Di++){
    for (int i = Ng; i < P[0].N-Ng; ++i){
      for (int Dj=0;Dj<NsubD;Dj++){
	D = Di*NsubD+Dj;

	for (int j = Ng; j < P[D].N-Ng; ++j){

	  data_fdm << std::setprecision(20)
		   << i + Di*(P[D].N+1) << " "
		   << j + Dj*(P[D].N+1) << " " 
  		   << P[D].x[i] << " " 
  		   << P[D].y[j] << " " 
  		   << std::abs(P[D].psi[i*(P[D].N+1)+j]) << " " 
  		   << std::abs(GWP_vals(P[D].x[i],P[D].y[j], tn*par.dt))<< " "
  		   << P[D].psi[i*(P[D].N+1)+j].real() << " " 
  		   << GWP_vals(P[D].x[i],P[D].y[j], tn*par.dt).real()<< " " 
  		   << P[D].psi[i*(P[D].N+1)+j].imag() << " " 
  		   << GWP_vals(P[D].x[i], P[D].y[j],tn*par.dt).imag()<< " "
  		   << std::abs(P[D].debug[i*(P[D].N+1)+j]) << "\n" ;
  	}
      }
      }
    }
    ffdm.write(data_fdm.str().c_str(), data_fdm.str().length());
    ffdm.close();
  }
  //  }
}
void Print_Patch(Params &par,
		 vector_patch &Patch,
		 unsigned long int NsubD,
		 int tn,
		 int D){

  int size;
  int Ng = 8;

  std::stringstream filename_fdm;
  filename_fdm << "output/P" << D<<"_"<< tn << ".dat";
  std::ofstream ffdm = std::ofstream(filename_fdm.str());

  if (ffdm)
  {
    std::stringstream data_fdm;

     for (int t=0; t<Patch[D].N+1; t++)
	for (int tt=0; tt<Patch[D].N+1; tt++)
	  data_fdm << std::setprecision(20)
		   << Patch[D].x[t] << " "
		    << Patch[D].y[tt] << " "
		    << t << " "
		    << tt << " "
		    << Patch[D].psi[t*(Patch[D].N+1)+tt].real() << " "
		    << cheby_interp(Patch[D].x[t],
				    Patch[D].y[tt],
				    Patch[D].coeff,
				    Patch[D].N,
				    Patch[D].x[0],
				    Patch[D].x.back(),
				    Patch[D].y[0],
				    Patch[D].y.back()).real() << " "
		    << Patch[D].coeff[t*(Patch[D].N+1)+tt].real() << " "
  		   << GWP_vals(Patch[D].x[t],Patch[D].y[tt], tn*par.dt).real()<< "\n"; 

    ffdm.write(data_fdm.str().c_str(), data_fdm.str().length());
    ffdm.close();
  }


  //  }
}
void Print_InterpolationPatch(Params &par,
		 vector_patch &Patch,
		 unsigned long int NsubD,
		 int tn,
		 int D){

  int size;
  int Ng = 8;

  std::stringstream filename_fdm;
  filename_fdm << "output/IP" << D<<"_"<< tn << ".dat";
  std::ofstream ffdm = std::ofstream(filename_fdm.str());

  int N = 1e2;
  double xmin = Patch[D].x[0];
  double xmax = Patch[D].x.back();
  double ymin = Patch[D].y[0];
  double ymax = Patch[D].y.back();
  double x,y;
  if (ffdm)
  {
    std::stringstream data_fdm;

     for (int t=0; t<N; t++){
	for (int tt=0; tt<N; tt++){
	  x = xmin + t*(xmax-xmin)/(N-1);
	  y = ymin + tt*(ymax-ymin)/(N-1);
	  data_fdm << std::setprecision(25)
		   << x << " "
		   << y << " "
		   << t << " "
		   << tt << " "
		   << cheby_interp(x,
				   y,
				   Patch[D].coeff,
				   Patch[D].N,
				   Patch[D].x[0],
				   Patch[D].x.back(),
				   Patch[D].y[0],
				   Patch[D].y.back()).real() << " "
   		   << GWP_vals(x,y, tn*par.dt).real()<< " " 
		   << cheby_interp(x,
				   y,
				   Patch[D].coeff,
				   Patch[D].N,
				   Patch[D].x[0],
				   Patch[D].x.back(),
				   Patch[D].y[0],
				   Patch[D].y.back()).imag() << " "
   		   << GWP_vals(x,y, tn*par.dt).imag()<< "\n"; 

	}}
    ffdm.write(data_fdm.str().c_str(), data_fdm.str().length());
    ffdm.close();
  }

  //  }
}

void Error(Params &par,
	   vector_real &x,
	   vector_real &y,
	   vector_complex  &psi,
	   double &Error){
  
  int Ng = 8;
  double size = x.size();
  for (size_t i = Ng; i < size-Ng; ++i){
  for (size_t j = Ng; j < size-Ng; ++j){
    Error += pow(std::abs(psi[i*(size)+j])
	       - std::abs(GWP_vals(x[i],y[j], par.timesteps*par.dt)),2);
  }}
  //  error = pow(error,1./2)/opr.size;  
  //  std::cout << "error = " << error << "\n";    
}
void RK2(int N, vector_real &D,vector_complex &psi, Params &par){

  vector_complex k1(N+1),tmp(N+1);
    for(int i=1; i<N;i++) {
      k1[i] = complex(0.0,0.0);
            for(int j=1; j<N; j++) k1[i] += D[i*(N+1)+j] * psi[j] ;
    }
    for(int i=1; i<N;i++) k1[i] = psi[i] + k1[i] * 0.25 * par.dt * complex(0.0,1.0);

    for(int i=1; i<N;i++) {
      tmp[i] = complex(0.0,0.0);
      for(int j=1; j<N; j++) tmp[i] += D[i*(N+1)+j] * k1[j] ;
    }
    for(int i=1; i<N;i++) psi[i] += tmp[i] * 0.5 * par.dt * complex(0.0,1.0);
}

void Refinement(Patches &P){
  
  // Perform Refinement: double cheby grid size

  // for (int i=0; i< P.N+1; i++){
  //   std::cout << P.x[i] << " " 
  // 	      << P.psi[i].real() << " " 
  // 	      << P.psi[i].imag() <<"\n";
  // }
  P.coeff = cheby_coeff(P.x,P.psi,P.N);
  P.x  = cheby_exgrid(2*P.N,P.x[0],P.x[P.N]);
//  P.psi = cheby_interp(P.x,P.coeff,P.N,P.x[0],P.x[2*P.N]); 
  cheby_DD_matrix(P.D, P.x, 2*P.N);
  P.level = 1;

  P.N *= 2.0; 
  P.tmp.resize(P.N);
  P.k1.resize(P.N);
  P.k2.resize(P.N);
  P.k3.resize(P.N);
  P.k4.resize(P.N);
  // std::cout <<"\n";

  // for (int i=0; i< P.N+1; i++){
  //   std::cout << P.x[i] << " " 
  // 	      << P.psi[i].real() << " " 
  // 	      << P.psi[i].imag() <<"\n";
  // }

  //   std::exit(0);
}

void Derefinement(Patches &P){
  
  // Perform De-Refinement: half cheby grid size

  // for (int i=0; i< P.N+1; i++){
  //   std::cout << P.x[i] << " " 
  // 	      << P.psi[i].real() << " " 
  // 	      << P.psi[i].imag() <<"\n";
  // }
  double half = 0.5 ; 
  P.coeff = cheby_coeff(P.x,P.psi,P.N);
  P.x  = cheby_exgrid(half*P.N,P.x[0],P.x[P.N]);
//  P.psi = cheby_interp(P.x,P.coeff,P.N,P.x[0],P.x[half*P.N]); 
  cheby_DD_matrix(P.D, P.x, half*P.N);
  P.level = 0;

  P.N *= half; 
  P.tmp.resize(P.N);
  P.k1.resize(P.N);
  P.k2.resize(P.N);
  P.k3.resize(P.N);
  P.k4.resize(P.N);
  //  std::cout <<"\n";

  // for (int i=0; i< P.N+1; i++){
  //   std::cout << P.x[i] << " " 
  // 	      << P.psi[i].real() << " " 
  // 	      << P.psi[i].imag() <<"\n";
  // }

  //   std::exit(0);
}

void Assign_BC(vector_patch &Patch,
	       signed long int NsubD,
	       Params &par,
	       double tn){


  //    std::cout << tn <<"\n";
  //  unsigned long int N = Patch[0].N;
  vector_complex psi1D(Patch[0].N+1);

  // get Cheby Coeffs
  for (int D=0; D<NsubD*NsubD; D++)
    cheby_coeff(Patch[D].psi,Patch[D].N,Patch[D].coeff,Patch[D]);

  for (int D=0; D<NsubD*NsubD; D++){
    for (int t=0; t<Patch[D].N_ghost; t++){
      for (int tt=0; tt<Patch[D].N+1; tt++){
	Patch[D].psi[t*(Patch[D].N+1)+tt] = complex(-1e10,-1e10);
	Patch[D].psi[tt*(Patch[D].N+1)+t] = complex(-1e10,-1e10);
	Patch[D].psi[(Patch[D].N-t)*(Patch[D].N+1)+tt] = complex(-1e10,-1e10);
	Patch[D].psi[tt*(Patch[D].N+1)+(Patch[D].N-t)] = complex(-1e10,-1e10);
	  }}}
  

  // compute interpolation at ghost zones
  int ghost_pts = Patch[0].N_ghost;
  int N , D, Dfrom;

  // int DD = 2;
  // int DD = 0;
  // if (tn == 1){
  //    for (int t=0; t<Patch[DD].N+1; t++)
  // 	for (int tt=0; tt<Patch[DD].N+1; tt++)
  // 	  std::cout << Patch[DD].x[t] << " "
  // 		    << Patch[DD].y[tt] << " "
  // 		    << t << " "
  // 		    << tt << " "
  // 		    << Patch[DD].psi[t*(Patch[DD].N+1)+tt].real() << " "
  // 		    << cheby_interp(Patch[DD].x[t],
  // 				    Patch[DD].y[tt],
  // 				    Patch[DD].coeff,
  // 				    Patch[DD].N,
  // 				    Patch[DD].x[0],
  // 				    Patch[DD].x.back(),
  // 				    Patch[DD].y[0],
  // 				    Patch[DD].y.back()).real() << " "
  // 		    << Patch[DD].coeff[t*(Patch[DD].N+1)+tt].real() << "\n";
  // 		  // << Patch[Dfrom].tmp[t*(N+1)+j].real() << " "
  // 		  // << Patch[D].psi[t*(N+1)+j].imag() << " "
  // 		  // << Patch[Dfrom].tmp[t*(N+1)+j].imag() << "\n";
  //     // 	std::cout << Patch[D].x[t] << " "
  //     // 		  << Patch[Dfrom].x[t] << " "
  //     // 		  << Patch[D].psi[t*(N+1)+j].real() << " "
  //     // 		  << Patch[Dfrom].tmp[t*(N+1)+j].real() << " "
  //     // 		  << Patch[D].psi[t*(N+1)+j].imag() << " "
  //     // 		  << Patch[Dfrom].tmp[t*(N+1)+j].imag() << " "
  //     // 		  << Patch[D].y[j] << "\n";
  // }
  // std::exit(0);}

  for (int Di=0;Di<NsubD;Di++) {
  for (int Dj=0;Dj<NsubD;Dj++) {
    D = Di*NsubD + Dj;
    N = Patch[D].N;

    // left face
    if (Di==0){      
      for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
      for (int ig=0; ig< ghost_pts; ig++){
	Patch[D].psi[ig*(N+1)+jg] =  GWP_vals(Patch[D].x[ig],Patch[D].y[jg],(tn-1.)*par.dt);
	Patch[D].tmp[ig*(N+1)+jg] =  GWP_vals(Patch[D].x[ig],Patch[D].y[jg],(tn-1.)*par.dt);
      }}
    }
    else if (Patch[(Di-1)*NsubD+Dj].level==0){
      Dfrom = (Di-1)*NsubD+Dj;

      for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
      for (int ig=0; ig< ghost_pts; ig++){
	  Patch[D].psi[ig*(N+1)+jg] =  cheby_interp(Patch[D].x[ig],
						    Patch[D].y[jg],
						    Patch[Dfrom].coeff,
						    Patch[Dfrom].N,
						    Patch[Dfrom].x[0],
						    Patch[Dfrom].x.back(),
						    Patch[Dfrom].y[0],
						    Patch[Dfrom].y.back());  
	  Patch[D].tmp[ig*(N+1)+jg] = Patch[D].psi[ig*(N+1)+jg];
      }}
    }
//     // right face
    if (Di==(NsubD-1)){
      for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
      for (int ig=0; ig< ghost_pts; ig++){
	Patch[D].psi[(N-ig)*(N+1)+jg] =  GWP_vals(Patch[D].x[N-ig],Patch[D].y[jg],(tn-1.)*par.dt);
	Patch[D].tmp[(N-ig)*(N+1)+jg] =  GWP_vals(Patch[D].x[N-ig],Patch[D].y[jg],(tn-1.)*par.dt);
	//	std::cout << "b"<<j<<"="<< Patch[D].b[j] <<"\n";
      }}
    }
    else if (Patch[(Di+1)*NsubD+Dj].level==0){
      Dfrom = (Di+1)*NsubD+Dj;

      for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
      for (int ig=0; ig< ghost_pts; ig++){

	  Patch[D].psi[(N-ig)*(N+1)+jg] =  cheby_interp(Patch[D].x[N-ig],
							Patch[D].y[jg],
							Patch[Dfrom].coeff,
							Patch[Dfrom].N,
							Patch[Dfrom].x[0],
							Patch[Dfrom].x.back(),
							Patch[Dfrom].y[0],
							Patch[Dfrom].y.back());  
	  Patch[D].tmp[(N-ig)*(N+1)+jg] = Patch[D].psi[(N-ig)*(N+1)+jg];
      } }
    }
//     // C face
    if (Dj==0){
      for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      for (int jg=0; jg< ghost_pts; jg++){
	Patch[D].psi[ig*(N+1)+jg] =  GWP_vals(Patch[D].x[ig],Patch[D].y[jg],(tn-1.)*par.dt);
	Patch[D].tmp[ig*(N+1)+jg] =  GWP_vals(Patch[D].x[ig],Patch[D].y[jg],(tn-1.)*par.dt);
	//	std::cout << "c"<<i<<"="<< Patch[D].c[i] <<"\n";
      }}
    }
    else if (Patch[Di*NsubD+Dj-1].level==0){
      Dfrom = Di*NsubD+Dj-1;
      for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      for (int jg=0; jg< ghost_pts; jg++){

	  Patch[D].psi[ig*(N+1)+jg] =  cheby_interp(Patch[D].x[ig],
						    Patch[D].y[jg],
						    Patch[Dfrom].coeff,
						    Patch[Dfrom].N,
						    Patch[Dfrom].x[0],
						    Patch[Dfrom].x.back(),
						    Patch[Dfrom].y[0],
						    Patch[Dfrom].y.back());  
	  Patch[D].tmp[ig*(N+1)+jg] = Patch[D].psi[ig*(N+1)+jg];
      }
      }
    }

//     // D face
    if (Dj==(NsubD-1)){
      for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      for (int jg=0; jg< ghost_pts; jg++){
	Patch[D].psi[ig*(N+1)+N-jg] =  GWP_vals(Patch[D].x[ig],Patch[D].y[N-jg],(tn-1.)*par.dt);
	Patch[D].tmp[ig*(N+1)+N-jg] =  GWP_vals(Patch[D].x[ig],Patch[D].y[N-jg],(tn-1.)*par.dt);
	//	std::cout << "d"<<i<<"="<< Patch[D].d[i] <<"\n";
      }}
    }
    else if (Patch[Di*NsubD+Dj+1].level==0){
      Dfrom = Di*NsubD+Dj+1;
      for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      for (int jg=0; jg< ghost_pts; jg++){

	  Patch[D].psi[ig*(N+1)+N-jg] =  cheby_interp(Patch[D].x[ig],
						    Patch[D].y[N-jg],
						    Patch[Dfrom].coeff,
						    Patch[Dfrom].N,
						    Patch[Dfrom].x[0],
						    Patch[Dfrom].x.back(),
						    Patch[Dfrom].y[0],
						    Patch[Dfrom].y.back());  
	  Patch[D].tmp[ig*(N+1)+N-jg] = Patch[D].psi[ig*(N+1)+N-jg];
      }}
    }

    // lower-left corner
   if (Di==0 || (Dj==0)){
	for (int jg=0; jg< ghost_pts; jg++){
	for (int ig=0; ig< ghost_pts; ig++){
	  Patch[D].psi[ig*(N+1)+jg] =  GWP_vals(Patch[D].x[ig],Patch[D].y[jg],(tn-1.)*par.dt);
	  Patch[D].tmp[ig*(N+1)+jg] =  GWP_vals(Patch[D].x[ig],Patch[D].y[jg],(tn-1.)*par.dt);
      }}
   }
   else if (Patch[(Di-1)*NsubD+Dj-1].level==0){
	Dfrom = (Di-1)*NsubD+Dj-1;
	for (int jg=0; jg< ghost_pts; jg++){
	for (int ig=0; ig< ghost_pts; ig++){
	    
	  Patch[D].psi[ig*(N+1)+jg] =  cheby_interp(Patch[D].x[ig],
						    Patch[D].y[jg],
						    Patch[Dfrom].coeff,
						    Patch[Dfrom].N,
						    Patch[Dfrom].x[0],
						    Patch[Dfrom].x.back(),
						    Patch[Dfrom].y[0],
						    Patch[Dfrom].y.back());
	  Patch[D].tmp[ig*(N+1)+jg] = Patch[D].psi[ig*(N+1)+jg];
	}}
   }

//     // upper-left corner
   if (Di==0 || (Dj==(NsubD-1))){
	for (int ig=0; ig< ghost_pts; ig++){
	for (int jg=0; jg< ghost_pts; jg++){	    
	  Patch[D].psi[ig*(N+1)+N-jg] =  GWP_vals(Patch[D].x[ig],Patch[D].y[N-jg],(tn-1.)*par.dt);
	  Patch[D].tmp[ig*(N+1)+N-jg] =  GWP_vals(Patch[D].x[ig],Patch[D].y[N-jg],(tn-1.)*par.dt);
      }}
   }
   else if (Patch[(Di-1)*NsubD+Dj+1].level==0){
	Dfrom = (Di-1)*NsubD+Dj+1;
	for (int ig=0; ig< ghost_pts; ig++){
	for (int jg=0; jg< ghost_pts; jg++){	    
	  Patch[D].psi[ig*(N+1)+N-jg] =  cheby_interp(Patch[D].x[ig],
						      Patch[D].y[N-jg],
						      Patch[Dfrom].coeff,
						      Patch[Dfrom].N,
						      Patch[Dfrom].x[0],
						      Patch[Dfrom].x.back(),
						      Patch[Dfrom].y[0],
						      Patch[Dfrom].y.back());  
	  Patch[D].tmp[ig*(N+1)+N-jg] = Patch[D].psi[ig*(N+1)+N-jg];
      }}
   }   

//     // lower-right corner
   if (Di==(NsubD-1) || (Dj==0)){
	for (int ig=0; ig< ghost_pts; ig++){
	for (int jg=0; jg< ghost_pts; jg++){	    
	  Patch[D].psi[(N-ig)*(N+1)+jg] =  GWP_vals(Patch[D].x[N-ig],Patch[D].y[jg],(tn-1.)*par.dt);
	  Patch[D].tmp[(N-ig)*(N+1)+jg] =  GWP_vals(Patch[D].x[N-ig],Patch[D].y[jg],(tn-1.)*par.dt);
      }}
   } 
  else if (Patch[(Di+1)*NsubD+Dj-1].level==0){
 	Dfrom = (Di+1)*NsubD+Dj-1;
	for (int ig=0; ig< ghost_pts; ig++){
	for (int jg=0; jg< ghost_pts; jg++){	    
	  Patch[D].psi[(N-ig)*(N+1)+jg] =  cheby_interp(Patch[D].x[(N-ig)],
							Patch[D].y[jg],
							Patch[Dfrom].coeff,
							Patch[Dfrom].N,
							Patch[Dfrom].x[0],
							Patch[Dfrom].x.back(),
							Patch[Dfrom].y[0],
							Patch[Dfrom].y.back());  
	  Patch[D].tmp[(N-ig)*(N+1)+jg] = Patch[D].psi[(N-ig)*(N+1)+jg];
      }}
  }
   

//     // upper-right corner
   if (Di==(NsubD-1) || (Dj==(NsubD-1))){
	for (int ig=0; ig< ghost_pts; ig++){
	for (int jg=0; jg< ghost_pts; jg++){	    
	  Patch[D].psi[(N-ig)*(N+1)+N-jg] =  GWP_vals(Patch[D].x[N-ig],Patch[D].y[N-jg],(tn-1.)*par.dt);
	  Patch[D].tmp[(N-ig)*(N+1)+N-jg] =  GWP_vals(Patch[D].x[N-ig],Patch[D].y[N-jg],(tn-1.)*par.dt);
      }}
   }
   else if (Patch[(Di+1)*NsubD+Dj+1].level==0){
 	Dfrom = (Di+1)*NsubD+Dj+1;
	for (int ig=0; ig< ghost_pts; ig++){
	for (int jg=0; jg< ghost_pts; jg++){	    
	  Patch[D].psi[(N-ig)*(N+1)+N-jg] =  cheby_interp(Patch[D].x[N-ig],
							  Patch[D].y[N-jg],
							  Patch[Dfrom].coeff,
							  Patch[Dfrom].N,
							  Patch[Dfrom].x[0],
							  Patch[Dfrom].x.back(),
							  Patch[Dfrom].y[0],
							  Patch[Dfrom].y.back());  
	  Patch[D].tmp[(N-ig)*(N+1)+N-jg] = Patch[D].psi[(N-ig)*(N+1)+N-jg];
      }}
   }
   }} //   for (int Di=0;Di<NsubD;Di++) { for (int Dj=0;Dj<NsubD;Dj++) {

//   int DD = 0;
//   if (tn == 1){
//      for (int t=0; t<Patch[DD].N+1; t++)
// 	for (int tt=0; tt<Patch[DD].N+1; tt++)
// 	  std::cout << Patch[DD].x[t] << " "
// 		    << Patch[DD].y[tt] << " "
// 		    << t << " "
// 		    << tt << " "
// 		    << Patch[DD].psi[t*(Patch[DD].N+1)+tt].real() << " "
// 		    << cheby_interp(Patch[DD].x[t],
// 				    Patch[DD].y[tt],
// 				    Patch[DD].coeff,
// 				    Patch[DD].N,
// 				    Patch[DD].x[0],
// 				    Patch[DD].x.back(),
// 				    Patch[DD].y[0],
// 				    Patch[DD].y.back()).real() << " "
// 		    << Patch[DD].coeff[t*(Patch[DD].N+1)+tt].real() << "\n";
// 		  // << Patch[Dfrom].tmp[t*(N+1)+j].real() << " "
// 		  // << Patch[D].psi[t*(N+1)+j].imag() << " "
// 		  // << Patch[Dfrom].tmp[t*(N+1)+j].imag() << "\n";
//       // 	std::cout << Patch[D].x[t] << " "
//       // 		  << Patch[Dfrom].x[t] << " "
//       // 		  << Patch[D].psi[t*(N+1)+j].real() << " "
//       // 		  << Patch[Dfrom].tmp[t*(N+1)+j].real() << " "
//       // 		  << Patch[D].psi[t*(N+1)+j].imag() << " "
//       // 		  << Patch[Dfrom].tmp[t*(N+1)+j].imag() << " "
//       // 		  << Patch[D].y[j] << "\n";

// //  std::exit(0);}
//   }
}
void Assign_CornerBC(vector_patch &Patch,
	       signed long int NsubD,
	       Params &par,
	       double tn){

  //    std::cout << tn <<"\n";
  //  unsigned long int N = Patch[0].N;
  vector_complex psi1D(Patch[0].N+1);

  // get Cheby Coeffs
  for (int D=0; D<NsubD*NsubD; D++)
    cheby_coeff(Patch[D].psi,Patch[D].N,Patch[D].coeff,Patch[D]);
  
  // compute interpolation at ghost zones
  int ghost_pts = Patch[0].N_ghost;
  int N , D, Dfrom;

  int LF;
  int RF;
  int UF;
  int LOF;
  int LLC;
  int ULC;
  int LRC;
  int URC;

  for (int Di=0;Di<NsubD;Di++) {
  for (int Dj=0;Dj<NsubD;Dj++) {
    LF=0;
    RF=0;
    UF=0;
    LOF=0;
    LLC=0;
    ULC=0;
    LRC=0;
    URC=0;

    D = Di*NsubD + Dj;
    N = Patch[D].N;

    // left face
    if (Di!=0){
      LF=1; 
      Dfrom = (Di-1)*NsubD+Dj;
      for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
      for (int ig=0; ig< ghost_pts; ig++){

	  Patch[D].psi[ig*(N+1)+jg] =  cheby_interp(Patch[D].x[ig],
						    Patch[D].y[jg],
						    Patch[Dfrom].coeff,
						    Patch[Dfrom].N,
						    Patch[Dfrom].x[0],
						    Patch[Dfrom].x.back(),
						    Patch[Dfrom].y[0],
						    Patch[Dfrom].y.back());  
	  Patch[D].tmp[ig*(N+1)+jg] = Patch[D].psi[ig*(N+1)+jg];
      }
      }
    }

    // right face
    if (Di!=(NsubD-1)){
      Dfrom = (Di+1)*NsubD+Dj;
      RF=1; 
      for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
      for (int ig=0; ig< ghost_pts; ig++){

	  Patch[D].psi[(N-ig)*(N+1)+jg] =  cheby_interp(Patch[D].x[N-ig],
							Patch[D].y[jg],
							Patch[Dfrom].coeff,
							Patch[Dfrom].N,
							Patch[Dfrom].x[0],
							Patch[Dfrom].x.back(),
							Patch[Dfrom].y[0],
							Patch[Dfrom].y.back());  
	  Patch[D].tmp[(N-ig)*(N+1)+jg] = Patch[D].psi[(N-ig)*(N+1)+jg];
      }
      }
    }

    // C face
    if (Dj!=0){

      LOF=1; 
      Dfrom = Di*NsubD+Dj-1;
      for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      for (int jg=0; jg< ghost_pts; jg++){

	  Patch[D].psi[ig*(N+1)+jg] =  cheby_interp(Patch[D].x[ig],
						    Patch[D].y[jg],
						    Patch[Dfrom].coeff,
						    Patch[Dfrom].N,
						    Patch[Dfrom].x[0],
						    Patch[Dfrom].x.back(),
						    Patch[Dfrom].y[0],
						    Patch[Dfrom].y.back());  
	  Patch[D].tmp[ig*(N+1)+jg] = Patch[D].psi[ig*(N+1)+jg];
      }
      }
    }

    // D face
    if (Dj!=(NsubD-1)){
      UF=1; 
      Dfrom = Di*NsubD+Dj+1;
      for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      for (int jg=0; jg< ghost_pts; jg++){

	  Patch[D].psi[ig*(N+1)+N-jg] =  cheby_interp(Patch[D].x[ig],
						    Patch[D].y[N-jg],
						    Patch[Dfrom].coeff,
						    Patch[Dfrom].N,
						    Patch[Dfrom].x[0],
						    Patch[Dfrom].x.back(),
						    Patch[Dfrom].y[0],
						    Patch[Dfrom].y.back());  
	  Patch[D].tmp[ig*(N+1)+N-jg] = Patch[D].psi[ig*(N+1)+N-jg];
      }}
    }

    // lower-left corner
   if (Di!=0 && (Dj!=0)){
   if (Patch[(Di-1)*NsubD+Dj-1].level==0){
     LLC=1; 
     Dfrom = (Di-1)*NsubD+Dj-1;

     for (int jg=0; jg< ghost_pts; jg++){
       for (int ig=0; ig< ghost_pts; ig++){
	    
	  Patch[D].psi[ig*(N+1)+jg] =  cheby_interp(Patch[D].x[ig],
						    Patch[D].y[jg],
						    Patch[Dfrom].coeff,
						    Patch[Dfrom].N,
						    Patch[Dfrom].x[0],
						    Patch[Dfrom].x.back(),
						    Patch[Dfrom].y[0],
						    Patch[Dfrom].y.back());
	  Patch[D].tmp[ig*(N+1)+jg] = Patch[D].psi[ig*(N+1)+jg];
      }
	}
   }}

//     // upper-left corner
   if (Di!=0 && (Dj!=(NsubD-1))){
   if (Patch[(Di-1)*NsubD+Dj+1].level==0){
     ULC=1; 
     Dfrom = (Di-1)*NsubD+Dj+1;
	for (int ig=0; ig< ghost_pts; ig++){
	for (int jg=0; jg< ghost_pts; jg++){	    
	  Patch[D].psi[ig*(N+1)+N-jg] =  cheby_interp(Patch[D].x[ig],
						      Patch[D].y[N-jg],
						      Patch[Dfrom].coeff,
						      Patch[Dfrom].N,
						      Patch[Dfrom].x[0],
						      Patch[Dfrom].x.back(),
						      Patch[Dfrom].y[0],
						      Patch[Dfrom].y.back());  
	  Patch[D].tmp[ig*(N+1)+N-jg] = Patch[D].psi[ig*(N+1)+N-jg];
      }
	}
   }}   

//     // lower-right corner
   if (Di!=(NsubD-1) && (Dj!=0)){
   if (Patch[(Di+1)*NsubD+Dj-1].level==0){
     LRC=1; 
 	Dfrom = (Di+1)*NsubD+Dj-1;
	for (int ig=0; ig< ghost_pts; ig++){
	for (int jg=0; jg< ghost_pts; jg++){	    
	  Patch[D].psi[(N-ig)*(N+1)+jg] =  cheby_interp(Patch[D].x[(N-ig)],
							Patch[D].y[jg],
							Patch[Dfrom].coeff,
							Patch[Dfrom].N,
							Patch[Dfrom].x[0],
							Patch[Dfrom].x.back(),
							Patch[Dfrom].y[0],
							Patch[Dfrom].y.back());  
	  Patch[D].tmp[(N-ig)*(N+1)+jg] = Patch[D].psi[(N-ig)*(N+1)+jg];
      }}
   }}

//     // upper-right corner
   if (Di!=(NsubD-1) && (Dj!=(NsubD-1))){
   if (Patch[(Di+1)*NsubD+Dj+1].level==0){
     URC=1; 
 	Dfrom = (Di+1)*NsubD+Dj+1;
	for (int ig=0; ig< ghost_pts; ig++){
	for (int jg=0; jg< ghost_pts; jg++){	    
	  Patch[D].psi[(N-ig)*(N+1)+N-jg] =  cheby_interp(Patch[D].x[N-ig],
							  Patch[D].y[N-jg],
							  Patch[Dfrom].coeff,
							  Patch[Dfrom].N,
							  Patch[Dfrom].x[0],
							  Patch[Dfrom].x.back(),
							  Patch[Dfrom].y[0],
							  Patch[Dfrom].y.back());  
	  Patch[D].tmp[(N-ig)*(N+1)+N-jg] = Patch[D].psi[(N-ig)*(N+1)+N-jg];
      }}
   }}
   std::cout << Di << Dj
	     << " LF=" << LF
	     << " RF=" << RF
	     << " UF=" << UF
	     << " LOF=" << LOF
	     << " LLC=" << LLC
	     << " LRC=" << LRC
	     << " ULC=" << ULC
	     << " URC=" << URC << "\n";

  }} //   for (int Di=0;Di<NsubD;Di++) { for (int Dj=0;Dj<NsubD;Dj++) {
}
void Assign_AnalyticBC(vector_patch &Patch,
	       signed long int NsubD,
	       Params &par,
	       double tn){


  //    std::cout << tn <<"\n";
  //  unsigned long int N = Patch[0].N;
  vector_complex psi1D(Patch[0].N+1);

  // get Cheby Coeffs
  for (int D=0; D<NsubD*NsubD; D++)
    cheby_coeff(Patch[D].psi,Patch[D].N,Patch[D].coeff,Patch[D]);

  for (int D=0; D<NsubD*NsubD; D++){
    for (int t=0; t<Patch[D].N_ghost; t++){
      for (int tt=0; tt<Patch[D].N+1; tt++){
	Patch[D].psi[t*(Patch[D].N+1)+tt] = complex(-1e10,-1e10);
	Patch[D].psi[tt*(Patch[D].N+1)+t] = complex(-1e10,-1e10);
	Patch[D].psi[(Patch[D].N-t)*(Patch[D].N+1)+tt] = complex(-1e10,-1e10);
	Patch[D].psi[tt*(Patch[D].N+1)+(Patch[D].N-t)] = complex(-1e10,-1e10);
	  }}}  

  // compute interpolation at ghost zones
  int ghost_pts = Patch[0].N_ghost;
  int N , D, Dfrom;

  for (int Di=0;Di<NsubD;Di++) {
  for (int Dj=0;Dj<NsubD;Dj++) {
    D = Di*NsubD + Dj;
    N = Patch[D].N;

    // left face
    for (int ig=0; ig< ghost_pts; ig++){
      for (int j = 0; j<N+1; j++){
	Patch[D].psi[ig*(N+1)+j] =  GWP_vals(Patch[D].x[ig],Patch[D].y[j],(tn)*par.dt);
	Patch[D].tmp[ig*(N+1)+j] =  GWP_vals(Patch[D].x[ig],Patch[D].y[j],(tn)*par.dt);
      }}
    // right face
    for (int ig=0; ig< ghost_pts; ig++){
      for (int j = 0; j<N+1; j++){
	Patch[D].psi[(N-ig)*(N+1)+j] =  GWP_vals(Patch[D].x[N-ig],Patch[D].y[j],(tn)*par.dt);
	Patch[D].tmp[(N-ig)*(N+1)+j] =  GWP_vals(Patch[D].x[N-ig],Patch[D].y[j],(tn)*par.dt);
	//	std::cout << "b"<<j<<"="<< Patch[D].b[j] <<"\n";
      }}
    
    // C face
    for (int ig=0; ig< ghost_pts; ig++){
      for (int i = 0; i<N+1; i++){
	Patch[D].psi[i*(N+1)+ig] =  GWP_vals(Patch[D].x[i],Patch[D].y[ig],(tn)*par.dt);
	Patch[D].tmp[i*(N+1)+ig] =  GWP_vals(Patch[D].x[i],Patch[D].y[ig],(tn)*par.dt);
	//	std::cout << "c"<<i<<"="<< Patch[D].c[i] <<"\n";
      }}
    
    // D face
    for (int ig=0; ig< ghost_pts; ig++){
      for (int i =0; i<N+1; i++){
	Patch[D].psi[i*(N+1)+N-ig] =  GWP_vals(Patch[D].x[i],Patch[D].y[N-ig],(tn)*par.dt);
	Patch[D].tmp[i*(N+1)+N-ig] =  GWP_vals(Patch[D].x[i],Patch[D].y[N-ig],(tn)*par.dt);
	//	std::cout << "d"<<i<<"="<< Patch[D].d[i] <<"\n";
      }}    
   }}
}

void RKstep(Patches &P, 
	    Params &par, 
	    vector_complex &k,
	    double cst)
{

  // Do chebyshev D2 multiplicaation
  for(unsigned long int i=1; i<P.N;i++)
    for(unsigned long int j=1; j<P.N;j++)
      k[i*(P.N+1)+j] = complex(0.0,0.0);

  for(unsigned long int i=1; i<P.N;i++){
    for(unsigned long int j=1; j<P.N;j++){
      for(unsigned long int t=0; t<P.N+1; t++){ 
		k[i*(P.N+1)+j] += P.D[j*(P.N+1)+t] * P.tmp[i*(P.N+1)+t] ;
		k[j*(P.N+1)+i] += P.D[j*(P.N+1)+t] * P.tmp[t*(P.N+1)+i] ;
    }}}

  for(int i=1; i<P.N;i++)
  for(int j=1; j<P.N;j++)
    P.tmp[i*(P.N+1)+j] = P.psi[i*(P.N+1)+j] + k[i*(P.N+1)+j] * cst * 0.5 * par.dt * complex(0.0,1.0);

}
void AMR_Scheme(Patches &P){
  double threshold = 0.1;

  double ave = 0.0;
  for (size_t j = 1; j < P.N-1; ++j)
    ave += std::abs(P.psi[j]);

  ave /= P.N;

  if ((ave>threshold) && (P.level==0)) {
    Refinement(P);
    std::cout <<"refined?\n";
  }
  else if ((ave<threshold) && (P.level==1)) {
    Derefinement(P);
    std::cout <<"Derefined?\n";
  }
  // if ((std::abs(P.psi[j])>threshold) && (P.level==0)) {
  //   Refinement(P);
  //   std::cout <<"refined?\n";
  //   break;
  // }
  // else if ((std::abs(P.psi[j])<threshold) && (P.level==1)) {
  //   Derefinement(P);
  //   std::cout <<"Derefined?\n";
  //   break;
  // }
}

void Time_evolv(vector_patch &Patch,
		signed long int NsubD,
		Params &par,
		double tn){

  // initialize
  for (int D=0;D<NsubD*NsubD;D++){
    for (size_t j = 0; j < (Patch[D].N+1)*(Patch[D].N+1) ; j++){
      Patch[D].tmp[j] = Patch[D].psi[j]; 
  }}

  Assign_BC(Patch,NsubD,par,tn);
   // for (int D=0; D<NsubD*NsubD; D++)
   //   cheby_coeff(Patch[D].psi,Patch[D].N,Patch[D].coeff,Patch[D]);

  if (tn == 2 ){
    std::cout << "rinting?>\n";
   for (int D=0;D<NsubD*NsubD;D++){
     Print_InterpolationPatch(par,Patch,NsubD,tn-1,D);
     Print_Patch(par,Patch,NsubD,tn-1,D);
   }
// Assign_AnalyticBC(Patch,NsubD,par,tn-1);
// Print_file(par,Patch,NsubD,tn-1);
  if ((int)tn == 2 )  std::exit(0);
 }
  //  for (int D=0;D<NsubD;D++) AMR_Scheme(Patch[D]);

  // record original BC
  for (int D=0;D<NsubD*NsubD;D++)
    for(int i=0; i<(Patch[D].N+1) * (Patch[D].N+1);i++)
      Patch[D].k1[i] = complex(0.0,0.0);

  for (int D=0;D<NsubD*NsubD;D++) {
    for(int i=0; i<Patch[D].N+1;i++){
    for(int j=0; j<Patch[D].N+1;j++){

      for(unsigned long int t=0; t<Patch[D].N+1; t++){ 
	Patch[D].k1[i*(Patch[D].N+1)+j] += Patch[D].Dc[j*(Patch[D].N+1)+t] * Patch[D].tmp[i*(Patch[D].N+1)+t] ;
//	Patch[D].k1[j*(Patch[D].N+1)+i] += Patch[D].Dc[j*(Patch[D].N+1)+t] * Patch[D].tmp[t*(Patch[D].N+1)+i] ;
    }
    }}}

  for (int D=0;D<NsubD*NsubD;D++)
    for(int i=0; i<(Patch[D].N+1) * (Patch[D].N+1);i++)
      Patch[D].psi[i] += Patch[D].k1[i];

  // initialize
  for (int D=0;D<NsubD*NsubD;D++){
    for (size_t j = 0; j < (Patch[D].N+1)*(Patch[D].N+1) ; ++j){
      Patch[D].tmp[j] = Patch[D].psi[j]; 
  }}

//  Assign_CornerBC(Patch,NsubD,par,tn);

  // record original BC
  for (int D=0;D<NsubD*NsubD;D++)
    for(int i=0; i<(Patch[D].N+1) * (Patch[D].N+1);i++)
      Patch[D].k1[i] = complex(0.0,0.0);

  for (int D=0;D<NsubD*NsubD;D++) {
    for(int i=0; i<Patch[D].N+1;i++){
    for(int j=0; j<Patch[D].N+1;j++){

      for(unsigned long int t=0; t<Patch[D].N+1; t++){ 
	//	Patch[D].k1[i*(Patch[D].N+1)+j] += Patch[D].Dc[j*(Patch[D].N+1)+t] * Patch[D].tmp[i*(Patch[D].N+1)+t] ;
	Patch[D].k1[j*(Patch[D].N+1)+i] += Patch[D].Dc[j*(Patch[D].N+1)+t] * Patch[D].tmp[t*(Patch[D].N+1)+i] ;
    }
    }}}

  for (int D=0;D<NsubD*NsubD;D++)
    for(int i=0; i<(Patch[D].N+1) * (Patch[D].N+1);i++)
      Patch[D].psi[i] += Patch[D].k1[i];
  
}

int main(int argc, char **argv)
{
  unsigned long int N,Nsub;
  signed long int NsubD;
  int ghost_pts;

  if (argc != 5) {
    std::cout << "ERROR: Arguments "<< argc << "!= 5 \n";
    std::exit(0);
  }
  else{
    inp_timestep = std::atof(argv[1]);
    N = std::atoi(argv[2]) ;
    NsubD = std::atoi(argv[3]) ;
    ghost_pts= std::atoi(argv[4]) ;
    //    NsubD = std::atoi(argv[4]) ;
    //    inp_dt = 0.5*0.05/inp_timestep; // 
    inp_dt = 0.1/inp_timestep; // 

    N += ghost_pts*2;
  }
  
  //  std::cout << N+1<<"\n";
  Params par = Params(inp_xmax, inp_res, inp_dt, inp_timestep, inp_np, inp_im, inp_xoffset); // xmax, res, dt, timestep, np, im

  //  std::cout << N+1<<"\n";
  Operators opr = Operators(par, 0.0, inp_xoffset);

  //  std::cout << N+1<<"\n";
    
  double Coeff = 0.5 * par.dt / (par.dx*par.dx);

  // sub domain parameters
  //  unsigned long int Nsub = 20 ; // doesnt matter if no refinement
  double a = 0.0;
  double b = 1.0;
  double L = b-a;

  vector_patch Patch;
  
  for (int i=0;i<NsubD*NsubD;i++) {
    Patch.push_back(Patches(N));
    Patch[i].L = L;
    Patch[i].N_ghost = ghost_pts;
    Patch[i].dt = par.dt;
  }
   
  // create subdomain x
  for (int Di=0;Di<NsubD;Di++){  
    for (int Dj=0;Dj<NsubD;Dj++){
    Patch[Di*NsubD+Dj].x = cheby_exgrid_oni(N,
					    a+L*Di/NsubD,
					    a+L*(Di+1)/NsubD,
					    ghost_pts);  
    Patch[Di*NsubD+Dj].y = cheby_exgrid_oni(N,
					    a+L*Dj/NsubD,
					    a+L*(Dj+1)/NsubD,
					    ghost_pts);  
  }}

  double dx_sub = Patch[0].x[1]-Patch[0].x[0];

  // make cheby D matrix
  for (int D=0;D<NsubD*NsubD;D++)  cheby_expDn_matrix(Patch[D].D, Patch[D].Dc,
						Patch[D].x,Patch[D].N, par.dt);

  // initialize GWP
  for (int D=0; D<NsubD*NsubD;D++)   initial_GWP(Patch[D].x,Patch[D].y,Patch[D].psi);
  
//  cheby_coeff();
//  std::exit(0);
  Print_file(par,Patch,NsubD,0);

  std::cout << inp_timestep << " CFL constant = "<< par.dt/(dx_sub * dx_sub) << "\n";
  std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
  
  //    std::cout<< Patch[0].N+1 <<"\n";
  Print_file(par,Patch,NsubD,0);
  for (int tn = 1; tn <= par.timesteps; ++tn){

//    Assign_BC(Patch,NsubD,par,tn);
//    Print_Patch(par,Patch,NsubD,tn-1,0);
//    if(tn==500) std::exit(0);
    // RK4
    Time_evolv(Patch,
	       NsubD,
	       par,
	       tn);

//    Print_file(par,Patch,NsubD,tn);
//    Print_Patch(par,Patch,NsubD,tn,0);
    //    if (tn%1000==0) Print_file(par,Patch,NsubD,tn);
   std::cout << tn << "\n";
 } //   for (int tn = 0; tn <= par.timesteps; ++tn){

  double err = 0.0;
  double NN = NsubD*(N-ghost_pts*2);
  //  double NN = N + (N-1)*(NsubD-1) ;

  //  for (int D=0;D<NsubD;D++)   Print(par,Patch[D].x,Patch[D].psi,par.timesteps);
  Print_file(par,Patch,NsubD,par.timesteps);
  for (int D=0;D<NsubD*NsubD;D++)   Error(par,Patch[D].x,Patch[D].y,Patch[D].psi,err);

  std::cout << "Error =" << pow(err,1./2)/pow(NN,2) << "\n";
  std::cout << NN <<"\n";

//  Print_file(par,Patch,NsubD,par.timesteps);
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  //  std::cout << "simulation takes = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin2).count()*1e-6<< " seconds"<< std::endl;

  return 0;
}
