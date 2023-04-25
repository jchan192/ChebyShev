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
//#include <fftw3.h>
#include <iomanip>

using complex = std::complex<double>; 
using vector_real = std::vector<double>;
using vector_complex = std::vector<std::complex<double>>;
using vector_patch = std::vector<Patches>;

complex linear_interp(double x, double x0, double x1,
		      complex psi1, complex psi2){

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

void Get_cheby_coeff(vector_real x,vector_complex y, double N, vector_complex &Coeff)
{

  //  std::cout << "start getting cheby ccoeff\n"; 
  //  vector_complex a;
  complex sum;
  vector_real x_chebgrid = cheby_exgrid(N,-1.0,1.0);

  //  std::cout << "start getting cheby ccoeff 2 ..\n"; 
  Coeff.clear();
  Coeff.shrink_to_fit();
  Coeff.reserve(N+1);

  //  Coeff.resize(N+1);
  //  std::cout << x.size() << " " << y.size() << " " << N << "\n";

   for (int j=0;j<N+1;j++){
     sum = y[0]*cheby_poly(x_chebgrid[0],j) * 0.5; 

     for (int j_sum=1;j_sum<N;j_sum++)
       sum += y[j_sum]*cheby_poly(x_chebgrid[j_sum],j);    

     sum += y[N]*cheby_poly(x_chebgrid[N],j) * 0.5;
     sum *= 2.0/N;
     Coeff.push_back(sum);
     //     std::cout << j << " " << Coeff[j] << " " << N+1 << "\n";
   }
   //   std::cout << Coeff.size() <<" " ;
   //   std::cout <<"out \n";

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
complex cheby_interp(double &x, 
		     vector_complex 
		     &coeff,double N,
		     double a,
		     double b)
{
  complex interp;
  double x_chebgrid;

  x_chebgrid = ( 2.0 * x - a  - b ) / ( b-a);

    interp  = coeff[0]*cheby_poly(x_chebgrid,0) * 0.5;

    for (int j=1;j<N;j++)
      interp += coeff[j]*cheby_poly(x_chebgrid,j);

    interp += coeff[N]*cheby_poly(x_chebgrid,N) * 0.5;

  return interp ;	  

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
complex GWP_vals(double x,double Time)
{
  
  double hbar = 1.;
  double m =1.;
  double  sigma_px = 8;
  //  double  vx0 =  v0_init ; //2e1;
  double  vx0 =  -4e1 ; //2e1;
  double  x0 = x0_init;

  double sigma_x, phase_x, M_x;

  sigma_x = hbar/(2*sigma_px) * pow(1+4*pow(sigma_px,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
  phase_x = 1./hbar * (m*vx0+sigma_px*sigma_px/(sigma_x*sigma_x)*Time/(2*m)*(x-x0-vx0*Time))*
    (x-x0-vx0*Time) + vx0*m/(2*hbar) * vx0*Time - atan(2*sigma_px*sigma_px*Time/(hbar*m))/ 2 ;
  M_x = 1./(pow(2*M_PI,1./4)*pow(sigma_x,1./2)) * exp(-(x-x0-vx0*Time)*(x-x0-vx0*Time)/(4*sigma_x*sigma_x)) ;

  return complex(M_x*cos(phase_x),M_x*sin(phase_x));
}
complex GWP_GAMER_vals(double x,double Time)
{
  
  // double hbar = 1.;
  // double m =1.;
  // double  sigma_px = 8;
  // //  double  vx0 =  v0_init ; //2e1;
  // double  vx0 =  -4e1 ; //2e1;
  // double  x0 = x0_init;

  // double sigma_x, phase_x, M_x;

  // sigma_x = hbar/(2*sigma_px) * pow(1+4*pow(sigma_px,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
  // phase_x = 1./hbar * (m*vx0+sigma_px*sigma_px/(sigma_x*sigma_x)*Time/(2*m)*(x-x0-vx0*Time))*
  //   (x-x0-vx0*Time) + vx0*m/(2*hbar) * vx0*Time - atan(2*sigma_px*sigma_px*Time/(hbar*m))/ 2 ;
  // M_x = 1./(pow(2*M_PI,1./4)*pow(sigma_x,1./2)) * exp(-(x-x0-vx0*Time)*(x-x0-vx0*Time)/(4*sigma_x*sigma_x)) ;

  // return complex(M_x*cos(phase_x),M_x*sin(phase_x));

  const double Gau_v0  = 32.0;                    //mean velocity [1.0]
  const double Gau_Width =  2.0;                 // Gaussian width [0.1]
  const double Gau_Center =  4.0;                // Gaussian center [box center]
  const double r = x ;
 
 double ELBDM_ETA = 1.0;
 double Re=0,Im=0;
 double Gau_PeriodicN = 0;

 for (int n=0; n<Gau_PeriodicN+1; n++) {
 for (int m=0; m<((n==0)?1:2); m++) {
 
 const double Center     = Gau_Center + n*(1-2*m)* 8.0;
 const double dr1        = r -     Gau_v0*Time - Center;
 const double dr2        = r - 0.5*Gau_v0*Time - Center;
 const double Gau_Const1 = 1.0 + pow(  Time / ( ELBDM_ETA*pow(Gau_Width,2.) ), 2.0  );
 const double Gau_Const2 = pow( pow(Gau_Width,2.)*M_PI*Gau_Const1, -0.25 )
                            *exp(  -0.5*pow( dr1/Gau_Width, 2.0 )/Gau_Const1);
 const double Gau_Theta1 = -0.5*acos(  pow( Gau_Const1, -0.5 )  );
 const double Gau_Theta2 = 0.5*pow( dr1, 2.0 )*ELBDM_ETA*Time/(  pow( ELBDM_ETA*pow(Gau_Width,2.), 2.0) + pow(Time,2.)  )
                             + Gau_v0*ELBDM_ETA*dr2;

 Re += Gau_Const2*cos( Gau_Theta1 + Gau_Theta2 );
 Im += Gau_Const2*sin( Gau_Theta1 + Gau_Theta2 );
 }}
// return complex(Gau_Const2*cos( Gau_Theta1 + Gau_Theta2 ),Gau_Const2*sin( Gau_Theta1 + Gau_Theta2));
 return complex(Re,Im);
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
void initial_GWP(vector_real &x, vector_complex &psi)
{
  
  unsigned long int N = x.size();
  double hbar = 1.;
  double m =1.;
  double  sigma_px = 10;
  //  double  vx0 =  v0_init ; //2e1;
  double  vx0 =  -4e1 ; //2e1;
  double  x0 = x0_init;

  double Time = 0;
  double sigma_x, phase_x, M_x;


  for (size_t i = 0; i < N; ++i){

//    psi[i] = GWP_vals(x[i], Time); //imag
    psi[i] = GWP_GAMER_vals(x[i], Time); //imag
    //    std::cout << x[i]<< " " << std::abs(psi[i]) <<  "\n"; 
  }
  //  write_data(par, opr, 0);

}
void Print(Params &par,
	   vector_real &x,
	   vector_complex  &psi,
	   int tn){

  double size = x.size()-1;

  for (int i=0; i< size+1; i++){
    std::cout << x[i] << " " 
  	      << std::abs(psi[i]) << " " 
  	      << std::abs(GWP_GAMER_vals(x[i], tn*par.dt))<< " "
  	      << psi[i].real() << " " 
  	      << GWP_GAMER_vals(x[i], tn*par.dt).real()<< " " 
  	      << psi[i].imag() << " " 
  	      << GWP_GAMER_vals(x[i], tn*par.dt).imag()<< "\n"; 
  }
}
void Print_file(Params &par,
		vector_patch &P,
		unsigned long int NsubD,
		int tn){

  int size;
  int Ng = 8;

  std::stringstream filename_fdm;
  filename_fdm << "output/" << tn << ".dat";
  std::ofstream ffdm = std::ofstream(filename_fdm.str());

  for (int D=0;D<NsubD;D++){
    for(int t=0; t<P[D].N+1; t++){
      P[D].k1[t].real(0.0);
      P[D].k1[t].imag(0.0);
  }}

  for (int D=0;D<NsubD;D++){
  for(int i=0; i<P[D].N+1;i++) {
    for(int t=0; t<P[D].N+1; t++)
      P[D].k1[i] += complex(P[D].D[i*(P[D].N+1)+t]*(pow(P[D].psi[t].real(),2)+pow(P[D].psi[t].imag(),2)),
			    P[D].D[i*(P[D].N+1)+t]*(pow(GWP_GAMER_vals(P[D].x[t], tn*par.dt).real(),2)+pow(GWP_GAMER_vals(P[D].x[t], tn*par.dt).imag(),2)));
    }}

  if (ffdm)
  {
    std::stringstream data_fdm;
    for (int D=0;D<NsubD;D++){
      //      for (size_t i = Ng; i < P[D].N-Ng; ++i){
      for (size_t i = Ng; i < P[D].N+1-Ng; ++i){
	data_fdm << std::setprecision(14)
		 << P[D].x[i] << " " 
  	      << std::abs(P[D].psi[i]) << " " 
  	      << std::abs(GWP_GAMER_vals(P[D].x[i], tn*par.dt))<< " "
  	      << P[D].psi[i].real() << " " 
	      << GWP_GAMER_vals(P[D].x[i], tn*par.dt).real()<< " " 
  	      << P[D].psi[i].imag() << " " 
	      << GWP_GAMER_vals(P[D].x[i], tn*par.dt).imag()<< " " 
  	      << P[D].k1[i].real() << " "
  	      << P[D].k1[i].imag() << "\n"; 

      }
    }
      ffdm.write(data_fdm.str().c_str(), data_fdm.str().length());
  }
  ffdm.close();

}

void Error(Params &par,
	   vector_real &x,
	   vector_complex  &psi,
	   double &Error){
  
  int Ng = 8;
  double size = x.size()-1;
  for (size_t i = Ng; i < size-Ng; ++i){
    Error += pow(std::abs(psi[i])
		 -std::abs(GWP_GAMER_vals(x[i], par.timesteps*par.dt)),2);
  }
  //  error = pow(error,1./2)/opr.size;  
  //  std::cout << "error = " << error << "\n";    
}
void RK2(int N, vector_real &D,vector_complex &psi, Params &par){

  vector_complex k1(N+1),tmp(N+1);
    for(int i=1; i<N;i++) {
      k1[i] = complex(0.0,0.0);
      //      for(int j=0; j<N+1; j++) k1[i] += D[i*(N+1)+j] * psi[j] ;
            for(int j=1; j<N; j++) k1[i] += D[i*(N+1)+j] * psi[j] ;
    }
    for(int i=1; i<N;i++) k1[i] = psi[i] + k1[i] * 0.25 * par.dt * complex(0.0,1.0);

    for(int i=1; i<N;i++) {
      tmp[i] = complex(0.0,0.0);
      //for(int j=0; j<N+1; j++) tmp[i] += D[i*(N+1)+j] * k1[j] ;
      for(int j=1; j<N; j++) tmp[i] += D[i*(N+1)+j] * k1[j] ;
    }
    for(int i=1; i<N;i++) psi[i] += tmp[i] * 0.5 * par.dt * complex(0.0,1.0);
}
void Refinement(Patches &P){
  
  // Perform Refinement: double cheby grid size

  int newN = 20;
  // for (int i=0; i< P.N+1; i++){
  //   std::cout << P.x[i] << " " 
  // 	      << P.psi[i].real() << " " 
  // 	      << P.psi[i].imag() <<"\n";
  // }
  
  //  std::cout << "Start refinement\n" << "\n";
  
  double a_old = P.x[0];
  double b_old = P.x.back();

  Get_cheby_coeff(P.x,P.psi,P.N, P.coeff);
  //  std::cout << "Start refinement 3 ...\n" << "\n";
  P.x  = cheby_exgrid_oni(newN+P.N_ghost*2,
			  P.x[P.N_ghost],
			  P.x[P.N-P.N_ghost+1],
			  P.N_ghost);
  P.psi = cheby_interp(P.x,P.coeff,P.N,a_old,b_old); 

  P.level = 1;

  //  P.N = 2*(P.N-P.N_ghost*2)+P.N_ghost*2; 
  P.N = newN+P.N_ghost*2;
  P.tmp.resize(P.N+1);
  P.k1.resize(P.N+1);
  cheby_expDn_matrix(P.D, 
		     P.Dc,
		     P.x,P.N, P.dt);

  //    std::cout << "done refinement" << "\n";
  //    std::cout << P.N << " " << P.x.size() << " " << P.Dc.size() << "\n";

  // for (int i=0; i< P.N+1; i++){
  //   std::cout << P.x[i] << " " 
  // 	      << P.psi[i].real() << " " 
  // 	      << P.psi[i].imag() <<" "
  // 	      << GWP_vals(P.x[i],0).real() << " "
  // 	      << GWP_vals(P.x[i],0).imag() << "\n";
  // }

  //  std::exit(0);
}

void Derefinement(Patches &P){
  
  // Perform De-Refinement: half cheby grid size

  // for (int i=0; i< P.N+1; i++){
  //   std::cout << P.x[i] << " " 
  // 	      << P.psi[i].real() << " " 
  // 	      << P.psi[i].imag() <<"\n";
  // }
  double half = 0.5 ; 
  Get_cheby_coeff(P.x,P.psi,P.N, P.coeff);
  //  P.x  = cheby_exgrid_oni(half*P.N,P.x[0],P.x[P.N]);
  P.psi = cheby_interp(P.x,P.coeff,P.N,P.x[0],P.x[half*P.N]); 
  //  cheby_DD_matrix(P.D, P.x, half*P.N);
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

  int ghost_pts = Patch[0].N_ghost;
  int N;
  //  vector_complex ;
  //  std::cout << "start Coeby Coeff?\n"; 

  // get Cheby Coeffs
  for (int D=0;D<NsubD;D++) {
    Get_cheby_coeff(Patch[D].x,Patch[D].psi,Patch[D].N,Patch[D].coeff);
  }

  // assign BC
  for (int D=0;D<NsubD;D++) {

    N = Patch[D].N;

    if (D==0){
      for (int ig=0; ig< ghost_pts; ig++){
	Patch[D].psi[ig] = GWP_GAMER_vals(Patch[D].x[ig],(tn-1.)*par.dt);
	Patch[D].tmp[ig] = GWP_GAMER_vals(Patch[D].x[ig],(tn-1.)*par.dt);
      }
    }
    //    else if (Patch[D-1].level==0){
    else {
      for (int ig=0; ig< ghost_pts; ig++){
	Patch[D].psi[ig] =  cheby_interp(Patch[D].x[ig],
					 Patch[D-1].coeff,
					 Patch[D-1].N, 
					 Patch[D-1].x[0],
					 Patch[D-1].x.back());
	Patch[D].tmp[ig] =  Patch[D].psi[ig];
      }
    }
    // else if (Patch[D-1].level==1){
    //   for (int ig=0; ig< ghost_pts; ig++){
    // 	Patch[D].psi[ig] =  cheby_interp(Patch[D].x[ig],
    // 					 Patch[D-1].coeff,
    // 					 Patch[D-1].N, 
    // 					 Patch[D-1].x[0],
    // 					 Patch[D-1].x.back());
    // 	Patch[D].tmp[ig] =  Patch[D].psi[ig];
    //   }
    // }

    if (D==(NsubD-1)){
      for (int ig=0; ig< ghost_pts; ig++){
	Patch[D].psi[N-ig] = GWP_GAMER_vals(Patch[D].x[N-ig],(tn-1.)*par.dt);
	Patch[D].tmp[N-ig] = GWP_GAMER_vals(Patch[D].x[N-ig],(tn-1.)*par.dt);
      }
    }
    //    else if (Patch[D+1].level==0){
    else {
      for (int ig=0; ig< ghost_pts; ig++){
	Patch[D].psi[N-ig] =  cheby_interp(Patch[D].x[N-ig],
					   Patch[D+1].coeff,
					   Patch[D+1].N, 
					   Patch[D+1].x[0],
					   Patch[D+1].x.back());
	Patch[D].tmp[N-ig] =  Patch[D].psi[N-ig];
      }
    }
    // else if (Patch[D+1].level==1){
    //   for (int ig=0; ig< ghost_pts; ig++){
    // 	Patch[D].psi[N-ig] =  cheby_interp(Patch[D].x[N-ig],
    // 					   Patch[D+1].coeff,
    // 					   Patch[D+1].N, 
    // 					   Patch[D+1].x[0],
    // 					   Patch[D+1].x.back());
    // 	Patch[D].tmp[N-ig] =  Patch[D].psi[N-ig];
    //   }
    // }
    //    std::cout << Patch[D].a  << " " << Patch[D].b <<"\n";
    
  }


  
  //    std::exit(0);
}

void RKstep(Patches &P, 
	    Params &par, 
	    vector_complex &k,
	    double cst,
	    complex aold,
	    complex bold)
{

  for(int i=1; i<P.N;i++){
    k[i] = complex(0.0,0.0);
    for(int j=0; j<P.N+1; j++) k[i] += P.D[i*(P.N+1)+j] * P.tmp[j] ;
  }
  for(int i=1; i<P.N;i++) P.tmp[i] = P.psi[i] + k[i] * cst * 0.5 * par.dt * complex(0.0,1.0);

}
void AMR_Scheme(Patches &P){
  double threshold = 0.1;

  double ave = 0.0;
  for (size_t j = 1; j < P.N-1; ++j)
    ave += std::abs(P.psi[j]);

  ave /= P.N;

  if ((ave>threshold) && (P.level==0)) {
    std::cout <<"before refined?\n";
    Refinement(P);
    std::cout <<"refined?\n";
  }
  // else if ((ave<threshold) && (P.level==1)) {
  //   Derefinement(P);
  //   std::cout <<"Derefined?\n";
  // }
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

  vector_complex a(NsubD),b(NsubD);
  complex Coeff = complex(0.0,1.0) * 0.5 * par.dt;
  int Ng = Patch[0].N_ghost;

  //  std::cout << tn <<"\n";
  // initialize

  Assign_BC(Patch,NsubD,par,tn);

  // for (int D=0;D<NsubD;D++) 
  //   if (D==2 || D==3) AMR_Scheme(Patch[D]);

  for (int D=0;D<NsubD;D++){
    for (size_t j = 0; j < Patch[D].N+1; ++j){
      Patch[D].tmp[j] = Patch[D].psi[j]; 
      //      std::cout << std::abs(Patch[D].tmp[j]) << "\n";
  }}

  // record original BC
  for (int D=0;D<NsubD;D++) {
    for(int i=0; i<Patch[D].N+1;i++){
      Patch[D].k1[i] = complex(0.0,0.0);

      for(int j=0; j<Patch[D].N+1; j++) Patch[D].k1[i] += (Patch[D].Dc[i*(Patch[D].N+1)+j]) *  Patch[D].tmp[j];

      Patch[D].psi[i] += Patch[D].k1[i]; 
  }}
  

  // std::exit(0);
}

int main(int argc, char **argv)
{
  unsigned long int N,Nsub;
  signed long int NsubD;
  int ghost_pts ; 

  if (argc != 5) {
    std::cout << "ERROR: Arguments "<< argc << "!= 6 \n";
    std::exit(0);
  }
  else{
    inp_timestep = std::atof(argv[1]);
    N = std::atoi(argv[2]) ;
    NsubD = std::atoi(argv[3]) ;
    ghost_pts= std::atoi(argv[4]) ;
    //    NsubD = std::atoi(argv[4]) ;
    //    inp_dt = 0.5*0.05/inp_timestep; // 
    inp_dt = 1e-2/inp_timestep; // 

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
  double b = 8.0;

  complex a1,b1,a2,b2,a3,b3,a4,b4;
  //  double L = par.xmax/((double) NsubD+(1-NsubD)*(1.0-cos(M_PI/N))/2.);
  double L = b-a;
  vector_patch Patch;
  
  for (int i=0;i<NsubD;i++) {
    Patch.push_back(Patches(N));
    Patch[i].L = L;
    Patch[i].N_ghost = ghost_pts;
    Patch[i].dt = par.dt;
  }

  // create subdomain x
  for (int D=0;D<NsubD;D++)  Patch[D].x = cheby_exgrid_oni(N,
							   a+L*D/NsubD,
							   a+L*(D+1)/NsubD,
							   ghost_pts);  
  double dx_sub = Patch[0].x[1]-Patch[0].x[0];

  // make cheby D matrix
  for (int D=0;D<NsubD;D++)  cheby_expDn_matrix(Patch[D].D, Patch[D].Dc,
						Patch[D].x,Patch[D].N, par.dt);

  // initialize GWP
  for (int D=0;D<NsubD;D++)   initial_GWP(Patch[D].x,Patch[D].psi);

  std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();

  Print_file(par,Patch,NsubD,0);

  for (int tn = 1; tn <= par.timesteps; ++tn){

    // RK4
    Time_evolv(Patch,
    	       NsubD,
    	       par,
    	       tn);

    //    dx_sub = Patch[2].x[1]-Patch[2].x[0];
//    Print_file(par,Patch,NsubD,tn);
    //    if (tn==1) std::exit(0);
    //    if (tn%2==0) Print_file(par,Patch,NsubD,tn);
    //    std::exit(0);
    //    if (tn==100) std::exit(0);
 } //   for (int tn = 0; tn <= par.timesteps; ++tn){

  double err = 0.0;
  //  double NN = N + (N-1)*(NsubD-1) ;
  double NN = NsubD*(N-ghost_pts*2);

  //  for (int D=0;D<NsubD;D++)   Print(par,Patch[D].x,Patch[D].psi,par.timesteps);
  Print_file(par,Patch,NsubD,par.timesteps);
  for (int D=0;D<NsubD;D++)   Error(par,Patch[D].x,Patch[D].psi,err);

  std::cout << "Error =" << pow(err,1./2)/NN << "\n";
  std::cout << NN <<"\n";
  std::cout << "CFL constant = "<< par.dt/(dx_sub * dx_sub) << "\n";

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "simulation takes = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin2).count()*1e-6<< " seconds"<< std::endl;

return 0;
}
