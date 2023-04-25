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
//include <fftw3.h>
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
complex cheby_interp(double &x, vector_complex &coeff,double N,double a,double b)
{
  //  int Nx = x.size();
  complex interp;
  double x_chebgrid;

    x_chebgrid = ( 2.0 * x - a  - b ) / ( b-a);

    interp  = coeff[0]*cheby_poly(x_chebgrid,0) * 0.5;

    for (int j=1;j<N;j++)
      interp += coeff[j]*cheby_poly(x_chebgrid,j);

    interp += coeff[N]*cheby_poly(x_chebgrid,N) * 0.5;

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
// void fft_r2c(vector_real &x, vector_complex &y, bool forward)
// {
//   fftw_complex *out = reinterpret_cast<fftw_complex*>(y.data());  // standard declaration of output fttw
//   fftw_plan p;

//   if (forward)
//   {
//     p = fftw_plan_dft_r2c_1d(x.size(), x.data(), out, FFTW_ESTIMATE);
//     fftw_execute(p);
//     fftw_destroy_plan(p);
//     for (size_t i = 0; i < x.size(); ++i)
//     {
//       y[i] = y[i]/sqrt(x.size());
//     }
//   }

//   else 
//   {
//     p = fftw_plan_dft_c2r_1d(x.size(), out, x.data(), FFTW_ESTIMATE);
//     fftw_execute(p);
//     fftw_destroy_plan(p);
//     for (size_t i = 0; i < x.size(); ++i)
//     {
//       x[i] = x[i]/sqrt(x.size());
//     }
//   }
// }  

// void fft(vector_complex &x, bool inverse)
// {
//   vector_complex y(x.size(), complex(0.0, 0.0));
//   fftw_plan p;

//   fftw_complex *in = reinterpret_cast<fftw_complex*>(x.data());   // standard declaration of input fttw
//   fftw_complex *out = reinterpret_cast<fftw_complex*>(y.data());  // standard declaration of output fttw

//   p = fftw_plan_dft_1d(x.size(), in, out, (inverse ? FFTW_BACKWARD : FFTW_FORWARD), FFTW_ESTIMATE);

//   fftw_execute(p);
//   fftw_destroy_plan(p);

//   for (size_t i = 0; i < x.size(); ++i)
//   {
//     x[i] = y[i] / sqrt(static_cast<double>(x.size()));   // fftw give unnormalized output, needs renormalization here.
//   }
// } 
complex GWP_vals(double x,double Time)
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
 double Gau_PeriodicN = 50;

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

    psi[i] = GWP_vals(x[i], Time); //imag
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
  	      << std::abs(GWP_vals(x[i], tn*par.dt))<< " "
  	      << psi[i].real() << " " 
  	      << GWP_vals(x[i], tn*par.dt).real()<< " " 
  	      << psi[i].imag() << " " 
  	      << GWP_vals(x[i], tn*par.dt).imag()<< "\n"; 
  }


}
void Print_file(Params &par,
		vector_patch &P,
		unsigned long int NsubD,
		int tn){

  int size;

  std::stringstream filename_fdm;
  filename_fdm << "output/" << tn << ".dat";
  std::ofstream ffdm = std::ofstream(filename_fdm.str());

  if (ffdm)
  {
    std::stringstream data_fdm;
    for (int D=0;D<NsubD;D++){
      for (int i = 1; i < P[D].N-1; ++i){
	  data_fdm << P[D].x[i] << " " 
  	      << std::abs(P[D].psi[i]) << " " 
  	      << std::abs(GWP_vals(P[D].x[i], tn*par.dt))<< " "
  	      << P[D].psi[i].real() << " " 
	      << GWP_vals(P[D].x[i], tn*par.dt).real()<< " " 
  	      << P[D].psi[i].imag() << " " 
	      << GWP_vals(P[D].x[i], tn*par.dt).imag()<< "\n"; 
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
  
  double size = x.size()-1;
  for (size_t i = 1; i < size; ++i){
    Error += pow(std::abs(psi[i])
		 -std::abs(GWP_vals(x[i], par.timesteps*par.dt)),2);
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

  // for (int i=0; i< P.N+1; i++){
  //   std::cout << P.x[i] << " " 
  // 	      << P.psi[i].real() << " " 
  // 	      << P.psi[i].imag() <<"\n";
  // }
  P.coeff = cheby_coeff(P.x,P.psi,P.N);
  P.x  = cheby_exgrid(2*P.N,P.x[0],P.x[P.N]);
  P.psi = cheby_interp(P.x,P.coeff,P.N,P.x[0],P.x[2*P.N]); 
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
  P.psi = cheby_interp(P.x,P.coeff,P.N,P.x[0],P.x[half*P.N]); 
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

  for (int D=0;D<NsubD;D++) {
    if (D==0){
      Patch[D].a = GWP_vals(Patch[D].x[0],(tn-1.)*par.dt); 
      Patch[D].psi[0] = Patch[D].a;
    }
    else if (Patch[D-1].level==0){
      Patch[D].a = Patch[D-1].tmp[Patch[D-1].N-1]; 
      Patch[D].psi[0] = Patch[D].a;
    }
    else if (Patch[D-1].level==1){
      Patch[D].a = Patch[D-1].tmp[Patch[D-1].N-2]; 
      Patch[D].psi[0] = Patch[D].a;
    }

    if (D==(NsubD-1)){
      Patch[D].b = GWP_vals(Patch[D].x[Patch[D].N],(tn-1.)*par.dt); 
      Patch[D].psi[Patch[D].N] = Patch[D].b;
    }
    else if (Patch[D+1].level==0){
      Patch[D].b = Patch[D+1].tmp[1]; 
      Patch[D].psi[Patch[D].N] = Patch[D].b;
    }
    else if (Patch[D+1].level==1){
      Patch[D].b = Patch[D+1].tmp[2]; 
      Patch[D].psi[Patch[D].N] = Patch[D].b;
    }
    //    std::cout << Patch[D].a  << " " << Patch[D].b <<"\n";
  }


  
  //  std::exit(0);
}

void RKstep(Patches &P, 
	    Params &par, 
	    vector_complex &k,
	    double cst,
	    complex aold,
	    complex bold)
{
  // transform homogeneous problem (based on new BC)
  for (size_t j = 0; j < P.N+1; ++j)
    P.tmp[j]  -= ((1.-(P.x[j] - P.x[0])/P.L)* P.a + (P.x[j]-P.x[0])/P.L* P.b);

  for(int i=1; i<P.N;i++){
    k[i] = complex(0.0,0.0);
    for(int j=1; j<P.N; j++) k[i] += P.D[i*(P.N+1)+j] * P.tmp[j] ;
  }
  for(int i=1; i<P.N;i++) P.tmp[i] = P.psi[i] + k[i] * cst * 0.5 * par.dt * complex(0.0,1.0);

  // de-transform (based on old BC) since its P.psi += RK
  for (size_t j = 0; j < P.N+1; ++j){
    P.tmp[j]  += ((1.-(P.x[j] - P.x[0])/P.L)* aold + (P.x[j]-P.x[0])/P.L* bold);  
    //    std::cout << P.tmp[j] <<"\n";
  }

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

  vector_complex a(NsubD),b(NsubD);


  // initialize
  for (int D=0;D<NsubD;D++){
    for (size_t j = 0; j < Patch[D].N+1; ++j){
      Patch[D].tmp[j] = Patch[D].psi[j]; 
  }}

  // transform homogeneous problem
  Assign_BC(Patch,NsubD,par,tn);

  //  for (int D=0;D<NsubD;D++) AMR_Scheme(Patch[D]);
  //  AMR_Scheme(Patch[6]);
  //  std::cout <<"tn="<<tn <<"\n";
  for (int D=0;D<NsubD;D++){
    for (size_t j = 0; j < Patch[D].N+1; ++j){
      Patch[D].psi[j]  -= ((1.-(Patch[D].x[j] - Patch[D].x[0])/Patch[D].L)* Patch[D].a + (Patch[D].x[j]-Patch[D].x[0])/Patch[D].L* Patch[D].b); 
  }}

  // record original BC
  for (int D=0;D<NsubD;D++){
    a[D] = Patch[D].a;
    b[D] = Patch[D].b;
  } 

  // Runge kutta 4
  for (int D=0;D<NsubD;D++)   RKstep(Patch[D],par,Patch[D].k1,0.5,a[D],b[D]);

  Assign_BC(Patch,NsubD,par,tn+0.5);
  for (int D=0;D<NsubD;D++)   RKstep(Patch[D],par,Patch[D].k2,0.5,a[D],b[D]);

  Assign_BC(Patch,NsubD,par,tn+0.5);
  for (int D=0;D<NsubD;D++)   RKstep(Patch[D],par,Patch[D].k3,1.0,a[D],b[D]);

  Assign_BC(Patch,NsubD,par,tn+1.0);
  for (int D=0;D<NsubD;D++)   RKstep(Patch[D],par,Patch[D].k4,1.0,a[D],b[D]);

  for (int D=0;D<NsubD;D++){
    for (unsigned long int i = 1; i < Patch[D].N; ++i){
	Patch[D].psi[i] += 1./6 * (Patch[D].k1[i] + 2.0*Patch[D].k2[i] + 2.0*Patch[D].k3[i] + Patch[D].k4[i]) * 0.5 * par.dt * complex(0.0,1.0);
       //      std::cout << Patch[D].psi[i] << "\n";
  }}
  
  // // de-transform
    for (int D=0;D<NsubD;D++){
      for (size_t j = 0; j < Patch[D].N+1; ++j){
	Patch[D].psi[j]  += ((1.-(Patch[D].x[j] - Patch[D].x[0])/Patch[D].L)* a[D] + (Patch[D].x[j]-Patch[D].x[0])/Patch[D].L*b[D]); 
  }}
}

int main(int argc, char **argv)
{
  unsigned long int N,Nsub;
  signed long int NsubD;

  if (argc != 4) {
    std::cout << "ERROR: Arguments "<< argc << "!= 5 \n";
    std::exit(0);
  }
  else{
    inp_timestep = std::atof(argv[1]);
    N = std::atoi(argv[2]) ;
    NsubD = std::atoi(argv[3]) ;
    //    NsubD = std::atoi(argv[4]) ;
    //    inp_dt = 0.5*0.05/inp_timestep; // 
    inp_dt = 0.25*0.05/inp_timestep; // 
  }
  
  //  std::cout << N+1<<"\n";
  Params par = Params(inp_xmax, inp_res, inp_dt, inp_timestep, inp_np, inp_im, inp_xoffset); // xmax, res, dt, timestep, np, im

  //  std::cout << N+1<<"\n";
  Operators opr = Operators(par, 0.0, inp_xoffset);

  //  std::cout << N+1<<"\n";
    
  double Coeff = 0.5 * par.dt / (par.dx*par.dx);

  // sub domain parameters
  //  unsigned long int Nsub = 20 ; // doesnt matter if no refinement
  complex a1,b1,a2,b2,a3,b3,a4,b4;
  double L = par.xmax/((double) NsubD+(1-NsubD)*(1.0-cos(M_PI/N))/2.);
  vector_patch Patch;
  
  for (int i=0;i<NsubD;i++) {
    Patch.push_back(Patches(N));
    Patch[i].L = L;
  }

  // create subdomain x
  //  std::cout << N+1<<"\n";
  for (size_t i = 0; i < N+1; ++i)
    for (int D=0;D<NsubD;D++) 
      Patch[D].x[i] = L/2.0*(-cos(i*M_PI/N)) + L/2.0 ;

  double dx_sub = Patch[0].x[1]-Patch[0].x[0];

  for (size_t i = 0; i < N+1; ++i){
    for (int D=0;D<NsubD;D++) Patch[D].x[i] += D * (L-dx_sub);
  }

  // make cheby D matrix
  for (int D=0;D<NsubD;D++)  cheby_DD_matrix(Patch[D].D, Patch[D].x,Patch[D].N);

  // initialize GWP
  for (int D=0;D<NsubD;D++)   initial_GWP(Patch[D].x,Patch[D].psi);

  //  std::exit(0);
  std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();

  std::cout << "CFL constant = "<< par.dt/(dx_sub * dx_sub) << "\n";
  Print_file(par,Patch,NsubD,0);
  for (int tn = 1; tn <= par.timesteps; ++tn){

    // RK4
    Time_evolv(Patch,
	       NsubD,
	       par,
	       tn);

    //    Print_file(par,Patch,NsubD,tn);
    if (tn%1000==0) Print_file(par,Patch,NsubD,tn);
    //    if (tn==100) std::exit(0);
 } //   for (int tn = 0; tn <= par.timesteps; ++tn){

  double err = 0.0;
  double NN = N + (N-1)*(NsubD-1) ;

  //  for (int D=0;D<NsubD;D++)   Print(par,Patch[D].x,Patch[D].psi,par.timesteps);
  Print_file(par,Patch,NsubD,par.timesteps);
  for (int D=0;D<NsubD;D++)   Error(par,Patch[D].x,Patch[D].psi,err);

  std::cout << "Error =" << pow(err/NN,1./2) << "\n";
  std::cout << NN <<"\n";

  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  //  std::cout << "simulation takes = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin2).count()*1e-6<< " seconds"<< std::endl;

return 0;
}
