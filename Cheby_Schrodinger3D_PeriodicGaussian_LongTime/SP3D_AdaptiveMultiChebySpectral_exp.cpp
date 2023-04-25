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
//   for(int i=0; i<N+1;i++) {
//     for(int j=0; j<N+1; j++){
//       dX[i*(N+1)+j] = 0.0;
//       for(int k=0; k<N+1; k++){
//   	dX[i*(N+1)+j] += D4[i*(N+1)+k].real()* D2[k*(N+1)+j].real() ;
//       }
//     }
//   }
//   for(int i=0; i<(N+1)*(N+1);i++) {
//     Dc[i] += dX[i] * Coeff*Coeff*Coeff /6.0;
//     D6[i] = dX[i];
//   }

//   // // make D^8
//   for(int i=0; i<N+1;i++) {
//     for(int j=0; j<N+1; j++){
//       dX[i*(N+1)+j] = 0.0;
//       for(int k=0; k<N+1; k++){
//   	dX[i*(N+1)+j] += D6[i*(N+1)+k].real()* D2[k*(N+1)+j].real() ;
//       }
//     }
//   }

//   for(int i=0; i<(N+1)*(N+1);i++) {
// //    Dc[i] += dX[i] * Coeff*Coeff*Coeff*Coeff / 24.0;
//     // D8[i] = dX[i] * Coeff*Coeff*Coeff*Coeff / 24.0;
//     // D6[i] *= Coeff*Coeff*Coeff / 6.0;
//     // D4[i] *= 0.5 * Coeff*Coeff;
//     // D2[i] *= Coeff;
//   }

}

vector_complex cheby_1Dcoeff(vector_real x,vector_complex y, double N)
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
void cheby_2Dcoeff(vector_complex &y, double N, vector_complex &Coeff,Patches &P)
{
  
  vector_real coeff_R((N+1)*(N+1));
  vector_real coeff_I((N+1)*(N+1));

  vector_real x = cheby_exgrid(N,0,1);

  for (int j=0;j<(N+1)*(N+1);j++){
    P.coeff2D_R[j] = y[j].real();
    P.coeff2D_I[j] = y[j].imag();
  }  

  fftw_execute(P.plan2D_R);
  fftw_execute(P.plan2D_I);

  for (int j=0;j<(N+1)*(N+1);j++){
    Coeff[j].real(P.coeff2D_R[j]) ;
    Coeff[j].imag(P.coeff2D_I[j]) ;
  }
  
  // re-scale coefficients
 for (int j=0;j<(N+1)*(N+1) ;j++)
   Coeff[j] /= N*N;

  for (int j=0;j<(N+1);j++){
    Coeff[0*(N+1)+j] *= 0.5;
    Coeff[j*(N+1)+0] *= 0.5;
  }
  for (int j=0;j<N+1;j++){
    Coeff[N*(N+1)+j] *= 0.5;
    Coeff[j*(N+1)+N] *= 0.5;
  }
}

void cheby_3Dcoeff(vector_complex &y, double N, vector_complex &Coeff,Patches &P)
{
  
  vector_real coeff_R((N+1)*(N+1)*(N+1));
  vector_real coeff_I((N+1)*(N+1)*(N+1));

//   for (int i=0;i<N+1 ;i++){
//   for (int j=0;j<N+1 ;j++){
//   for (int k=0;k<N+1 ;k++){
//     coeff[k+(N+1)*(i*(N+1)+j)] = x[i]*x[i] + x[j]*x[j] + x[k]*x[k] ;
// //    std::cout << y[i] <<"\n";
//   }}

  // vector_real x = cheby_exgrid(N,0,1);

  for (int j=0;j<(N+1)*(N+1)*(N+1);j++){
    P.coeff3D_R[j] = y[j].real();
    P.coeff3D_I[j] = y[j].imag();
  }  

  fftw_execute(P.plan3D_R);
  fftw_execute(P.plan3D_I);

  for (int j=0;j<(N+1)*(N+1)*(N+1);j++){
    Coeff[j].real(P.coeff3D_R[j]) ;
    Coeff[j].imag(P.coeff3D_I[j]) ;
  }
  
  // re-scale coefficients
 for (int j=0;j<(N+1)*(N+1)*(N+1) ;j++)
   Coeff[j] /= N*N*N;

  for (int i=0;i<(N+1);i++){
  for (int j=0;j<(N+1);j++){
    Coeff[j+(N+1)*(i+0*(N+1))] *= 0.5;
    Coeff[j+(N+1)*(0+i*(N+1))] *= 0.5;
    Coeff[0+(N+1)*(j+i*(N+1))] *= 0.5;
  }}

  for (int i=0;i<(N+1);i++){
  for (int j=0;j<(N+1);j++){
    Coeff[j+(N+1)*(i+N*(N+1))] *= 0.5;
    Coeff[j+(N+1)*(N+i*(N+1))] *= 0.5;
    Coeff[N+(N+1)*(j+i*(N+1))] *= 0.5;
  }}

 //  double x_chebgrid = -( 2.0 * x[1] - 0  -1 ) / ( 1);
 //  double y_chebgrid = -( 2.0 * x[1] - 0  -1 ) / ( 1);

 //  double interp = 0.0 ;
 //  for (int i=0;i<N+1;i++)
 //  for (int j=0;j<N+1;j++)
 //      interp += Coeff[i*(N+1)+j].real()*cheby_poly(x_chebgrid,j)*cheby_poly(y_chebgrid,i);

 // std::cout << interp <<"\n";
 // std::exit(0);
}
complex cheby_3Dinterp(double &x, double &y,double &z, 
		       vector_complex &coeff,double N,
		       double ax,double bx, double ay, 
		       double by, double az, double bz)
{
  //  int Nx = x.size();
  complex interp;
  double x_chebgrid,y_chebgrid, z_chebgrid;
  double Tx0, Ty0, Tz0, Tx1, Ty1,Tz1;
  vector_real polyx(N+1),polyy(N+1),polyz(N+1);

  x_chebgrid = -( 2.0 * x - ax  - bx ) / ( bx-ax);
  y_chebgrid = -( 2.0 * y - ay  - by ) / ( by-ay);
  z_chebgrid = -( 2.0 * z - az  - bz ) / ( bz-az);

  for (int i=0;i<N+1;i++){
    polyx[i] = cheby_poly(x_chebgrid,i);
    polyy[i] = cheby_poly(y_chebgrid,i);
    polyz[i] = cheby_poly(z_chebgrid,i);
  }    

  interp = complex(0.0,0.0) ;
  for (int i=0;i<N+1;i++)
  for (int j=0;j<N+1;j++)
  for (int k=0;k<N+1;k++)
      interp += coeff[k+(N+1)*(j+i*(N+1))]*polyy[j]*polyx[i]*polyz[k];
//      interp += coeff[k+(N+1)*(j+i*(N+1))]*polyy[j];
//      interp += coeff[i*(N+1)+j]*cheby_poly(y_chebgrid,j)*cheby_poly(x_chebgrid,i);
//      interp += coeff[i*(N+1)+j]*polyy[j]*polyx[i];

  return interp ;	  
}
void cheby_3Dinterp_Optimize(vector_patch &Patch,
			     signed long int NsubD,
			     int Di,
			     int Dj,
			     int Dk)
{
  //  int Nx = x.size();
  int Dfrom, Dto1, Dto2, Dto3, Dto4,Dto5, Dto6,Dto7,Dto8;
  int N, ghost_pts;
  complex interp1,interp2,interp3,interp4,interp5,interp6,interp7,interp8;
  double Tx0, Ty0, Tz0, Tx1, Ty1,Tz1;
  vector_real polyx(N+1),polyy(N+1),polyz(N+1);
  vector_real polyxyz((N+1)*(N+1)*(N+1),0.0);
  complex cTxyz;
  double ax, ay, az, bx, by,bz;
 
  
  Dfrom = Dk + NsubD*(Dj+(Di)*NsubD);
  N = Patch[Dfrom].N;
  ghost_pts = Patch[Dfrom].N_ghost;
  vector_real x_chebgrid(ghost_pts),y_chebgrid(ghost_pts),z_chebgrid(ghost_pts);

  Dto1  = Dk + 1 + NsubD*(Dj + 1+(Di + 1)*NsubD);
  Dto2  = Dk - 1 + NsubD*(Dj + 1+(Di + 1)*NsubD);
  Dto3  = Dk + 1 + NsubD*(Dj - 1+(Di + 1)*NsubD);
  Dto4  = Dk + 1 + NsubD*(Dj + 1+(Di - 1)*NsubD);
  Dto5  = Dk - 1 + NsubD*(Dj - 1+(Di + 1)*NsubD);
  Dto6  = Dk + 1 + NsubD*(Dj - 1+(Di - 1)*NsubD);
  Dto7  = Dk - 1 + NsubD*(Dj + 1+(Di - 1)*NsubD);
  Dto8  = Dk - 1 + NsubD*(Dj - 1+(Di - 1)*NsubD);

  ax = Patch[NsubD*NsubD*NsubD-1].x[0];
  ay = Patch[NsubD*NsubD*NsubD-1].y[0];
  az = Patch[NsubD*NsubD*NsubD-1].z[0];
  bx = Patch[NsubD*NsubD*NsubD-1].x.back();
  by = Patch[NsubD*NsubD*NsubD-1].y.back();
  bz = Patch[NsubD*NsubD*NsubD-1].z.back();
    
  cheby_3Dcoeff(Patch[Dfrom].psi,Patch[Dfrom].N,Patch[Dfrom].coeff,Patch[Dfrom]);  

  for (int ig=0; ig< ghost_pts; ig++){
    x_chebgrid[ig] = -( 2.0 * Patch[NsubD - 2 + NsubD*(NsubD - 2+(NsubD - 2)*NsubD)].x[N-ig] - ax  - bx ) / ( bx-ax);
    y_chebgrid[ig] = -( 2.0 * Patch[NsubD - 2 + NsubD*(NsubD - 2+(NsubD - 2)*NsubD)].y[N-ig] - ay  - by ) / ( by-ay);
    z_chebgrid[ig] = -( 2.0 * Patch[NsubD - 2 + NsubD*(NsubD - 2+(NsubD - 2)*NsubD)].z[N-ig] - az  - bz ) / ( bz-az);
  }
  
  for (int ig=0; ig< ghost_pts; ig++){
  for (int jg=0; jg< ghost_pts; jg++){
  for (int kg=0; kg< ghost_pts; kg++){


      for (int i=0;i<N+1;i++){
	polyx[i] = cheby_poly(x_chebgrid[ig],i);
	polyy[i] = cheby_poly(y_chebgrid[jg],i);
	polyz[i] = cheby_poly(z_chebgrid[kg],i);
      }    

//      std::exit(0);
      for (int i=0;i<N+1;i++){
      for (int j=0;j<N+1;j++){
      for (int k=0;k<N+1;k++){
//	polyxyz.push_back(polyx[i]*polyy[j]*polyz[k]);
	polyxyz[k+(N+1)*(j+i*(N+1))] = polyx[i]*polyy[j]*polyz[k];
      }}}

      interp1 = complex(0.0,0.0);
      interp2 = complex(0.0,0.0);
      interp3 = complex(0.0,0.0);
      interp4 = complex(0.0,0.0);
      interp5 = complex(0.0,0.0);
      interp6 = complex(0.0,0.0);
      interp7 = complex(0.0,0.0);
      interp8 = complex(0.0,0.0);

      for (int i=0;i<N+1;i++){
      for (int j=0;j<N+1;j++){
      for (int k=0;k<N+1;k++){

	cTxyz = Patch[Dfrom].coeff[k+(N+1)*(j+i*(N+1))]*polyxyz[k+(N+1)*(j+(N+1)*i)];
//	std::cout << i << j << k <<"\n";
	// if ((Di!=(NsubD-1)) && (Dj!=(NsubD-1)) && (Dk!=(NsubD-1)) && (Patch[Dto1].level==0))
	//   interp1     -= cTxyz;
	
	// if ((Di!=(NsubD-1)) && (Dj!=(NsubD-1)) && (Dk!=0) && (Patch[Dto2].level==0))
	//   interp2     -= cTxyz;

	// if ((Di!=(NsubD-1)) && (Dj!=0) && (Dk!=(NsubD-1)) && (Patch[Dto3].level==0))
	//   interp3     -= cTxyz;

	// if ((Di!=0) && (Dj!=(NsubD-1)) && (Dk!=(NsubD-1)) && (Patch[Dto4].level==0))
	//   interp4     -= cTxyz;

	// if ((Di!=(NsubD-1)) && (Dj!=0) && (Dk!=0) && (Patch[Dto5].level==0))
	//   interp5     -= cTxyz;

	// if ((Di!=0) && (Dj!=0) && (Dk!=(NsubD-1)) && (Patch[Dto6].level==0))
	//   interp6 -= cTxyz ;  

	// if ((Di!=0) && (Dj!=(NsubD-1)) && (Dk!=0) && (Patch[Dto7].level==0))
	//   interp7  -= cTxyz;   

	// if ((Di!=0) && (Dj!=0) && (Dk!=0) && (Patch[Dto8].level==0))
	//   interp8 +=  cTxyz; 

      }}} // i,j,k loop

	if ((Di!=(NsubD-1)) && (Dj!=(NsubD-1)) && (Dk!=(NsubD-1)) && (Patch[Dto1].level==0))
	  Patch[Dto1].psi[kg+(N+1)*(jg+ig*(N+1))] = interp1;
	
	if ((Di!=(NsubD-1)) && (Dj!=(NsubD-1)) && (Dk!=0) && (Patch[Dto2].level==0))
	  Patch[Dto2].psi[N-kg+(N+1)*(jg+ig*(N+1))]  =  interp2;  

	if ((Di!=(NsubD-1)) && (Dj!=0) && (Dk!=(NsubD-1)) && (Patch[Dto3].level==0))
	  Patch[Dto3].psi[kg+(N+1)*(N-jg+ig*(N+1))]   =  interp3;  

	if ((Di!=0) && (Dj!=(NsubD-1)) && (Dk!=(NsubD-1)) && (Patch[Dto4].level==0))
	  Patch[Dto4].psi[kg+(N+1)*(jg+(N-ig)*(N+1))] =  interp4;  

	if ((Di!=(NsubD-1)) && (Dj!=0) && (Dk!=0) && (Patch[Dto5].level==0))
	  Patch[Dto5].psi[N-kg+(N+1)*(N-jg+ig*(N+1))] = interp5;  

	if ((Di!=0) && (Dj!=0) && (Dk!=(NsubD-1)) && (Patch[Dto6].level==0))
	  Patch[Dto6].psi[kg+(N+1)*(N-jg+(N-ig)*(N+1))] = interp6 ;  

	if ((Di!=0) && (Dj!=(NsubD-1)) && (Dk!=0) && (Patch[Dto7].level==0))
	  Patch[Dto7].psi[N-kg+(N+1)*(jg+(N-ig)*(N+1))]  = interp7;   

	if ((Di!=0) && (Dj!=0) && (Dk!=0) && (Patch[Dto8].level==0))
	  Patch[Dto8].psi[N-kg+(N+1)*(N-jg+(N-ig)*(N+1))] =  interp8; 
      }}} // ig,jg,kg loop

//  polyxyz.clear();
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
complex cheby_1Dinterp(double &x, vector_complex &coeff,double N,double a,double b)
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
complex cheby_2Dinterp(double &x, double &y, vector_complex &coeff,double N,
		     double ax,double bx, double ay, double by)
{
  //  int Nx = x.size();
  complex interp;
  double x_chebgrid,y_chebgrid;
  double Tx0, Ty0, Tx1, Ty1;
  vector_real polyx(N+1),polyy(N+1);

  x_chebgrid = -( 2.0 * x - ax  - bx ) / ( bx-ax);
  y_chebgrid = -( 2.0 * y - ay  - by ) / ( by-ay);

//  std::cout << x_chebgrid << " " << y_chebgrid <<"\n";
  for (int i=0;i<N+1;i++){
    polyx[i] = cheby_poly(x_chebgrid,i);
    polyy[i] = cheby_poly(y_chebgrid,i);
  }

  interp = complex(0.0,0.0) ;
  for (int i=0;i<N+1;i++)
  for (int j=0;j<N+1;j++)
      interp += coeff[i*(N+1)+j]*polyy[j]*polyx[i];

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
void fft_r2c(vector_real &x, vector_complex &y, bool forward)
{
  fftw_complex *out = reinterpret_cast<fftw_complex*>(y.data());  // standard declaration of output fttw
  fftw_plan p;

  if (forward)
  {
    p = fftw_plan_dft_r2c_1d(x.size(), x.data(), out, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    for (size_t i = 0; i < x.size(); ++i)
    {
      y[i] = y[i]/sqrt(x.size());
    }
  }

  else 
  {
    p = fftw_plan_dft_c2r_1d(x.size(), out, x.data(), FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    for (size_t i = 0; i < x.size(); ++i)
    {
      x[i] = x[i]/sqrt(x.size());
    }
  }
}  

void fft(vector_complex &x, bool inverse)
{
  vector_complex y(x.size(), complex(0.0, 0.0));
  fftw_plan p;

  fftw_complex *in = reinterpret_cast<fftw_complex*>(x.data());   // standard declaration of input fttw
  fftw_complex *out = reinterpret_cast<fftw_complex*>(y.data());  // standard declaration of output fttw

  p = fftw_plan_dft_1d(x.size(), in, out, (inverse ? FFTW_BACKWARD : FFTW_FORWARD), FFTW_ESTIMATE);

  fftw_execute(p);
  fftw_destroy_plan(p);

  for (size_t i = 0; i < x.size(); ++i)
  {
    x[i] = y[i] / sqrt(static_cast<double>(x.size()));   // fftw give unnormalized output, needs renormalization here.
  }
}
// complex GWP_GAMER_vals(double x,double Time)
// {

//   const double Gau_v0  = 32.0;                    //mean velocity [1.0]
//   const double Gau_Width =  2.0;                 // Gaussian width [0.1]
//   const double Gau_Center =  4.0;                // Gaussian center [box center]
//   const double r = x ;
 
//  double ELBDM_ETA = 1.0;
//  double Re=0,Im=0;
//  double Gau_PeriodicN = 0;

//  for (int n=0; n<Gau_PeriodicN+1; n++) {
//  for (int m=0; m<((n==0)?1:2); m++) {
 
// const double Center     = Gau_Center + n*(1-2*m)* 8.0;
//  const double dr1        = r -     Gau_v0*Time - Center;
//  const double dr2        = r - 0.5*Gau_v0*Time - Center;
//  const double Gau_Const1 = 1.0 + pow(  Time / ( ELBDM_ETA*pow(Gau_Width,2.) ), 2.0  );
//  const double Gau_Const2 = pow( pow(Gau_Width,2.)*M_PI*Gau_Const1, -0.25 )
//                             *exp(  -0.5*pow( dr1/Gau_Width, 2.0 )/Gau_Const1);
//  const double Gau_Theta1 = -0.5*acos(  pow( Gau_Const1, -0.5 )  );
//  const double Gau_Theta2 = 0.5*pow( dr1, 2.0 )*ELBDM_ETA*Time/(  pow( ELBDM_ETA*pow(Gau_Width,2.), 2.0) + pow(Time,2.)  )
//                              + Gau_v0*ELBDM_ETA*dr2;

//  Re += Gau_Const2*cos( Gau_Theta1 + Gau_Theta2 );
//  Im += Gau_Const2*sin( Gau_Theta1 + Gau_Theta2 );
//  }}
// // return complex(Gau_Const2*cos( Gau_Theta1 + Gau_Theta2 ),Gau_Const2*sin( Gau_Theta1 + Gau_Theta2));
//  return complex(Re,Im);
// }

complex GWP_vals(double x,double y,double z, double Time)
{
  
  double hbar = 1.;
  double m =1.;
  double  sigma_px = 9, sigma_py=9,sigma_pz=9;
  //  double  vx0 =  v0_init ; //2e1;
  double  vx0 =  0,  vy0 = 0, vz0 = 0e0; //2e1;
  double  x0 = x0_init;
  double  y0 = x0_init;
  double  z0 = x0_init;

  double sigma_x, phase_x, M_x;
  double sigma_y, phase_y, M_y;
  double sigma_z, phase_z, M_z;

  // periodic version parameters
  int Gau_PeriodicN = 5;
  complex sum;
  double L = 1.0;

  sigma_x = hbar/(2*sigma_px) * pow(1+4*pow(sigma_px,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
  phase_x = 1./hbar * (m*vx0+sigma_px*sigma_px/(sigma_x*sigma_x)*Time/(2*m)*(x-x0-vx0*Time))*
    (x-x0-vx0*Time) + vx0*m/(2*hbar) * vx0*Time - atan(2*sigma_px*sigma_px*Time/(hbar*m))/ 2 ;
  M_x = 1./(pow(2*M_PI,1./4)*pow(sigma_x,1./2)) * exp(-(x-x0-vx0*Time)*(x-x0-vx0*Time)/(4*sigma_x*sigma_x)) ;
  
  sigma_y = hbar/(2*sigma_py) * pow(1+4*pow(sigma_py,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
  phase_y = 1./hbar * (m*vy0+sigma_py*sigma_py/(sigma_y*sigma_y)*Time/(2*m)*(y-y0-vy0*Time))*
    (y-y0-vy0*Time) + vy0*m/(2*hbar) * vy0*Time - atan(2*sigma_py*sigma_py*Time/(hbar*m))/ 2 ;
  M_y = 1./(pow(2*M_PI,1./4)*pow(sigma_y,1./2)) * exp(-(y-y0-vy0*Time)*(y-y0-vy0*Time)/(4*sigma_y*sigma_y));
  
  sigma_z = hbar/(2*sigma_pz) * pow(1+4*pow(sigma_pz,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
  phase_z = 1./hbar * (m*vz0+sigma_pz*sigma_pz/(sigma_z*sigma_z)*Time/(2*m)*(z-z0-vz0*Time))*
    (z-z0-vz0*Time) + vz0*m/(2*hbar) * vz0*Time - atan(2*sigma_pz*sigma_pz*Time/(hbar*m))/ 2 ;
  M_z = 1./(pow(2*M_PI,1./4)*pow(sigma_z,1./2)) * exp(-(z-z0-vz0*Time)*(z-z0-vz0*Time)/(4*sigma_z*sigma_z));

  sum =  complex(M_x*M_y*M_z*cos(phase_x+phase_y+phase_z),M_x*M_y*M_z*sin(phase_x+phase_y+phase_z));

  for (int n=1; n<Gau_PeriodicN+1; n++) {
    for (int mx=-1; mx<=1; mx+=2) {
    for (int my=-1; my<=1; my+=2) {
    for (int mz=-1; mz<=1; mz+=2) {

      sigma_x = hbar/(2*sigma_px) * pow(1+4*pow(sigma_px,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
      phase_x = 1./hbar * (m*vx0+sigma_px*sigma_px/(sigma_x*sigma_x)*Time/(2*m)*(x-(x0+n*L*mx)-vx0*Time))*
	(x-(x0+n*L*mx)-vx0*Time) + vx0*m/(2*hbar) * vx0*Time - atan(2*sigma_px*sigma_px*Time/(hbar*m))/ 2 ;
      M_x = 1./(pow(2*M_PI,1./4)*pow(sigma_x,1./2)) * exp(-(x-(x0+n*L*mx)-vx0*Time)*(x-(x0+n*L*mx)-vx0*Time)/(4*sigma_x*sigma_x)) ;
      
      sigma_y = hbar/(2*sigma_py) * pow(1+4*pow(sigma_py,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
      phase_y = 1./hbar * (m*vy0+sigma_py*sigma_py/(sigma_y*sigma_y)*Time/(2*m)*(y-(y0+n*L*my)-vy0*Time))*
	(y-(y0+n*L*my)-vy0*Time) + vy0*m/(2*hbar) * vy0*Time - atan(2*sigma_py*sigma_py*Time/(hbar*m))/ 2 ;
      M_y = 1./(pow(2*M_PI,1./4)*pow(sigma_y,1./2)) * exp(-(y-(y0+n*L*my)-vy0*Time)*(y-(y0+n*L*my)-vy0*Time)/(4*sigma_y*sigma_y));

      sigma_z = hbar/(2*sigma_pz) * pow(1+4*pow(sigma_pz,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
      phase_z = 1./hbar * (m*vz0+sigma_pz*sigma_pz/(sigma_z*sigma_z)*Time/(2*m)*(z-(z0+n*L*mz)-vz0*Time))*
	(z-(z0+n*L*mz)-vz0*Time) + vz0*m/(2*hbar) * vz0*Time - atan(2*sigma_pz*sigma_pz*Time/(hbar*m))/ 2 ;
      M_z = 1./(pow(2*M_PI,1./4)*pow(sigma_z,1./2)) * exp(-(z-(z0+n*L*mz)-vz0*Time)*(z-(z0+n*L*mz)-vz0*Time)/(4*sigma_z*sigma_z));

      sum += complex(M_x*M_y*M_z*cos(phase_x+phase_y+phase_z),M_x*M_y*M_z*sin(phase_x+phase_y+phase_z));
    }}}}

  return sum;
//  return complex(M_x*M_y*M_z*cos(phase_x+phase_y+phase_z),M_x*M_y*M_z*sin(phase_x+phase_y+phase_z));

}

void initial_GWP(vector_real &x, 
		 vector_real &y,
		 vector_real &z,
		 vector_complex &psi)
{
  
  unsigned long int N = x.size();
  double Time = 0;

  for (size_t i = 0; i < N; ++i){
  for (size_t j = 0; j < N; ++j){
  for (size_t k = 0; k < N; ++k){
    psi[k+N*(i*N+j)] = GWP_vals(x[i], y[j],z[k], Time); //imag
//    psi[k+N*(j+i*N)] = GWP_2Dvals(x[i], z[k], Time); //imag
//    psi[k+N*(j+i*N)] = GWP_2Dvals(y[j], z[k], Time); //imag
    //    std::cout << x[i]<< " " << std::abs(psi[i]) <<  "\n"; 
    }}}
  //  write_data(par, opr, 0);

}
void Print(Params &par,
	   vector_real &x,
	   vector_real &y,
	   vector_real &z,
	   vector_complex  &psi,
	   int tn){

  double size = x.size()-1;

  for (int i=0; i< size+1; i++){
  for (int j=0; j< size+1; j++){
  for (int k=0; k< size+1; k++){
    std::cout << x[i] << " " 
  	      << std::abs(psi[i]) << " " 
  	      << std::abs(GWP_vals(x[i],y[j],z[k], tn*par.dt))<< " "
  	      << psi[i*(size+1)+j].real() << " " 
  	      << GWP_vals(x[i],y[j],z[k], tn*par.dt).real()<< " " 
  	      << psi[i*(size+1)+j].imag() << " " 
  	      << GWP_vals(x[i],y[j],z[k], tn*par.dt).imag()<< "\n"; 
  }}}
}
void Print_XZPatch(Params &par,
		 vector_patch &Patch,
		 unsigned long int NsubD,
		 int tn,
		 int D){

  int size;
  int Ng = 8;

  std::stringstream filename_fdm;
  filename_fdm << "output/XZP" << D<<"_"<< tn << ".dat";
  std::ofstream ffdm = std::ofstream(filename_fdm.str());

  if (ffdm)
  {
    std::stringstream data_fdm;
    int j = Patch[D].N-2;	
     for (int t=0; t<Patch[D].N+1; t++)
	for (int tt=0; tt<Patch[D].N+1; tt++){

	  data_fdm << std::setprecision(25)
		   << Patch[D].x[t] << " "
		    << Patch[D].z[tt] << " "
		    << t << " "
		    << tt << " "
		    << Patch[D].psi[tt+(Patch[D].N+1)*(j+t*(Patch[D].N+1))].real() << " "
		    // << cheby_2Dinterp(Patch[D].x[t],
		    // 		    Patch[D].y[tt],
		    // 		    Patch[D].coeff,
		    // 		    Patch[D].N,
		    // 		    Patch[D].x[0],
		    // 		    Patch[D].x.back(),
		    // 		    Patch[D].y[0],
		    // 		    Patch[D].y.back()).real() << " "
		   //<< Patch[D].coeff[t*(Patch[D].N+1)+tt].real() << " "
  		   << GWP_vals(Patch[D].x[t],Patch[D].y[j], Patch[D].z[tt],tn*par.dt).real()<< "\n"; }

    ffdm.write(data_fdm.str().c_str(), data_fdm.str().length());
    ffdm.close();
  }
}
void Print_XYPatch(Params &par,
		 vector_patch &Patch,
		 unsigned long int NsubD,
		 int tn,
		 int D){

  int size;
  int Ng = 8;

  std::stringstream filename_fdm;
  filename_fdm << "output/XYP" << D<<"_"<< tn << ".dat";
  std::ofstream ffdm = std::ofstream(filename_fdm.str());

  if (ffdm)
  {
    std::stringstream data_fdm;
    int k = Patch[D].N/2;	
     for (int t=0; t<Patch[D].N+1; t++)
	for (int tt=0; tt<Patch[D].N+1; tt++){

	  data_fdm << std::setprecision(25)
		   << Patch[D].x[t] << " "
		    << Patch[D].y[tt] << " "
		    << t << " "
		    << tt << " "
		    << Patch[D].psi[k+(Patch[D].N+1)*(tt+t*(Patch[D].N+1))].real() << " "
		    // << cheby_2Dinterp(Patch[D].x[t],
		    // 		    Patch[D].y[tt],
		    // 		    Patch[D].coeff,
		    // 		    Patch[D].N,
		    // 		    Patch[D].x[0],
		    // 		    Patch[D].x.back(),
		    // 		    Patch[D].y[0],
		    // 		    Patch[D].y.back()).real() << " "
		   //<< Patch[D].coeff[t*(Patch[D].N+1)+tt].real() << " "
  		   << GWP_vals(Patch[D].x[t],Patch[D].y[tt], Patch[D].z[k],tn*par.dt).real()<< "\n"; }

//     std::cout << data_fdm.str();
    ffdm.write(data_fdm.str().c_str(), data_fdm.str().length());
    ffdm.close();
  }
  //  }
}
void Print_YZPatch(Params &par,
		 vector_patch &Patch,
		 unsigned long int NsubD,
		 int tn,
		 int D){

  int size;
  int Ng = 8;

  std::stringstream filename_fdm;
  filename_fdm << "output/YZP" << D<<"_"<< tn << ".dat";
  std::ofstream ffdm = std::ofstream(filename_fdm.str());

  if (ffdm)
  {
    std::stringstream data_fdm;
//    int i = Patch[D].N/2;	
    int i = Patch[D].N-2;	
     for (int t=0; t<Patch[D].N+1; t++)
	for (int tt=0; tt<Patch[D].N+1; tt++){

	  data_fdm << std::setprecision(25)
		   << Patch[D].y[t] << " "
		    << Patch[D].z[tt] << " "
		    << t << " "
		    << tt << " "
		    << Patch[D].psi[tt+(Patch[D].N+1)*(t+i*(Patch[D].N+1))].real() << " "
		    // << cheby_2Dinterp(Patch[D].x[t],
		    // 		    Patch[D].y[tt],
		    // 		    Patch[D].coeff,
		    // 		    Patch[D].N,
		    // 		    Patch[D].x[0],
		    // 		    Patch[D].x.back(),
		    // 		    Patch[D].y[0],
		    // 		    Patch[D].y.back()).real() << " "
		   //<< Patch[D].coeff[t*(Patch[D].N+1)+tt].real() << " "
  		   << GWP_vals(Patch[D].x[i],Patch[D].y[t], Patch[D].z[tt],tn*par.dt).real()<< "\n"; }

//     std::cout << data_fdm.str();
    ffdm.write(data_fdm.str().c_str(), data_fdm.str().length());
    ffdm.close();
  }
  //  }
}

void Print_XYPlane(Params &par,
		vector_patch &P,
		unsigned long int NsubD,
		int tn){

  int size,D,k;
  int Ng = 8;

  std::stringstream filename_fdm;
  filename_fdm << "output/XY" << tn << ".dat";
  std::ofstream ffdm = std::ofstream(filename_fdm.str());

  if (ffdm)
  {
    std::stringstream data_fdm;

    for (int Di=0;Di<NsubD;Di++){
      for (int Dj=0;Dj<NsubD;Dj++){
	if (NsubD !=1){
	  D = NsubD/2 -1 + NsubD*(Di*NsubD+Dj);
	  k = P[D].N-1;
	}
	else {
	  D = NsubD*(Di*NsubD+Dj);
	  k = P[D].N/2;	
	}
//	for (int i = Ng; i < P[D].N+1-Ng; ++i){
//	for (int j = Ng; j < P[D].N+1-Ng; ++j){
	for (int i = 0; i < P[D].N+1; ++i){
	for (int j = 0; j < P[D].N+1; ++j){
   	  data_fdm << i + Di*(P[D].N+1) << " "
		   << j + Dj*(P[D].N+1) << " " 
  		   << P[D].x[i] << " " 
  		   << P[D].y[j] << " " 
  		   << std::abs(P[D].psi[k+(P[D].N+1)*(i*(P[D].N+1)+j)]) << " " 
  		   << std::abs(GWP_vals(P[D].x[i],P[D].y[j],P[D].z[k], tn*par.dt))<< " "
  		   << P[D].psi[k+(P[D].N+1)*(i*(P[D].N+1)+j)].real() << " " 
  		   << GWP_vals(P[D].x[i],P[D].y[j],P[D].z[k], tn*par.dt).real()<< " " 
//  		   << GWP_2Dvals(P[D].x[i],P[D].y[j],P[D].z[k], tn*par.dt).real()<< " " 
  		   << P[D].psi[k+(P[D].N+1)*(i*(P[D].N+1)+j)].imag() << " " 
  		   << GWP_vals(P[D].x[i],P[D].y[j],P[D].z[k], tn*par.dt).imag()<< "\n"; 
//  		   << GWP_2Dvals(P[D].x[i], P[D].y[j],tn*par.dt).imag()<< "\n";
  	}} // for i,j
      }} // for Di, Dj 
    ffdm.write(data_fdm.str().c_str(), data_fdm.str().length());
    ffdm.close();
  }
  //  }
  //std::exit(0);
}
void Print_XZPlane(Params &par,
		vector_patch &P,
		unsigned long int NsubD,
		int tn)
{
  int size,D,j;
  int Ng = 8;
  std::stringstream filename_fdm;
  filename_fdm << "output/XZ" << tn << ".dat";
  std::ofstream ffdm = std::ofstream(filename_fdm.str());

  if (ffdm)
  {
    std::stringstream data_fdm;
    for (int Di=0;Di<NsubD;Di++){
      for (int Dj=0;Dj<NsubD;Dj++){
	if (NsubD !=1){
	  D = NsubD/2 -1 + NsubD*(Di*NsubD+Dj);
	  j = P[D].N-1;
	}
	else {
	  D = NsubD*(Di*NsubD+Dj);
	  j = P[D].N/2;	
	}
	for (int i = Ng; i < P[D].N+1-Ng; ++i){
	for (int k = Ng; k < P[D].N+1-Ng; ++k){
   	  data_fdm << i + Di*(P[D].N+1) << " "
		   << k + Dj*(P[D].N+1) << " " 
  		   << P[D].x[i] << " " 
  		   << P[D].z[k] << " "
  		   << std::abs(P[D].psi[k+(P[D].N+1)*(i*(P[D].N+1)+j)]) << " " 
  		   << std::abs(GWP_vals(P[D].x[i],P[D].y[j],P[D].z[k], tn*par.dt))<< " "
  		   << P[D].psi[k+(P[D].N+1)*(i*(P[D].N+1)+j)].real() << " " 
  		   << GWP_vals(P[D].x[i], P[D].y[j],P[D].z[k],tn*par.dt).real()<< " "
//  		   << GWP_2Dvals(P[D].x[i],P[D].y[j],P[D].z[k], tn*par.dt).real()<< " " 
  		   << P[D].psi[k+(P[D].N+1)*(i*(P[D].N+1)+j)].imag() << " " 
  		   << GWP_vals(P[D].x[i], P[D].y[j],P[D].z[k],tn*par.dt).imag()<< "\n";
//  		   << GWP_2Dvals(P[D].x[i], P[D].z[k],tn*par.dt).imag()<< "\n";
  	}} // for i,k
      }} // for Di, Dj
    ffdm.write(data_fdm.str().c_str(), data_fdm.str().length());
    ffdm.close();
  }
}
void Print_YZPlane(Params &par,
		vector_patch &P,
		unsigned long int NsubD,
		int tn){

  int size,D,i;
  int Ng=8;
  std::stringstream filename_fdm;
  filename_fdm << "output/YZ" << tn << ".dat";
  std::ofstream ffdm = std::ofstream(filename_fdm.str());

  if (ffdm)
  {
    std::stringstream data_fdm;

    for (int Di=0;Di<NsubD;Di++){
      for (int Dj=0;Dj<NsubD;Dj++){
	if (NsubD !=1){
	  D = NsubD/2 -1 + NsubD*(Di*NsubD+Dj);
	  i = P[D].N-1;
	}
	else {
	  D = NsubD*(Di*NsubD+Dj);
	  i = P[D].N/2;	
	}
	for (int j = Ng; j < P[D].N+1-Ng; ++j){
	for (int k = Ng; k < P[D].N+1-Ng; ++k){
   	  data_fdm << j + Di*(P[D].N+1) << " "
		   << k + Dj*(P[D].N+1) << " " 
  		   << P[D].y[j] << " " 
  		   << P[D].z[k] << " " 
  		   << std::abs(P[D].psi[k+(P[D].N+1)*(i*(P[D].N+1)+j)]) << " " 
  		   << std::abs(GWP_vals(P[D].x[i],P[D].y[j],P[D].z[k], tn*par.dt))<< " "
  		   << P[D].psi[k+(P[D].N+1)*(i*(P[D].N+1)+j)].real() << " " 
  		   << GWP_vals(P[D].x[i], P[D].y[j],P[D].z[k],tn*par.dt).real()<< " "
//  		   << GWP_2Dvals(P[D].x[i],P[D].y[j],P[D].z[k], tn*par.dt).real()<< " " 
  		   << P[D].psi[k+(P[D].N+1)*(i*(P[D].N+1)+j)].imag() << " " 
  		   << GWP_vals(P[D].x[i], P[D].y[j],P[D].z[k],tn*par.dt).imag()<< "\n";
//  		   << GWP_2Dvals(P[D].y[j], P[D].z[k],tn*par.dt).imag()<< "\n";
  	}
      }
      }
    }
    ffdm.write(data_fdm.str().c_str(), data_fdm.str().length());
    ffdm.close();
  }
}
void Print_XYZ(Params &par,
		vector_patch &P,
		unsigned long int NsubD,
		int tn){

  int size,D;
  int Ng = P[0].N_ghost;

  std::stringstream filename_fdm;
  filename_fdm << "output6/XYZ" << tn << ".dat";
  std::ofstream ffdm = std::ofstream(filename_fdm.str());

  if (ffdm)
  {
    std::stringstream data_fdm;

    for (int Di=0;Di<NsubD;Di++){
      for (int Dj=0;Dj<NsubD;Dj++){
	for (int Dk=0;Dk<NsubD;Dk++){
	  D = Dk + NsubD*(Di*NsubD+Dj);

	for (int i = Ng; i < P[D].N+1-Ng; ++i){
	for (int j = Ng; j < P[D].N+1-Ng; ++j){
	for (int k = Ng; k < P[D].N+1-Ng; ++k){
	  data_fdm << std::setprecision(25)
		   << i + Di*(P[D].N+1) << " "
		   << j + Dj*(P[D].N+1) << " " 
  		   << P[D].x[i] << " " 
  		   << P[D].y[j] << " " 
  		   << std::abs(P[D].psi[k+(P[D].N+1)*(i*(P[D].N+1)+j)]) << " " 
  		   << std::abs(GWP_vals(P[D].x[i],P[D].y[j],P[D].z[k], tn*par.dt))<< " "
  		   << P[D].psi[k+(P[D].N+1)*(i*(P[D].N+1)+j)].real() << " " 
  		   << GWP_vals(P[D].x[i],P[D].y[j],P[D].z[k], tn*par.dt).real()<< " " 
//  		   << GWP_2Dvals(P[D].x[i],P[D].y[j],P[D].z[k], tn*par.dt).real()<< " " 
  		   << P[D].psi[k+(P[D].N+1)*(i*(P[D].N+1)+j)].imag() << " " 
  		   << GWP_vals(P[D].x[i],P[D].y[j],P[D].z[k], tn*par.dt).imag()<< "\n"; 
//  		   << GWP_2Dvals(P[D].x[i], P[D].y[j],tn*par.dt).imag()<< "\n";
  	}}}
    }}}
    ffdm.write(data_fdm.str().c_str(), data_fdm.str().length());
    ffdm.close();
  }
  //  }
  //std::exit(0);
}

void Error(Params &par,
	   vector_real &x,
	   vector_real &y,
	   vector_real &z,
	   vector_complex  &psi,
	   double &Error){
  
  int Ng = 8;
  double size = x.size();
  for (size_t i = Ng; i < size-Ng; ++i){
  for (size_t j = Ng; j < size-Ng; ++j){
  for (size_t k = Ng; k < size-Ng; ++k){
    // Error += pow(std::abs(psi[k+size*(j+i*(size))])
    // 		 - std::abs(GWP_vals(x[i],y[j],z[k], par.timesteps*par.dt)),2);
   Error += pow((psi[k+size*(j+i*(size))].real() * psi[k+size*(j+i*(size))].real() + psi[k+size*(j+i*(size))].imag() * psi[k+size*(j+i*(size))].imag())
   		 - (GWP_vals(x[i],y[j],z[k], par.timesteps*par.dt).real()*GWP_vals(x[i],y[j],z[k], par.timesteps*par.dt).real() + 
		    GWP_vals(x[i],y[j],z[k], par.timesteps*par.dt).imag()*GWP_vals(x[i],y[j],z[k], par.timesteps*par.dt).imag()),2);
  }}}
  //  error = pow(error,1./2)/opr.size;  
  //  std::cout << "error = " << error << "\n";    
}
void MaxError(Params &par,
	   vector_real &x,
	   vector_real &y,
	   vector_real &z,
	   vector_complex  &psi,
	   double &Error){
  double err;

  double size = x.size();
  for (size_t i = 1; i < size-1; ++i){
  for (size_t j = 1; j < size-1; ++j){
  for (size_t k = 1; k < size-1; ++k){
    err = pow((psi[k+size*(j+i*(size))].real() * psi[k+size*(j+i*(size))].real() + psi[k+size*(j+i*(size))].imag() * psi[k+size*(j+i*(size))].imag())
    		 - (GWP_vals(x[i],y[j],z[k], par.timesteps*par.dt).real()*GWP_vals(x[i],y[j],z[k], par.timesteps*par.dt).real()+GWP_vals(x[i],y[j],z[k], par.timesteps*par.dt).imag()*GWP_vals(x[i],y[j],z[k], par.timesteps*par.dt).imag()),2);
    if (std::abs(err) >std::abs(Error)) 
      Error = std::abs(err);
  }}}
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

// void Refinement(Patches &P){
  
//   // Perform Refinement: double cheby grid size

//   // for (int i=0; i< P.N+1; i++){
//   //   std::cout << P.x[i] << " " 
//   // 	      << P.psi[i].real() << " " 
//   // 	      << P.psi[i].imag() <<"\n";
//   // }
//   P.coeff = cheby_1Dcoeff(P.x,P.psi,P.N);
//   P.x  = cheby_exgrid(2*P.N,P.x[0],P.x[P.N]);
//   P.psi = cheby_interp(P.x,P.coeff,P.N,P.x[0],P.x[2*P.N]); 
//   cheby_DD_matrix(P.D, P.x, 2*P.N);
//   P.level = 1;

//   P.N *= 2.0; 
//   P.tmp.resize(P.N);
//   P.k1.resize(P.N);
//   P.k2.resize(P.N);
//   P.k3.resize(P.N);
//   P.k4.resize(P.N);
//   // std::cout <<"\n";

//   // for (int i=0; i< P.N+1; i++){
//   //   std::cout << P.x[i] << " " 
//   // 	      << P.psi[i].real() << " " 
//   // 	      << P.psi[i].imag() <<"\n";
//   // }

//   //   std::exit(0);
// }

// void Derefinement(Patches &P){
  
//   // Perform De-Refinement: half cheby grid size

//   // for (int i=0; i< P.N+1; i++){
//   //   std::cout << P.x[i] << " " 
//   // 	      << P.psi[i].real() << " " 
//   // 	      << P.psi[i].imag() <<"\n";
//   // }
//   double half = 0.5 ; 
//   P.coeff = cheby_1Dcoeff(P.x,P.psi,P.N);
//   P.x  = cheby_exgrid(half*P.N,P.x[0],P.x[P.N]);
//   P.psi = cheby_interp(P.x,P.coeff,P.N,P.x[0],P.x[half*P.N]); 
//   cheby_DD_matrix(P.D, P.x, half*P.N);
//   P.level = 0;

//   P.N *= half; 
//   P.tmp.resize(P.N);
//   P.k1.resize(P.N);
//   P.k2.resize(P.N);
//   P.k3.resize(P.N);
//   P.k4.resize(P.N);
//   //  std::cout <<"\n";

//   // for (int i=0; i< P.N+1; i++){
//   //   std::cout << P.x[i] << " " 
//   // 	      << P.psi[i].real() << " " 
//   // 	      << P.psi[i].imag() <<"\n";
//   // }

//   //   std::exit(0);
// }
void Assign_1D2D3DBC(vector_patch &Patch,
		     signed long int NsubD,
		     Params &par,
		     double tn){
 

  fftw_plan plan2D_R, plan2D_I;
  int N , D, Dfrom,Dto1,Dto2,Dto3,Dto4, Dto5,Dto6,Dto7,Dto8;
  vector_complex psi1D(Patch[0].N+1),coeff1D(Patch[0].N+1);
  vector_complex psi2D((Patch[0].N+1) * (Patch[0].N+1)),coeff2D((Patch[0].N+1) * (Patch[0].N+1));
  int ghost_pts = Patch[0].N_ghost;
  //    std::cout << tn <<"\n";
  //  unsigned long int N = Patch[0].N;
//  std::cout << "start 1D2D3D\n";

  for (int Di=0;Di<NsubD;Di++) {
  for (int Dj=0;Dj<NsubD;Dj++) {
  for (int Dk=0;Dk<NsubD;Dk++) {
    D = Dk + NsubD*(Dj+Di*NsubD);
    N = Patch[D].N;

    // X0 face
    if (Di==0){
      for (int j = 0; j<N+1; j++){
      for (int k = 0; k<N+1; k++){
      for (int ig=0; ig< ghost_pts; ig++){

	Patch[D].psi[k+(N+1)*(j+(N+1)*ig)] = GWP_vals(Patch[D].x[ig],Patch[D].y[j],Patch[D].z[k],(tn-1.)*par.dt);
      }}}
    }

//    XN face
    if (Di==(NsubD-1)){
      for (int j = 0; j<N+1; j++){
      for (int k = 0; k<N+1; k++){
      for (int ig=0; ig< ghost_pts; ig++){
	Patch[D].psi[k+(N+1)*(j+(N-ig)*(N+1))] = GWP_vals(Patch[D].x[N-ig],Patch[D].y[j],Patch[D].z[k],(tn-1.)*par.dt); 
	//	std::cout << "b"<<j<<"="<< Patch[D].b[j] <<"\n";
      }}}
    }

//    Y0 face
    if (Dj==0){
      for (int i = 0; i<N+1; i++){
      for (int k = 0; k<N+1; k++){
      for (int jg=0; jg< ghost_pts; jg++){

	Patch[D].psi[k+(N+1)*(jg+i*(N+1))] = GWP_vals(Patch[D].x[i],Patch[D].y[jg],Patch[D].z[k],(tn-1.)*par.dt); 
	//	std::cout << "c"<<i<<"="<< Patch[D].c[i] <<"\n";
      }}}
    }

//    YN face
    if (Dj==(NsubD-1)){
      for (int i = 0; i<N+1; i++){
      for (int k = 0; k<N+1; k++){
      for (int jg=0; jg< ghost_pts; jg++){
	Patch[D].psi[k+(N+1)*(i*(N+1)+N-jg)] = GWP_vals(Patch[D].x[i],Patch[D].y[N-jg],Patch[D].z[k],(tn-1.)*par.dt); 
	//	std::cout << "d"<<i<<"="<< Patch[D].d[i] <<"\n";
      }}}
    }

//    Z0 face
    if (Dk==0){
      for (int i = 0; i<N+1; i++){
      for (int j = 0; j<N+1; j++){
      for (int kg=0; kg< ghost_pts; kg++){

	Patch[D].psi[kg+(N+1)*(j+i*(N+1))] = GWP_vals(Patch[D].x[i],Patch[D].y[j],Patch[D].z[kg],(tn-1.)*par.dt); 
	//	std::cout << "c"<<i<<"="<< Patch[D].c[i] <<"\n";
      }}}
    }

//    ZN face
    if (Dk==(NsubD-1)){
      for (int i = 0; i<N+1; i++){
      for (int j = 0; j<N+1; j++){
      for (int kg=0; kg< ghost_pts; kg++){
	Patch[D].psi[N-kg+(N+1)*(j+i*(N+1))] = GWP_vals(Patch[D].x[i],Patch[D].y[j],Patch[D].z[N-kg],(tn-1.)*par.dt); 
	//	std::cout << "d"<<i<<"="<< Patch[D].d[i] <<"\n";
      }}}
    }


//    std::chrono::steady_clock::time_point start1D = std::chrono::steady_clock::now();
    // X +- Face 
    Dfrom = Dk + NsubD*(Dj+Di*NsubD);
    Dto1  = Dk + NsubD*(Dj+(Di-1)*NsubD);
    Dto2  = Dk + NsubD*(Dj+(Di+1)*NsubD);

    for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){	
    for (int kg=ghost_pts; kg<N+1-ghost_pts; kg++){

      for (int i=0; i< N+1; i++)
	psi1D[i] = Patch[Dfrom].psi[kg+(N+1)*(jg+i*(N+1))];
      
      coeff1D = cheby_1Dcoeff(Patch[Dfrom].x, 
			      psi1D, 
			      Patch[Dfrom].N);

      for (int ig=0; ig< ghost_pts; ig++){
	if ((Di!=0) && (Patch[Dto1].level==0))
	  Patch[Dto1].psi[kg+(N+1)*(jg+(N-ig)*(N+1))] = cheby_1Dinterp(Patch[Dto1].x[N-ig],
									 coeff1D,
									 Patch[Dfrom].N,
									 Patch[Dfrom].x[0],
									 Patch[Dfrom].x.back());
	
	if ((Di!=(NsubD-1)) && (Patch[Dto2].level==0))
	  Patch[Dto2].psi[kg+(N+1)*(jg+ig*(N+1))] = cheby_1Dinterp(Patch[Dto2].x[ig],
								   coeff1D,
								   Patch[Dfrom].N,
								   Patch[Dfrom].x[0],
								   Patch[Dfrom].x.back());
	}
    }} // kg, jg loop

    // // Y +- Face 
    Dfrom = Dk + NsubD*(Dj+Di*NsubD);
    Dto1  = Dk + NsubD*(Dj-1+Di*NsubD);
    Dto2  = Dk + NsubD*(Dj+1+Di*NsubD);

    for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){	
    for (int kg=ghost_pts; kg<N+1-ghost_pts; kg++){

      for (int j=0; j< N+1; j++)
	psi1D[j] = Patch[Dfrom].psi[kg+(N+1)*(j+ig*(N+1))];
      
      coeff1D = cheby_1Dcoeff(Patch[Dfrom].y, 
			      psi1D, 
			      Patch[Dfrom].N);

      for (int jg=0; jg< ghost_pts; jg++){
	if ((Dj!=0) && (Patch[Dto1].level==0))
	  Patch[Dto1].psi[kg+(N+1)*(N-jg+ig*(N+1))] = cheby_1Dinterp(Patch[Dto1].y[N-jg],
								     coeff1D,
								     Patch[Dfrom].N,
								     Patch[Dfrom].y[0],
								     Patch[Dfrom].y.back());
	
	if ((Dj!=(NsubD-1)) && (Patch[Dto2].level==0))
	  Patch[Dto2].psi[kg+(N+1)*(jg+ig*(N+1))] = cheby_1Dinterp(Patch[Dto2].y[jg],
								   coeff1D,
								   Patch[Dfrom].N,
								   Patch[Dfrom].y[0],
								   Patch[Dfrom].y.back());
	}
    }} // kg, ig loop

    // // Z +- Face 
    Dfrom = Dk + NsubD*(Dj+Di*NsubD);
    Dto1  = Dk-1 + NsubD*(Dj+Di*NsubD);
    Dto2  = Dk+1 + NsubD*(Dj+Di*NsubD);

    for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
    for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){	

      for (int k=0; k< N+1; k++)
	psi1D[k] = Patch[Dfrom].psi[k+(N+1)*(jg+ig*(N+1))];
      
      coeff1D = cheby_1Dcoeff(Patch[Dfrom].z, 
				psi1D, 
				Patch[Dfrom].N);

      for (int kg=0; kg< ghost_pts; kg++){

	if ((Dk!=0) && (Patch[Dto1].level==0))
	  Patch[Dto1].psi[N-kg+(N+1)*(jg+ig*(N+1))] = cheby_1Dinterp(Patch[Dto1].z[N-kg],
								     coeff1D,
								     Patch[Dfrom].N,
								     Patch[Dfrom].z[0],
								     Patch[Dfrom].z.back());
	
	if ((Dk!=(NsubD-1)) && (Patch[Dto2].level==0))
	  Patch[Dto2].psi[kg+(N+1)*(jg+ig*(N+1))] = cheby_1Dinterp(Patch[Dto2].z[kg],
								   coeff1D,
								   Patch[Dfrom].N,
								   Patch[Dfrom].z[0],
								   Patch[Dfrom].z.back());
	}
    }} // ig, jg loop

    // std::chrono::steady_clock::time_point end1D = std::chrono::steady_clock::now();
    // std::cout << "1D interpolations take "  << std::chrono::duration_cast<std::chrono::microseconds>(end1D - start1D).count()*1e-6<< " seconds"<< std::endl;

    // std::chrono::steady_clock::time_point start2D = std::chrono::steady_clock::now();
    // XY slices
    Dfrom = Dk + NsubD*(Dj+Di*NsubD);
    Dto1  = Dk + NsubD*(Dj-1+(Di-1)*NsubD);
    Dto2  = Dk + NsubD*(Dj-1+(Di+1)*NsubD);
    Dto3  = Dk + NsubD*(Dj+1+(Di-1)*NsubD);
    Dto4  = Dk + NsubD*(Dj+1+(Di+1)*NsubD);
    
    for (int kg=ghost_pts; kg<N+1-ghost_pts; kg++){

      for (int i=0; i< N+1; i++)
      for (int j=0; j< N+1; j++)
	psi2D[j+(N+1)*i] = Patch[Dfrom].psi[kg+(N+1)*(j+i*(N+1))];
      
      cheby_2Dcoeff(psi2D, 
		    Patch[Dfrom].N,
		    coeff2D,
		    Patch[Dfrom]);
      
      for (int ig=0; ig < ghost_pts; ig++){
      for (int jg=0; jg < ghost_pts; jg++){	

	if ((Di!=0) && (Dj!=0) && (Patch[Dto1].level==0))
	  Patch[Dto1].psi[kg+(N+1)*(N-jg+(N-ig)*(N+1))] =  cheby_2Dinterp(Patch[Dto1].x[N-ig],
									  Patch[Dto1].y[N-jg],
									  coeff2D,
									  Patch[Dfrom].N,
									  Patch[Dfrom].x[0],
									  Patch[Dfrom].x.back(),
									  Patch[Dfrom].y[0],
									  Patch[Dfrom].y.back());
	
	if ((Di!=(NsubD-1)) && (Dj!=0) && (Patch[Dto2].level==0))
	  Patch[Dto2].psi[kg+(N+1)*(N-jg+ig*(N+1))] =  cheby_2Dinterp(Patch[Dto2].x[ig],
								      Patch[Dto2].y[N-jg],
								      coeff2D,
								      Patch[Dfrom].N,
								      Patch[Dfrom].x[0],
								      Patch[Dfrom].x.back(),
								      Patch[Dfrom].y[0],
								      Patch[Dfrom].y.back());

	if ((Di!=0) && (Dj!=(NsubD-1)) && (Patch[Dto3].level==0))
	  Patch[Dto3].psi[kg+(N+1)*(jg+(N-ig)*(N+1))] =  cheby_2Dinterp(Patch[Dto3].x[N-ig],
									Patch[Dto3].y[jg],
									coeff2D,
									Patch[Dfrom].N,
									Patch[Dfrom].x[0],
									Patch[Dfrom].x.back(),
									Patch[Dfrom].y[0],
									Patch[Dfrom].y.back());
	
	if ((Di!=(NsubD-1)) && (Dj!=(NsubD-1)) && (Patch[Dto4].level==0))
	  Patch[Dto4].psi[kg+(N+1)*(jg+ig*(N+1))] =  cheby_2Dinterp(Patch[Dto4].x[ig],
								    Patch[Dto4].y[jg],
								    coeff2D,
								    Patch[Dfrom].N,
								    Patch[Dfrom].x[0],
								    Patch[Dfrom].x.back(),
								    Patch[Dfrom].y[0],
								    Patch[Dfrom].y.back());
	}} // jg, kg loop
    } // ig loop

    // // XZ slices
    Dfrom = Dk + NsubD*(Dj+Di*NsubD);
    Dto1  = Dk-1 + NsubD*(Dj+(Di-1)*NsubD);
    Dto2  = Dk-1 + NsubD*(Dj+(Di+1)*NsubD);
    Dto3  = Dk+1 + NsubD*(Dj+(Di-1)*NsubD);
    Dto4  = Dk+1 + NsubD*(Dj+(Di+1)*NsubD);

    for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){

      for (int i=0; i< N+1; i++)
      for (int k=0; k< N+1; k++)
	psi2D[k+(N+1)*i] = Patch[Dfrom].psi[k+(N+1)*(jg+i*(N+1))];
      
      cheby_2Dcoeff(psi2D, 
		    Patch[Dfrom].N,
		    coeff2D,
		    Patch[Dfrom]);
      
      for (int ig=0; ig < ghost_pts; ig++){	
      for (int kg=0; kg < ghost_pts; kg++){

	if ((Di!=0) && (Dk!=0) && (Patch[Dto1].level==0))
	  Patch[Dto1].psi[N-kg+(N+1)*(jg+(N-ig)*(N+1))] =  cheby_2Dinterp(Patch[Dto1].x[N-ig],
									  Patch[Dto1].z[N-kg],
									  coeff2D,
									  Patch[Dfrom].N,
									  Patch[Dfrom].x[0],
									  Patch[Dfrom].x.back(),
									  Patch[Dfrom].z[0],
									  Patch[Dfrom].z.back());
	
	if ((Di!=(NsubD-1)) && (Dk!=0) && (Patch[Dto2].level==0))
	  Patch[Dto2].psi[N-kg+(N+1)*(jg+ig*(N+1))] =  cheby_2Dinterp(Patch[Dto2].x[ig],
								      Patch[Dto2].z[N-kg],
								      coeff2D,
								      Patch[Dfrom].N,
								      Patch[Dfrom].x[0],
								      Patch[Dfrom].x.back(),
								      Patch[Dfrom].z[0],
								      Patch[Dfrom].z.back());

	if ((Di!=0) && (Dk!=(NsubD-1)) && (Patch[Dto3].level==0))
	  Patch[Dto3].psi[kg+(N+1)*(jg+(N-ig)*(N+1))] =  cheby_2Dinterp(Patch[Dto3].x[N-ig],
								      Patch[Dto3].z[kg],
								      coeff2D,
								      Patch[Dfrom].N,
								      Patch[Dfrom].x[0],
								      Patch[Dfrom].x.back(),
								      Patch[Dfrom].z[0],
								      Patch[Dfrom].z.back());
	
	if ((Di!=(NsubD-1)) && (Dk!=(NsubD-1)) && (Patch[Dto4].level==0))
	  Patch[Dto4].psi[kg+(N+1)*(jg+ig*(N+1))] =  cheby_2Dinterp(Patch[Dto4].x[ig],
								    Patch[Dto4].z[kg],
								    coeff2D,
								    Patch[Dfrom].N,
								    Patch[Dfrom].x[0],
								    Patch[Dfrom].x.back(),
								    Patch[Dfrom].z[0],
								    Patch[Dfrom].z.back());
	}} // jg, kg loop
    } // ig loop

    // // YZ slices
    Dfrom = Dk + NsubD*(Dj+Di*NsubD);
    Dto1  = Dk-1 + NsubD*(Dj-1+Di*NsubD);
    Dto2  = Dk-1 + NsubD*(Dj+1+Di*NsubD);
    Dto3  = Dk+1 + NsubD*(Dj-1+Di*NsubD);
    Dto4  = Dk+1 + NsubD*(Dj+1+Di*NsubD);

    for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){

      for (int j=0; j< N+1; j++)
      for (int k=0; k< N+1; k++)
	psi2D[k+(N+1)*j] = Patch[Dfrom].psi[k+(N+1)*(j+ig*(N+1))];
      
      cheby_2Dcoeff(psi2D, 
		    Patch[Dfrom].N,
		    coeff2D,
		    Patch[Dfrom]);
      
      for (int jg=0; jg < ghost_pts; jg++){	
      for (int kg=0; kg < ghost_pts; kg++){

	if ((Dj!=0) && (Dk!=0) && (Patch[Dto1].level==0))
	  Patch[Dto1].psi[N-kg+(N+1)*(N-jg+ig*(N+1))] =  cheby_2Dinterp(Patch[Dto1].y[N-jg],
									Patch[Dto1].z[N-kg],
									coeff2D,
									Patch[Dfrom].N,
									Patch[Dfrom].y[0],
									Patch[Dfrom].y.back(),
									Patch[Dfrom].z[0],
									Patch[Dfrom].z.back());
	
	if ((Dj!=(NsubD-1)) && (Dk!=0) && (Patch[Dto2].level==0))
	  Patch[Dto2].psi[N-kg+(N+1)*(jg+ig*(N+1))] =  cheby_2Dinterp(Patch[Dto2].y[jg],
								      Patch[Dto2].z[N-kg],
								      coeff2D,
								      Patch[Dfrom].N,
								      Patch[Dfrom].y[0],
								      Patch[Dfrom].y.back(),
								      Patch[Dfrom].z[0],
								      Patch[Dfrom].z.back());

	if ((Dj!=0) && (Dk!=(NsubD-1)) && (Patch[Dto3].level==0))
	  Patch[Dto3].psi[kg+(N+1)*(N-jg+ig*(N+1))] =  cheby_2Dinterp(Patch[Dto3].y[N-jg],
								      Patch[Dto3].z[kg],
								      coeff2D,
								      Patch[Dfrom].N,
								      Patch[Dfrom].y[0],
								      Patch[Dfrom].y.back(),
								      Patch[Dfrom].z[0],
								      Patch[Dfrom].z.back());
	
	if ((Dj!=(NsubD-1)) && (Dk!=(NsubD-1)) && (Patch[Dto4].level==0))
	  Patch[Dto4].psi[kg+(N+1)*(jg+ig*(N+1))] =  cheby_2Dinterp(Patch[Dto4].y[jg],
								    Patch[Dto4].z[kg],
								    coeff2D,
								    Patch[Dfrom].N,
								    Patch[Dfrom].y[0],
								    Patch[Dfrom].y.back(),
								    Patch[Dfrom].z[0],
								    Patch[Dfrom].z.back());
	}} // jg, kg loop
    } // ig loop


    // std::chrono::steady_clock::time_point end2D = std::chrono::steady_clock::now();
    // std::cout << "2D interpolations take "  << std::chrono::duration_cast<std::chrono::microseconds>(end2D - start2D).count()*1e-6<< " seconds"<< std::endl;

    // std::chrono::steady_clock::time_point start3D = std::chrono::steady_clock::now();
    // 3D Corners 
//    cheby_3Dinterp_Optimize(Patch,NsubD,Di,Dj,Dk);

    // get Cheby Coeffs
 
    cheby_3Dcoeff(Patch[D].psi,Patch[D].N,Patch[D].coeff,Patch[D]);

    Dfrom = Dk + NsubD*(Dj+(Di)*NsubD);
    Dto1  = Dk + 1 + NsubD*(Dj + 1+(Di + 1)*NsubD);
    Dto2  = Dk - 1 + NsubD*(Dj + 1+(Di + 1)*NsubD);
    Dto3  = Dk + 1 + NsubD*(Dj - 1+(Di + 1)*NsubD);
    Dto4  = Dk + 1 + NsubD*(Dj + 1+(Di - 1)*NsubD);
    Dto5  = Dk - 1 + NsubD*(Dj - 1+(Di + 1)*NsubD);
    Dto6  = Dk + 1 + NsubD*(Dj - 1+(Di - 1)*NsubD);
    Dto7  = Dk - 1 + NsubD*(Dj + 1+(Di - 1)*NsubD);
    Dto8  = Dk - 1 + NsubD*(Dj - 1+(Di - 1)*NsubD);
      
    for (int ig=0; ig< ghost_pts; ig++){
    for (int jg=0; jg< ghost_pts; jg++){
    for (int kg=0; kg< ghost_pts; kg++){

      if ((Di!=(NsubD-1)) && (Dj!=(NsubD-1)) && (Dk!=(NsubD-1)) && (Patch[Dto1].level==0))
	Patch[Dto1].psi[kg+(N+1)*(jg+ig*(N+1))] =  cheby_3Dinterp(Patch[Dto1].x[ig],
								 Patch[Dto1].y[jg],
								 Patch[Dto1].z[kg],
								 Patch[Dfrom].coeff,
								 Patch[Dfrom].N,
								 Patch[Dfrom].x[0],
								 Patch[Dfrom].x.back(),
								 Patch[Dfrom].y[0],
								 Patch[Dfrom].y.back(),
								 Patch[Dfrom].z[0],
								 Patch[Dfrom].z.back());  

      if ((Di!=(NsubD-1)) && (Dj!=(NsubD-1)) && (Dk!=0) && (Patch[Dto2].level==0))
	Patch[Dto2].psi[N-kg+(N+1)*(jg+ig*(N+1))] =  cheby_3Dinterp(Patch[Dto2].x[ig],
								 Patch[Dto2].y[jg],
								 Patch[Dto2].z[N-kg],
								 Patch[Dfrom].coeff,
								 Patch[Dfrom].N,
								 Patch[Dfrom].x[0],
								 Patch[Dfrom].x.back(),
								 Patch[Dfrom].y[0],
								 Patch[Dfrom].y.back(),
								 Patch[Dfrom].z[0],
								 Patch[Dfrom].z.back());  

      if ((Di!=(NsubD-1)) && (Dj!=0) && (Dk!=(NsubD-1)) && (Patch[Dto3].level==0))
	Patch[Dto3].psi[kg+(N+1)*(N-jg+ig*(N+1))] =  cheby_3Dinterp(Patch[Dto3].x[ig],
								 Patch[Dto3].y[N-jg],
								 Patch[Dto3].z[kg],
								 Patch[Dfrom].coeff,
								 Patch[Dfrom].N,
								 Patch[Dfrom].x[0],
								 Patch[Dfrom].x.back(),
								 Patch[Dfrom].y[0],
								 Patch[Dfrom].y.back(),
								 Patch[Dfrom].z[0],
								 Patch[Dfrom].z.back());  

      if ((Di!=0) && (Dj!=(NsubD-1)) && (Dk!=(NsubD-1)) && (Patch[Dto4].level==0))
	Patch[Dto4].psi[kg+(N+1)*(jg+(N-ig)*(N+1))] =  cheby_3Dinterp(Patch[Dto4].x[N-ig],
								 Patch[Dto4].y[jg],
								 Patch[Dto4].z[kg],
								 Patch[Dfrom].coeff,
								 Patch[Dfrom].N,
								 Patch[Dfrom].x[0],
								 Patch[Dfrom].x.back(),
								 Patch[Dfrom].y[0],
								 Patch[Dfrom].y.back(),
								 Patch[Dfrom].z[0],
								 Patch[Dfrom].z.back());  

      if ((Di!=(NsubD-1)) && (Dj!=0) && (Dk!=0) && (Patch[Dto5].level==0))
	Patch[Dto5].psi[N-kg+(N+1)*(N-jg+ig*(N+1))] =  cheby_3Dinterp(Patch[Dto5].x[ig],
								 Patch[Dto5].y[N-jg],
								 Patch[Dto5].z[N-kg],
								 Patch[Dfrom].coeff,
								 Patch[Dfrom].N,
								 Patch[Dfrom].x[0],
								 Patch[Dfrom].x.back(),
								 Patch[Dfrom].y[0],
								 Patch[Dfrom].y.back(),
								 Patch[Dfrom].z[0],
								 Patch[Dfrom].z.back());  

      if ((Di!=0) && (Dj!=0) && (Dk!=(NsubD-1)) && (Patch[Dto6].level==0))
	Patch[Dto6].psi[kg+(N+1)*(N-jg+(N-ig)*(N+1))] =  cheby_3Dinterp(Patch[Dto6].x[N-ig],
								 Patch[Dto6].y[N-jg],
								 Patch[Dto6].z[kg],
								 Patch[Dfrom].coeff,
								 Patch[Dfrom].N,
								 Patch[Dfrom].x[0],
								 Patch[Dfrom].x.back(),
								 Patch[Dfrom].y[0],
								 Patch[Dfrom].y.back(),
								 Patch[Dfrom].z[0],
								 Patch[Dfrom].z.back());  

      if ((Di!=0) && (Dj!=(NsubD-1)) && (Dk!=0) && (Patch[Dto7].level==0))
	Patch[Dto7].psi[N-kg+(N+1)*(jg+(N-ig)*(N+1))] =  cheby_3Dinterp(Patch[Dto7].x[N-ig],
								 Patch[Dto7].y[jg],
								 Patch[Dto7].z[N-kg],
								 Patch[Dfrom].coeff,
								 Patch[Dfrom].N,
								 Patch[Dfrom].x[0],
								 Patch[Dfrom].x.back(),
								 Patch[Dfrom].y[0],
								 Patch[Dfrom].y.back(),
								 Patch[Dfrom].z[0],
								 Patch[Dfrom].z.back());  

      if ((Di!=0) && (Dj!=0) && (Dk!=0) && (Patch[Dto8].level==0))
	Patch[Dto8].psi[N-kg+(N+1)*(N-jg+(N-ig)*(N+1))] =  cheby_3Dinterp(Patch[Dto8].x[N-ig],
								 Patch[Dto8].y[N-jg],
								 Patch[Dto8].z[N-kg],
								 Patch[Dfrom].coeff,
								 Patch[Dfrom].N,
								 Patch[Dfrom].x[0],
								 Patch[Dfrom].x.back(),
								 Patch[Dfrom].y[0],
								 Patch[Dfrom].y.back(),
								 Patch[Dfrom].z[0],
								 Patch[Dfrom].z.back());  

      }}}

    // std::chrono::steady_clock::time_point end3D = std::chrono::steady_clock::now();
    // std::cout << "3D interpolations take "  << std::chrono::duration_cast<std::chrono::microseconds>(end3D - start3D).count()*1e-6<< " seconds"<< std::endl;

    }}} // Di, Dj, Dk loop
//  std::exit(0);
  
//  std::cout << "start printing\n";
  // Print_XYPatch(par,Patch,NsubD,tn-1,0);
  // Print_XYPatch(par,Patch,NsubD,tn-1,1);
  // Print_XYPatch(par,Patch,NsubD,tn-1,2);
  // Print_XYPatch(par,Patch,NsubD,tn-1,3);
  // Print_XYPatch(par,Patch,NsubD,tn-1,4);
  // Print_XYPatch(par,Patch,NsubD,tn-1,5);
  // Print_XYPatch(par,Patch,NsubD,tn-1,6);
  // Print_XYPatch(par,Patch,NsubD,tn-1,7);
//  std::exit(0);
 
}

void Assign_1D3DBC(vector_patch &Patch,
	       signed long int NsubD,
	       Params &par,
	       double tn){
 
  //    std::cout << tn <<"\n";
  //  unsigned long int N = Patch[0].N;

  int N , D, Dfrom;
  vector_complex psi1D(Patch[0].N+1),coeff1D(Patch[0].N+1);
  int ghost_pts = Patch[0].N_ghost;

  // get Cheby Coeffs
  for (int D=0; D<NsubD*NsubD*NsubD; D++)
    cheby_3Dcoeff(Patch[D].psi,Patch[D].N,Patch[D].coeff,Patch[D]);

//  std::exit(0);
  // for (int D=0; D<NsubD*NsubD*NsubD; D++){
  //   N = Patch[D].N;
  //   for (int j = 0; j<N+1; j++){
  //   for (int k = 0; k<N+1; k++){
  //   for (int ig=0; ig< ghost_pts; ig++){
  // 	Patch[D].psi[k+(N+1)*(j+(N+1)*ig)]     = complex(-1e10,-1e10);
  // 	Patch[D].psi[k+(N+1)*(j+(N+1)*ig)]     = complex(-1e10,-1e10);
  // 	Patch[D].psi[k+(N+1)*(j+(N+1)*(N-ig))] = complex(-1e10,-1e10);
  // 	Patch[D].psi[k+(N+1)*(ig+(N+1)*j)]     = complex(-1e10,-1e10);
  // 	Patch[D].psi[k+(N+1)*(N-ig+(N+1)*j)]   = complex(-1e10,-1e10);
  // 	Patch[D].psi[ig+(N+1)*(j+(N+1)*k)]     = complex(-1e10,-1e10);
  // 	Patch[D].psi[N-ig+(N+1)*(j+(N+1)*k)]   = complex(-1e10,-1e10);
  //     }}}}


  // Print_XYPatch(par,Patch,NsubD,tn-1,0);
  // Print_XYPatch(par,Patch,NsubD,tn-1,1);
  // Print_XYPatch(par,Patch,NsubD,tn-1,2);
  // Print_XYPatch(par,Patch,NsubD,tn-1,3);
  // Print_XYPatch(par,Patch,NsubD,tn-1,4);
  // Print_XYPatch(par,Patch,NsubD,tn-1,5);
  // Print_XYPatch(par,Patch,NsubD,tn-1,6);
  // Print_XYPatch(par,Patch,NsubD,tn-1,7);
  // std::exit(0);

  for (int Di=0;Di<NsubD;Di++) {
  for (int Dj=0;Dj<NsubD;Dj++) {
  for (int Dk=0;Dk<NsubD;Dk++) {
    D = Dk + NsubD*(Dj+Di*NsubD);
    N = Patch[D].N;

    // X0 face
    if (Di==0){
      // for (int j = 0; j<N+1; j++){
      // for (int k = 0; k<N+1; k++){
      // for (int ig=0; ig< ghost_pts; ig++){

      // 	Patch[D].psi[k+(N+1)*(j+(N+1)*ig)] = GWP_vals(Patch[D].x[ig],Patch[D].y[j],Patch[D].z[k],(tn-1.)*par.dt);
      // }}}
    }
    else if (Patch[Dk + NsubD*((Di-1)*NsubD+Dj)].level==0){
      Dfrom = Dk + NsubD*(Dj+(Di-1)*NsubD);
      for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
      for (int kg=ghost_pts; kg<N+1-ghost_pts; kg++){
	
	for (int t=0; t< N+1; t++)
	  psi1D[t] = Patch[Dfrom].psi[kg+(N+1)*(jg+t*(N+1))];

	coeff1D = cheby_1Dcoeff(Patch[Dfrom].x, 
				psi1D, 
				Patch[Dfrom].N);

	for (int ig=0; ig< ghost_pts; ig++)
	  Patch[D].psi[kg+(N+1)*(jg+ig*(N+1))] = cheby_1Dinterp(Patch[D].x[ig],
							       coeff1D,
							       Patch[Dfrom].N,
							       Patch[Dfrom].x[0],
							       Patch[Dfrom].x.back());
      }}
    }
    // }
    // else if (Patch[D-1].level==1){
    //   Patch[D].a = Patch[D-1].tmp[Patch[D-1].N-2]; 
    //   Patch[D].psi[0] = Patch[D].a;
    // }

    // XN face
    if (Di==(NsubD-1)){
      // for (int j = 0; j<N+1; j++){
      // for (int k = 0; k<N+1; k++){
      // for (int ig=0; ig< ghost_pts; ig++){
      // 	Patch[D].psi[k+(N+1)*(j+(N-ig)*(N+1))] = GWP_vals(Patch[D].x[N-ig],Patch[D].y[j],Patch[D].z[k],(tn-1.)*par.dt); 
      // 	//	std::cout << "b"<<j<<"="<< Patch[D].b[j] <<"\n";
      // }}}
    }
    else if (Patch[Dk + NsubD*((Di+1)*NsubD+Dj)].level==0){
      Dfrom = Dk + NsubD*(Dj+(Di+1)*NsubD);
      for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
      for (int kg=ghost_pts; kg<N+1-ghost_pts; kg++){
	
	for (int t=0; t< N+1; t++)
	  psi1D[t] = Patch[Dfrom].psi[kg+(N+1)*(jg+t*(N+1))];

	coeff1D = cheby_1Dcoeff(Patch[Dfrom].x, 
				psi1D, 
				Patch[Dfrom].N);

	for (int ig=0; ig< ghost_pts; ig++)
	  Patch[D].psi[kg+(N+1)*(jg+(N-ig)*(N+1))] = cheby_1Dinterp(Patch[D].x[N-ig],
								    coeff1D,
								    Patch[Dfrom].N,
								    Patch[Dfrom].x[0],
								    Patch[Dfrom].x.back());
      }}
    }
    // else if (Patch[D+1].level==1){
    //   Patch[D].b = Patch[D+1].tmp[2]; 
    //   Patch[D].psi[Patch[D].N] = Patch[D].b;
    // }
    //    std::cout << Patch[D].a  << " " << Patch[D].b <<"\n";

    // Y0 face
    if (Dj==0){
      // for (int i = 0; i<N+1; i++){
      // for (int k = 0; k<N+1; k++){
      // for (int jg=0; jg< ghost_pts; jg++){

      // 	Patch[D].psi[k+(N+1)*(jg+i*(N+1))] = GWP_vals(Patch[D].x[i],Patch[D].y[jg],Patch[D].z[k],(tn-1.)*par.dt); 
      // 	//	std::cout << "c"<<i<<"="<< Patch[D].c[i] <<"\n";
      // }}}
    }
    else if (Patch[Dk + NsubD*(Di*NsubD+Dj-1)].level==0){
      Dfrom = Dk + NsubD*(Di*NsubD+Dj-1);
      for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      for (int kg=ghost_pts; kg<N+1-ghost_pts; kg++){
	
	for (int t=0; t< N+1; t++)
	  psi1D[t] = Patch[Dfrom].psi[kg+(N+1)*(t+ig*(N+1))];

	coeff1D = cheby_1Dcoeff(Patch[Dfrom].y, 
				psi1D, 
				Patch[Dfrom].N);

	for (int jg=0; jg< ghost_pts; jg++)
	  Patch[D].psi[kg+(N+1)*(jg+ig*(N+1))] = cheby_1Dinterp(Patch[D].y[jg],
								    coeff1D,
								    Patch[Dfrom].N,
								    Patch[Dfrom].y[0],
								    Patch[Dfrom].y.back());
      }}
    }

    // YN face
    if (Dj==(NsubD-1)){
      // for (int i = 0; i<N+1; i++){
      // for (int k = 0; k<N+1; k++){
      // for (int jg=0; jg< ghost_pts; jg++){
      // 	Patch[D].psi[k+(N+1)*(i*(N+1)+N-jg)] = GWP_vals(Patch[D].x[i],Patch[D].y[N-jg],Patch[D].z[k],(tn-1.)*par.dt); 
      // 	//	std::cout << "d"<<i<<"="<< Patch[D].d[i] <<"\n";
      // }}}
    }
    else if (Patch[Dk + NsubD*(Di*NsubD+Dj+1)].level==0){
      Dfrom = Dk + NsubD*(Di*NsubD+Dj+1);
      for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      for (int kg=ghost_pts; kg<N+1-ghost_pts; kg++){
	
	for (int t=0; t< N+1; t++)
	  psi1D[t] = Patch[Dfrom].psi[kg+(N+1)*(t+ig*(N+1))];

	coeff1D = cheby_1Dcoeff(Patch[Dfrom].y, 
				psi1D, 
				Patch[Dfrom].N);

	for (int jg=0; jg< ghost_pts; jg++)
	  Patch[D].psi[kg+(N+1)*(N-jg+ig*(N+1))] = cheby_1Dinterp(Patch[D].y[N-jg],
								    coeff1D,
								    Patch[Dfrom].N,
								    Patch[Dfrom].y[0],
								    Patch[Dfrom].y.back());
      }}
    }
    // else if (Patch[Di*NsubD+Dj+1].level==0){
    //   Dfrom = Di*NsubD+Dj+1;
    //   for (int i = 0; i<N+1; i++){
    // 	Patch[D].d[i] = Patch[Dfrom].tmp[i*(Patch[Dfrom].N+1)+1]; 
    // 	Patch[D].psi[i*(N+1)+N] = Patch[D].d[i];
    // 	Patch[D].tmp[i*(N+1)+N] = Patch[D].d[i];
    //   }
    // }

    // Z0 face
    if (Dk==0){
      // for (int i = 0; i<N+1; i++){
      // for (int j = 0; j<N+1; j++){
      // for (int kg=0; kg< ghost_pts; kg++){

      // 	Patch[D].psi[kg+(N+1)*(j+i*(N+1))] = GWP_vals(Patch[D].x[i],Patch[D].y[j],Patch[D].z[kg],(tn-1.)*par.dt); 
      // 	//	std::cout << "c"<<i<<"="<< Patch[D].c[i] <<"\n";
      // }}}
    }
    else if (Patch[Dk-1 + NsubD*(Di*NsubD+Dj)].level==0){
      Dfrom = Dk-1 + NsubD*(Di*NsubD+Dj);
      for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
	
	for (int t=0; t< N+1; t++)
	  psi1D[t] = Patch[Dfrom].psi[t+(N+1)*(jg+ig*(N+1))];

	coeff1D = cheby_1Dcoeff(Patch[Dfrom].z, 
				psi1D, 
				Patch[Dfrom].N);

	for (int kg=0; kg< ghost_pts; kg++)
	  Patch[D].psi[kg+(N+1)*(jg+ig*(N+1))] = cheby_1Dinterp(Patch[D].z[kg],
								    coeff1D,
								    Patch[Dfrom].N,
								    Patch[Dfrom].z[0],
								    Patch[Dfrom].z.back());
      }}
    }

    // ZN face
    if (Dk==(NsubD-1)){
      // for (int i = 0; i<N+1; i++){
      // for (int j = 0; j<N+1; j++){
      // for (int kg=0; kg< ghost_pts; kg++){
      // 	Patch[D].psi[N-kg+(N+1)*(j+i*(N+1))] = GWP_vals(Patch[D].x[i],Patch[D].y[j],Patch[D].z[N-kg],(tn-1.)*par.dt); 
      // 	//	std::cout << "d"<<i<<"="<< Patch[D].d[i] <<"\n";
      // }}}
    }
    else if (Patch[Dk+1 + NsubD*(Di*NsubD+Dj)].level==0){
      Dfrom = Dk+1 + NsubD*(Di*NsubD+Dj);
      for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
	
	for (int t=0; t< N+1; t++)
	  psi1D[t] = Patch[Dfrom].psi[t+(N+1)*(jg+ig*(N+1))];

	coeff1D = cheby_1Dcoeff(Patch[Dfrom].z, 
				psi1D, 
				Patch[Dfrom].N);

	for (int kg=0; kg< ghost_pts; kg++)
	  Patch[D].psi[N-kg+(N+1)*(jg+ig*(N+1))] = cheby_1Dinterp(Patch[D].z[N-kg],
								  coeff1D,
								  Patch[Dfrom].N,
								  Patch[Dfrom].z[0],
								  Patch[Dfrom].z.back());
      }}
    }

    }}} // Di, Dj, Dk loop
  //  std::exit(0);
  

  Print_XYPatch(par,Patch,NsubD,tn-1,0);
  Print_XYPatch(par,Patch,NsubD,tn-1,1);
  Print_XYPatch(par,Patch,NsubD,tn-1,2);
  Print_XYPatch(par,Patch,NsubD,tn-1,3);
  Print_XYPatch(par,Patch,NsubD,tn-1,4);
  Print_XYPatch(par,Patch,NsubD,tn-1,5);
  Print_XYPatch(par,Patch,NsubD,tn-1,6);
  Print_XYPatch(par,Patch,NsubD,tn-1,7);
  std::exit(0);
 
}
void Assign_3DBC(vector_patch &Patch,
	       signed long int NsubD,
	       Params &par,
	       double tn){
 
  //    std::cout << tn <<"\n";
  //  unsigned long int N = Patch[0].N;

  int N , D, Dfrom;
  vector_complex psi1D(Patch[0].N+1),coeff1D(Patch[0].N+1);
  int ghost_pts = Patch[0].N_ghost;

  // get Cheby Coeffs
  for (int D=0; D<NsubD*NsubD*NsubD; D++)
    cheby_3Dcoeff(Patch[D].psi,Patch[D].N,Patch[D].coeff,Patch[D]);

  for (int D=0; D<NsubD*NsubD*NsubD; D++){
    N = Patch[D].N;
    for (int j = 0; j<N+1; j++){
    for (int k = 0; k<N+1; k++){
    for (int ig=0; ig< ghost_pts; ig++){
	Patch[D].psi[k+(N+1)*(j+(N+1)*ig)]     = complex(-1e10,-1e10);
	Patch[D].psi[k+(N+1)*(j+(N+1)*ig)]     = complex(-1e10,-1e10);
	Patch[D].psi[k+(N+1)*(j+(N+1)*(N-ig))] = complex(-1e10,-1e10);
	Patch[D].psi[k+(N+1)*(ig+(N+1)*j)]     = complex(-1e10,-1e10);
	Patch[D].psi[k+(N+1)*(N-ig+(N+1)*j)]   = complex(-1e10,-1e10);
	Patch[D].psi[ig+(N+1)*(j+(N+1)*k)]     = complex(-1e10,-1e10);
	Patch[D].psi[N-ig+(N+1)*(j+(N+1)*k)]   = complex(-1e10,-1e10);
    }}}}


  for (int Di=0;Di<NsubD;Di++) {
  for (int Dj=0;Dj<NsubD;Dj++) {
  for (int Dk=0;Dk<NsubD;Dk++) {
    D = Dk + NsubD*(Dj+Di*NsubD);
    N = Patch[D].N;

    // X0 face
    if (Di==0){
      for (int j = 0; j<N+1; j++){
      for (int k = 0; k<N+1; k++){
      for (int ig=0; ig< ghost_pts; ig++){

	Patch[D].psi[k+(N+1)*(j+(N+1)*ig)] = GWP_vals(Patch[D].x[ig],Patch[D].y[j],Patch[D].z[k],(tn-1.)*par.dt);
      }}}
    }
    else if (Patch[Dk + NsubD*((Di-1)*NsubD+Dj)].level==0){
      Dfrom = Dk + NsubD*(Dj+(Di-1)*NsubD);
      for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
      for (int kg=ghost_pts; kg<N+1-ghost_pts; kg++){
      for (int ig=0; ig< ghost_pts; ig++){

	  Patch[D].psi[kg+(N+1)*(jg+ig*(N+1))] =  cheby_3Dinterp(Patch[D].x[ig],
							Patch[D].y[jg],
							Patch[D].z[kg],
							Patch[Dfrom].coeff,
							Patch[Dfrom].N,
							Patch[Dfrom].x[0],
							Patch[Dfrom].x.back(),
							Patch[Dfrom].y[0],
							Patch[Dfrom].y.back(),
							Patch[Dfrom].z[0],
							Patch[Dfrom].z.back());  
      }}}
    }
    // }
    // else if (Patch[D-1].level==1){
    //   Patch[D].a = Patch[D-1].tmp[Patch[D-1].N-2]; 
    //   Patch[D].psi[0] = Patch[D].a;
    // }

    // XN face
    if (Di==(NsubD-1)){
      for (int j = 0; j<N+1; j++){
      for (int k = 0; k<N+1; k++){
      for (int ig=0; ig< ghost_pts; ig++){
	Patch[D].psi[k+(N+1)*(j+(N-ig)*(N+1))] = GWP_vals(Patch[D].x[N-ig],Patch[D].y[j],Patch[D].z[k],(tn-1.)*par.dt); 
//	Patch[D].psi[k+(N+1)*(j+(N-ig)*(N+1))] = GWP_2Dvals(Patch[D].x[N-ig],Patch[D].y[j],(tn-1.)*par.dt); 
//	Patch[D].tmp[k+(N+1)*(j+N*(N+1))] = Patch[D].psi[k+(N+1)*(j+N*(N+1))];
	//	std::cout << "b"<<j<<"="<< Patch[D].b[j] <<"\n";
      }}}
    }
    else if (Patch[Dk + NsubD*((Di+1)*NsubD+Dj)].level==0){
      Dfrom = Dk + NsubD*(Dj+(Di+1)*NsubD);
      // for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
      // for (int kg=ghost_pts; kg<N+1-ghost_pts; kg++){
	
      // 	for (int t=0; t< N+1; t++)
      // 	  psi1D[t] = Patch[Dfrom].psi[kg+(N+1)*(jg+t*(N+1))];

      // 	coeff1D = cheby_1Dcoeff(Patch[Dfrom].x, 
      // 				psi1D, 
      // 				Patch[Dfrom].N);

      // 	for (int ig=0; ig< ghost_pts; ig++)
      // 	  Patch[D].psi[kg+(N+1)*(jg+(N-ig)*(N+1))] = cheby_1Dinterp(Patch[D].x[N-ig],
      // 								    coeff1D,
      // 								    Patch[Dfrom].N,
      // 								    Patch[Dfrom].x[0],
      // 								    Patch[Dfrom].x.back());
      // }}
      for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
      for (int kg=ghost_pts; kg<N+1-ghost_pts; kg++){
      for (int ig=0; ig< ghost_pts; ig++){

	  Patch[D].psi[kg+(N+1)*(jg+(N-ig)*(N+1))] =  cheby_3Dinterp(Patch[D].x[N-ig],
							Patch[D].y[jg],
							Patch[D].z[kg],
							Patch[Dfrom].coeff,
							Patch[Dfrom].N,
							Patch[Dfrom].x[0],
							Patch[Dfrom].x.back(),
							Patch[Dfrom].y[0],
							Patch[Dfrom].y.back(),
							Patch[Dfrom].z[0],
							Patch[Dfrom].z.back());  
      }}}
    }
    // else if (Patch[D+1].level==1){
    //   Patch[D].b = Patch[D+1].tmp[2]; 
    //   Patch[D].psi[Patch[D].N] = Patch[D].b;
    // }
    //    std::cout << Patch[D].a  << " " << Patch[D].b <<"\n";

    // Y0 face
    if (Dj==0){
      for (int i = 0; i<N+1; i++){
      for (int k = 0; k<N+1; k++){
      for (int jg=0; jg< ghost_pts; jg++){

	Patch[D].psi[k+(N+1)*(jg+i*(N+1))] = GWP_vals(Patch[D].x[i],Patch[D].y[jg],Patch[D].z[k],(tn-1.)*par.dt); 
	//	std::cout << "c"<<i<<"="<< Patch[D].c[i] <<"\n";
      }}}
    }
    else if (Patch[Dk + NsubD*(Di*NsubD+Dj-1)].level==0){
      Dfrom = Dk + NsubD*(Di*NsubD+Dj-1);
      // for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      // for (int kg=ghost_pts; kg<N+1-ghost_pts; kg++){
	
      // 	for (int t=0; t< N+1; t++)
      // 	  psi1D[t] = Patch[Dfrom].psi[kg+(N+1)*(t+ig*(N+1))];

      // 	coeff1D = cheby_1Dcoeff(Patch[Dfrom].y, 
      // 				psi1D, 
      // 				Patch[Dfrom].N);

      // 	for (int jg=0; jg< ghost_pts; jg++)
      // 	  Patch[D].psi[kg+(N+1)*(jg+ig*(N+1))] = cheby_1Dinterp(Patch[D].y[jg],
      // 								    coeff1D,
      // 								    Patch[Dfrom].N,
      // 								    Patch[Dfrom].y[0],
      // 								    Patch[Dfrom].y.back());
      // }}
      for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      for (int kg=ghost_pts; kg<N+1-ghost_pts; kg++){
	for (int jg=0; jg< ghost_pts; jg++)
	  Patch[D].psi[kg+(N+1)*(jg+ig*(N+1))] =  cheby_3Dinterp(Patch[D].x[ig],
							Patch[D].y[jg],
							Patch[D].z[kg],
							Patch[Dfrom].coeff,
							Patch[Dfrom].N,
							Patch[Dfrom].x[0],
							Patch[Dfrom].x.back(),
							Patch[Dfrom].y[0],
							Patch[Dfrom].y.back(),
							Patch[Dfrom].z[0],
							Patch[Dfrom].z.back());  
      }}
    }

    // YN face
    if (Dj==(NsubD-1)){
      for (int i = 0; i<N+1; i++){
      for (int k = 0; k<N+1; k++){
      for (int jg=0; jg< ghost_pts; jg++){
	Patch[D].psi[k+(N+1)*(i*(N+1)+N-jg)] = GWP_vals(Patch[D].x[i],Patch[D].y[N-jg],Patch[D].z[k],(tn-1.)*par.dt); 
	//	std::cout << "d"<<i<<"="<< Patch[D].d[i] <<"\n";
      }}}
    }
    else if (Patch[Dk + NsubD*(Di*NsubD+Dj+1)].level==0){
      Dfrom = Dk + NsubD*(Di*NsubD+Dj+1);
      // for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      // for (int kg=ghost_pts; kg<N+1-ghost_pts; kg++){
	
      // 	for (int t=0; t< N+1; t++)
      // 	  psi1D[t] = Patch[Dfrom].psi[kg+(N+1)*(t+ig*(N+1))];

      // 	coeff1D = cheby_1Dcoeff(Patch[Dfrom].y, 
      // 				psi1D, 
      // 				Patch[Dfrom].N);

      // 	for (int jg=0; jg< ghost_pts; jg++)
      // 	  Patch[D].psi[kg+(N+1)*(N-jg+ig*(N+1))] = cheby_1Dinterp(Patch[D].y[N-jg],
      // 								    coeff1D,
      // 								    Patch[Dfrom].N,
      // 								    Patch[Dfrom].y[0],
      // 								    Patch[Dfrom].y.back());
      // }}
      for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      for (int kg=ghost_pts; kg<N+1-ghost_pts; kg++){
      for (int jg=0; jg< ghost_pts; jg++)
	  Patch[D].psi[kg+(N+1)*(N-jg+ig*(N+1))] =  cheby_3Dinterp(Patch[D].x[ig],
							Patch[D].y[N-jg],
							Patch[D].z[kg],
							Patch[Dfrom].coeff,
							Patch[Dfrom].N,
							Patch[Dfrom].x[0],
							Patch[Dfrom].x.back(),
							Patch[Dfrom].y[0],
							Patch[Dfrom].y.back(),
							Patch[Dfrom].z[0],
							Patch[Dfrom].z.back());  
      }}
    }
    // else if (Patch[Di*NsubD+Dj+1].level==0){
    //   Dfrom = Di*NsubD+Dj+1;
    //   for (int i = 0; i<N+1; i++){
    // 	Patch[D].d[i] = Patch[Dfrom].tmp[i*(Patch[Dfrom].N+1)+1]; 
    // 	Patch[D].psi[i*(N+1)+N] = Patch[D].d[i];
    // 	Patch[D].tmp[i*(N+1)+N] = Patch[D].d[i];
    //   }
    // }

    // Z0 face
    if (Dk==0){
      for (int i = 0; i<N+1; i++){
      for (int j = 0; j<N+1; j++){
      for (int kg=0; kg< ghost_pts; kg++){

	Patch[D].psi[kg+(N+1)*(j+i*(N+1))] = GWP_vals(Patch[D].x[i],Patch[D].y[j],Patch[D].z[kg],(tn-1.)*par.dt); 
	//	std::cout << "c"<<i<<"="<< Patch[D].c[i] <<"\n";
      }}}
    }
    else if (Patch[Dk-1 + NsubD*(Di*NsubD+Dj)].level==0){
      Dfrom = Dk-1 + NsubD*(Di*NsubD+Dj);
      // for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      // for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
	
      // 	for (int t=0; t< N+1; t++)
      // 	  psi1D[t] = Patch[Dfrom].psi[t+(N+1)*(jg+ig*(N+1))];

      // 	coeff1D = cheby_1Dcoeff(Patch[Dfrom].z, 
      // 				psi1D, 
      // 				Patch[Dfrom].N);

      // 	for (int kg=0; kg< ghost_pts; kg++)
      // 	  Patch[D].psi[kg+(N+1)*(jg+ig*(N+1))] = cheby_1Dinterp(Patch[D].z[kg],
      // 								    coeff1D,
      // 								    Patch[Dfrom].N,
      // 								    Patch[Dfrom].z[0],
      // 								    Patch[Dfrom].z.back());
      // }}
      for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
      for (int kg=0; kg< ghost_pts; kg++)
	  Patch[D].psi[kg+(N+1)*(jg+ig*(N+1))] =  cheby_3Dinterp(Patch[D].x[ig],
							Patch[D].y[jg],
							Patch[D].z[kg],
							Patch[Dfrom].coeff,
							Patch[Dfrom].N,
							Patch[Dfrom].x[0],
							Patch[Dfrom].x.back(),
							Patch[Dfrom].y[0],
							Patch[Dfrom].y.back(),
							Patch[Dfrom].z[0],
							Patch[Dfrom].z.back());  
      }}
    }

    // ZN face
    if (Dk==(NsubD-1)){
      for (int i = 0; i<N+1; i++){
      for (int j = 0; j<N+1; j++){
      for (int kg=0; kg< ghost_pts; kg++){
	Patch[D].psi[N-kg+(N+1)*(j+i*(N+1))] = GWP_vals(Patch[D].x[i],Patch[D].y[j],Patch[D].z[N-kg],(tn-1.)*par.dt); 
	//	std::cout << "d"<<i<<"="<< Patch[D].d[i] <<"\n";
      }}}
    }
    else if (Patch[Dk+1 + NsubD*(Di*NsubD+Dj)].level==0){
      Dfrom = Dk+1 + NsubD*(Di*NsubD+Dj);
      // for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      // for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
	
      // 	for (int t=0; t< N+1; t++)
      // 	  psi1D[t] = Patch[Dfrom].psi[t+(N+1)*(jg+ig*(N+1))];

      // 	coeff1D = cheby_1Dcoeff(Patch[Dfrom].z, 
      // 				psi1D, 
      // 				Patch[Dfrom].N);

      // 	for (int kg=0; kg< ghost_pts; kg++)
      // 	  Patch[D].psi[N-kg+(N+1)*(jg+ig*(N+1))] = cheby_1Dinterp(Patch[D].z[N-kg],
      // 								  coeff1D,
      // 								  Patch[Dfrom].N,
      // 								  Patch[Dfrom].z[0],
      // 								  Patch[Dfrom].z.back());
      // }}
      for (int ig=ghost_pts; ig<N+1-ghost_pts; ig++){
      for (int jg=ghost_pts; jg<N+1-ghost_pts; jg++){
      for (int kg=0; kg< ghost_pts; kg++)
	  Patch[D].psi[N-kg+(N+1)*(jg+ig*(N+1))] =  cheby_3Dinterp(Patch[D].x[ig],
							Patch[D].y[jg],
							Patch[D].z[N-kg],
							Patch[Dfrom].coeff,
							Patch[Dfrom].N,
							Patch[Dfrom].x[0],
							Patch[Dfrom].x.back(),
							Patch[Dfrom].y[0],
							Patch[Dfrom].y.back(),
							Patch[Dfrom].z[0],
							Patch[Dfrom].z.back());  
      }}
    }

    }}} // Di, Dj, Dk loop
  //  std::exit(0);
  
  Print_XYPatch(par,Patch,NsubD,tn-1,0);
  Print_XYPatch(par,Patch,NsubD,tn-1,1);
  Print_XYPatch(par,Patch,NsubD,tn-1,2);
  Print_XYPatch(par,Patch,NsubD,tn-1,3);
  Print_XYPatch(par,Patch,NsubD,tn-1,4);
  Print_XYPatch(par,Patch,NsubD,tn-1,5);
  Print_XYPatch(par,Patch,NsubD,tn-1,6);
  Print_XYPatch(par,Patch,NsubD,tn-1,7);
  std::exit(0);
 
}
void Assign_AnalyticalBC(vector_patch &Patch,
	       signed long int NsubD,
	       Params &par,
	       double tn){

  //    std::cout << tn <<"\n";
  //  unsigned long int N = Patch[0].N;

  int N , D, Dfrom;
  int ghost_pts = Patch[0].N_ghost;

  for (int Di=0;Di<NsubD;Di++) {
  for (int Dj=0;Dj<NsubD;Dj++) {
  for (int Dk=0;Dk<NsubD;Dk++) {
    D = Dk + NsubD*(Dj+Di*NsubD);
    N = Patch[D].N;

    // X0 face
    if (Di==0){
      for (int j = 0; j<N+1; j++){
      for (int k = 0; k<N+1; k++){
      for (int ig=0; ig< ghost_pts; ig++){
	Patch[D].psi[k+(N+1)*(j+(N+1)*ig)] = GWP_vals(Patch[D].x[ig],Patch[D].y[j],Patch[D].z[k],(tn-1.)*par.dt);
      }}}
    }

    // XN face
    if (Di==(NsubD-1)){
      for (int j = 0; j<N+1; j++){
      for (int k = 0; k<N+1; k++){
      for (int ig=0; ig< ghost_pts; ig++){
	Patch[D].psi[k+(N+1)*(j+(N-ig)*(N+1))] = GWP_vals(Patch[D].x[N-ig],Patch[D].y[j],Patch[D].z[k],(tn-1.)*par.dt); 
//	Patch[D].psi[k+(N+1)*(j+(N-ig)*(N+1))] = GWP_2Dvals(Patch[D].x[N-ig],Patch[D].y[j],(tn-1.)*par.dt); 
//	Patch[D].tmp[k+(N+1)*(j+N*(N+1))] = Patch[D].psi[k+(N+1)*(j+N*(N+1))];
	//	std::cout << "b"<<j<<"="<< Patch[D].b[j] <<"\n";
      }}}
    }

    // Y0 face
    if (Dj==0){
      for (int i = 0; i<N+1; i++){
      for (int k = 0; k<N+1; k++){
      for (int jg=0; jg< ghost_pts; jg++){

	Patch[D].psi[k+(N+1)*(jg+i*(N+1))] = GWP_vals(Patch[D].x[i],Patch[D].y[jg],Patch[D].z[k],(tn-1.)*par.dt); 
	//	std::cout << "c"<<i<<"="<< Patch[D].c[i] <<"\n";
      }}}
    }

    // YN face
    if (Dj==(NsubD-1)){
      for (int i = 0; i<N+1; i++){
      for (int k = 0; k<N+1; k++){
      for (int jg=0; jg< ghost_pts; jg++){
	Patch[D].psi[k+(N+1)*(i*(N+1)+N-jg)] = GWP_vals(Patch[D].x[i],Patch[D].y[N-jg],Patch[D].z[k],(tn-1.)*par.dt); 
	//	std::cout << "d"<<i<<"="<< Patch[D].d[i] <<"\n";
      }}}
    }

    // Z0 face
    if (Dk==0){
      for (int i = 0; i<N+1; i++){
      for (int j = 0; j<N+1; j++){
      for (int kg=0; kg< ghost_pts; kg++){

	Patch[D].psi[kg+(N+1)*(j+i*(N+1))] = GWP_vals(Patch[D].x[i],Patch[D].y[j],Patch[D].z[kg],(tn-1.)*par.dt); 
	//	std::cout << "c"<<i<<"="<< Patch[D].c[i] <<"\n";
      }}}
    }

    // ZN face
    if (Dk==(NsubD-1)){
      for (int i = 0; i<N+1; i++){
      for (int j = 0; j<N+1; j++){
      for (int kg=0; kg< ghost_pts; kg++){
	Patch[D].psi[N-kg+(N+1)*(j+i*(N+1))] = GWP_vals(Patch[D].x[i],Patch[D].y[j],Patch[D].z[N-kg],(tn-1.)*par.dt); 
	//	std::cout << "d"<<i<<"="<< Patch[D].d[i] <<"\n";
      }}}
    }

    }}} // Di, Dj, Dk loop
  //  std::exit(0);
}
void AMR_Scheme(Patches &P){
  double threshold = 0.1;

  double ave = 0.0;
  for (size_t j = 1; j < P.N-1; ++j)
    ave += std::abs(P.psi[j]);

  ave /= P.N;

  // if ((ave>threshold) && (P.level==0)) {
  //   Refinement(P);
  //   std::cout <<"refined?\n";
  // }
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

  // initialize
//  Assign_3DBC(Patch,NsubD,par,tn);
  
  Assign_1D2D3DBC(Patch,NsubD,par,tn);
//  Assign_AnalyticalBC(Patch,NsubD,par,tn);
//  std::exit(0);
  //  for (int D=0;D<NsubD;D++) AMR_Scheme(Patch[D]);

//  if ((int)tn%10==0) Print_XYPlane(par,Patch,NsubD,tn-1);

    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

  // x direction
  for (int D=0;D<NsubD*NsubD*NsubD;D++)
    for(int i=0; i<(Patch[D].N+1) * (Patch[D].N+1) * (Patch[D].N+1);i++)
      Patch[D].k1[i] = complex(0.0,0.0);

  for (int D=0;D<NsubD*NsubD*NsubD;D++) {
    for(int i=0; i<Patch[D].N+1;i++){ 
    for(int j=0; j<Patch[D].N+1;j++){
    for(int k=0; k<Patch[D].N+1;k++){
 
      for(unsigned long int t=0; t<Patch[D].N+1; t++)
	Patch[D].k1[j+(Patch[D].N+1)*(i+k*(Patch[D].N+1))] += Patch[D].Dc[k*(Patch[D].N+1)+t] * Patch[D].psi[j+(Patch[D].N+1)*(i+t*(Patch[D].N+1))];
    
    }}}
  }

  for (int D=0;D<NsubD*NsubD*NsubD;D++)
    for(int i=0; i<(Patch[D].N+1) * (Patch[D].N+1) * (Patch[D].N+1);i++)
      Patch[D].psi[i] += Patch[D].k1[i];

  // y direction
  for (int D=0;D<NsubD*NsubD*NsubD;D++)
    for(int i=0; i<(Patch[D].N+1) * (Patch[D].N+1) * (Patch[D].N+1);i++)
      Patch[D].k1[i] = complex(0.0,0.0);

  for (int D=0;D<NsubD*NsubD*NsubD;D++) {
    for(int i=0; i<Patch[D].N+1;i++){
    for(int j=0; j<Patch[D].N+1;j++){
    for(int k=0; k<Patch[D].N+1;k++){

      for(unsigned long int t=0; t<Patch[D].N+1; t++)
	Patch[D].k1[j+(Patch[D].N+1)*(k+i*(Patch[D].N+1))] += Patch[D].Dc[k*(Patch[D].N+1)+t] * Patch[D].psi[j+(Patch[D].N+1)*(t+i*(Patch[D].N+1))];
    
    }}}
  }

  for (int D=0;D<NsubD*NsubD*NsubD;D++)
    for(int i=0; i<(Patch[D].N+1) * (Patch[D].N+1) * (Patch[D].N+1);i++)
      Patch[D].psi[i] += Patch[D].k1[i];

  // z direction
  for (int D=0;D<NsubD*NsubD*NsubD;D++)
    for(int i=0; i<(Patch[D].N+1) * (Patch[D].N+1)* (Patch[D].N+1);i++)
      Patch[D].k1[i] = complex(0.0,0.0);

  for (int D=0;D<NsubD*NsubD*NsubD;D++){
    for(int i=0; i<Patch[D].N+1;i++){
    for(int j=0; j<Patch[D].N+1;j++){
    for(int k=0; k<Patch[D].N+1;k++){

      for(unsigned long int t=0; t<Patch[D].N+1; t++)
	Patch[D].k1[k + (Patch[D].N+1)*(j + i*(Patch[D].N+1))] += Patch[D].Dc[k*(Patch[D].N+1)+t] * Patch[D].psi[ t + (Patch[D].N+1)*(j + i*(Patch[D].N+1))];    
    }}}
    }

  for (int D=0;D<NsubD*NsubD*NsubD;D++)
    for(int i=0; i<(Patch[D].N+1) * (Patch[D].N+1)* (Patch[D].N+1);i++)
      Patch[D].psi[i] += Patch[D].k1[i];


    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//    std::cout << " exponential evolv take "  << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()*1e-6<< " seconds \n\n"<< std::endl;

    // Kn[k+(P.N+1)*(j+i*(P.N+1))] += P.D[k*(P.N+1)+t] * P.tmp[t+(P.N+1)*(j+i*(P.N+1))] ;
    // Kn[j+(P.N+1)*(k+i*(P.N+1))] += P.D[k*(P.N+1)+t] * P.tmp[j+(P.N+1)*(t+i*(P.N+1))] ;
    // Kn[j+(P.N+1)*(i+k*(P.N+1))] += P.D[k*(P.N+1)+t] * P.tmp[j+(P.N+1)*(i+t*(P.N+1))] ;
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
    //    inp_dt =5e-3/inp_timestep; // for v=1e1 
    inp_dt =1e1/inp_timestep; // 
//    inp_dt =1e-1/inp_timestep; //  for 2D

    N += ghost_pts*2;
  }
  std::cout << inp_timestep << " " << N << " " << NsubD << " " << ghost_pts << "\n";
  
  Params par = Params(inp_xmax, inp_res, inp_dt, inp_timestep, inp_np, inp_im, inp_xoffset); // xmax, res, dt, timestep, np, im

  Operators opr = Operators(par, 0.0, inp_xoffset);

    
  double Coeff = 0.5 * par.dt / (par.dx*par.dx);

  double a = 0.0;
  double b = 1.0;
  double L = b-a;

  // sub domain parameters
  //  unsigned long int Nsub = 20 ; // doesnt matter if no refinement
  vector_patch Patch;
  
  for (int i=0;i<NsubD*NsubD*NsubD;i++) {
    Patch.push_back(Patches(N));
    Patch[i].L = L;
    Patch[i].N_ghost = ghost_pts;
    Patch[i].dt = par.dt;
  }

  // create subdomain x
  for (int Di=0;Di<NsubD;Di++){  
  for (int Dj=0;Dj<NsubD;Dj++){
  for (int Dk=0;Dk<NsubD;Dk++){
    Patch[Dk+NsubD*(Dj+NsubD*Di)].x = cheby_exgrid_oni(N,
					    a+L*Di/NsubD,
					    a+L*(Di+1)/NsubD,
					    ghost_pts);  
    Patch[Dk+NsubD*(Dj+NsubD*Di)].y = cheby_exgrid_oni(N,
					    a+L*Dj/NsubD,
					    a+L*(Dj+1)/NsubD,
					    ghost_pts);  
    Patch[Dk+NsubD*(Dj+NsubD*Di)].z = cheby_exgrid_oni(N,
					    a+L*Dk/NsubD,
					    a+L*(Dk+1)/NsubD,
					    ghost_pts);  
  }}}


  double dx_sub = Patch[0].x[1]-Patch[0].x[0];

  // make cheby D matrix
  for (int D=0;D<NsubD*NsubD*NsubD;D++)  cheby_expDn_matrix(Patch[D].D, 
							    Patch[D].Dc,
							    Patch[D].x,
							    Patch[D].N, 
							    par.dt);

  // initialize GWP
  for (int D=0; D<NsubD*NsubD*NsubD;D++)   initial_GWP(Patch[D].x,Patch[D].y,Patch[D].z,Patch[D].psi);

//  Print_XYPlane(par,Patch,NsubD,0);
//  Print_XZPlane(par,Patch,NsubD,0);
//  Print_YZPlane(par,Patch,NsubD,0);
  std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();

  std::cout << inp_timestep << " CFL constant = "<< par.dt/(dx_sub * dx_sub) << "\n";
  Print_XYZ(par,Patch,NsubD,0);
  for (int tn = 1; tn <= par.timesteps; ++tn){
  std::chrono::steady_clock::time_point begin3 = std::chrono::steady_clock::now();

    // RK4
    Time_evolv(Patch,
	       NsubD,
	       par,
	       tn);
//    std::exit(0);

   // if (tn%10==0){ 
   //   Print_XYPlane(par,Patch,NsubD,tn);
   //   Print_XZPlane(par,Patch,NsubD,tn);
   //   Print_YZPlane(par,Patch,NsubD,tn);
   // }
    //    if (tn%1000==0) Print_file(par,Patch,NsubD,tn);
    std::chrono::steady_clock::time_point end3 = std::chrono::steady_clock::now();
    std::cout << tn << " ";
    std::cout << "step takes = " << std::chrono::duration_cast<std::chrono::microseconds>(end3 - begin3).count()*1e-6<< " seconds"<< std::endl;
    if (tn%100==0) Print_XYZ(par,Patch,NsubD,tn);
    Print_XYPatch(par,Patch,NsubD,tn,0);
//    Print_XYPatch(par,Patch,NsubD,tn,1);
//    Print_XYPatch(par,Patch,NsubD,tn,2);
//    Print_XYPatch(par,Patch,NsubD,tn,3);
//    Print_XYPatch(par,Patch,NsubD,tn,4);
    Print_XYPatch(par,Patch,NsubD,tn,5);
//    Print_XYPatch(par,Patch,NsubD,tn,6);
//    Print_XYPatch(par,Patch,NsubD,tn,7);
  
 } //   for (int tn = 0; tn <= par.timesteps; ++tn){
  std::cout << "";

  double err = 0.0;
  double Maxerr = 0.0;
  double NN = NsubD*(N-ghost_pts*2); 

  //  for (int D=0;D<NsubD;D++)   Print(par,Patch[D].x,Patch[D].psi,par.timesteps);
//  Print_XYPlane(par,Patch,NsubD,par.timesteps);
//  Print_XZPlane(par,Patch,NsubD,par.timesteps);
//  Print_YZPlane(par,Patch,NsubD,par.timesteps);

  for (int D=0;D<NsubD*NsubD*NsubD;D++)   Error(par,Patch[D].x,Patch[D].y,Patch[D].z,Patch[D].psi,err);
  for (int D=0;D<NsubD*NsubD*NsubD;D++)   MaxError(par,Patch[D].x,Patch[D].y,Patch[D].z,Patch[D].psi,Maxerr);

  std::cout << "Error =" << pow(err/pow(NN,3),1./2) << "\n" << std::flush;
  std::cout << "MaxError =" << Maxerr << "\n" << std::flush;
  std::cout << NN <<"\n" << std::flush;;

  Print_XYZ(par,Patch,NsubD,par.timesteps);
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "simulation takes = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin2).count()*1e-6<< " seconds"<< std::endl;

return 0;
}
