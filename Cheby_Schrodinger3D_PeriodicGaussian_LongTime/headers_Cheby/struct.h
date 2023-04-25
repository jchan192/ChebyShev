#include <fftw3.h>

extern const double mfdm, bigG, ivR, ivL, Mh1, Mh2;
struct Params
{      
  double xmax;
  unsigned int res;
  double dt;
  unsigned int timesteps;
  double dx;
  vector_real x,y;
  double dk;
  vector_real k;
  bool im_time;
  int np;
  vector_real px;
  vector_real pv;
  double rn, rp, rs, del, v0; // Sponge parameter
  double a_R,a_I,b_R,b_I;
  double aold_R,aold_I,bold_R,bold_I;
  double as_R, bs_R, asold_R,bsold_R; 
  double as_I, bs_I, asold_I,bsold_I; 
  // initialization
  Params(double _xmax, unsigned int _res, double _dt, unsigned int _timesteps, int _np, bool im, double inp_xoffset)
  {

    xmax = _xmax;
    res = _res;
    dt = _dt;
    timesteps = _timesteps;
    dx = xmax / res;
    x.reserve(res);
    np = _np;
    px.reserve(np);
    pv.reserve(np);
    dk = M_PI / xmax;
    k.reserve(res);
    im_time = im;

    // standard chebyshev is -1 to 1
    double a = 0.0;
    double b = 1.0;

    for (size_t i = 0; i < res+1; ++i)
    {
      //      x.emplace_back(xmax / res - xmax + i * (2.0 * xmax / res));
      //      x.emplace_back(0.5 * xmax / res + i * (xmax / res));
      //      x.emplace_back(-cos(i*M_PI/res));
      x.emplace_back(-cos(i*M_PI/res)*(b-a)*0.5 + (b+a)*0.5);
      y.emplace_back(-cos(i*M_PI/res)*(b-a)*0.5 + (b+a)*0.5);

      if (i < res / 2)
      {
	  k.push_back(i * M_PI / xmax);
      }
      else
      {
	  k.push_back((static_cast<double>(i) - res) * M_PI / xmax);
      }
    }

    // N-body particles x,v
    srand(time(NULL));
    for (int i = 0; i < np; ++i)
    {
//      px.push_back(-xmax + i * (2.0*xmax)/(double)np);
//      px.push_back(sqrt(pow(fRand(0,1.0),-1./8)-1.0)/0.091);
      pv.push_back(0.0);
    }
  }

};

double SolitonDensity (double x, double M )
{
  double rc = 1.6 * pow(M/1e9, -1./3) * (1e-22/mfdm); // kpc
  double rho0 = 3.1e15 * pow(2.5e-22/mfdm,2) * pow(rc,-4); // Msolar/Mpc**3
  return rho0 * pow(1+0.091*pow(x/rc, 2), -8); // Msolar/Mpc**3
}

struct Operators
{
  size_t size;
  vector_real phi;  // potential
  vector_real oldphi; 
  vector_complex pe; // opr.r
  vector_complex ke; // opr.k
  //  vector_complex wfc; // wavefunction
  vector_complex wfc,k1,k2,k3,k4,k5,k6;

  public:
  Operators(Params &par, double voffset, double wfcoffset)
  {
    size = par.res;
    phi.reserve(size);
    oldphi.reserve(size);
    pe.reserve(size);
    ke.reserve(size);
    wfc.reserve(size+1);
    for (size_t i = 0; i < size+1; ++i)
    {
      phi.push_back(0.0);
      oldphi.push_back(0.0);
      k1.push_back(complex(0.0,0.0));
      k2.push_back(complex(0.0,0.0));            
      k3.push_back(complex(0.0,0.0));            
      k4.push_back(complex(0.0,0.0));            
      k5.push_back(complex(0.0,0.0));            
      k6.push_back(complex(0.0,0.0));            
      wfc.push_back(complex(0.0,0.0));            
      //      ke[i] = exp(-0.5 * par.dt * hbar/(mfdm/pow(c,2)) * pow(par.k[i], 2) * complex(0.0, 1.0));

//      wfc.push_back(exp(-pow(par.x[i] - wfcoffset, 2) / 2.0));
//      wfc.push_back(exp(-pow(par.x[i] - wfcoffset, 2) / 1.0) * exp(ivR*(par.x[i] - wfcoffset)* complex(0.0, 1.0))/1.0 + exp(-pow(par.x[i] + wfcoffset, 2) / 1.0) * exp(ivL*(par.x[i] + wfcoffset)* complex(0.0, 1.0)/1.0)); 
     // wfc.push_back(sqrt(SolitonDensity(par.x[i] - wfcoffset, Mh1)) * exp(ivL*(par.x[i]) * complex(0.0, 1.0))+ sqrt(SolitonDensity(par.x[i] + wfcoffset, Mh2)) * exp(ivR*(par.x[i]) * complex(0.0, 1.0)));
//      std::cout << pow(abs(wfc[i]),2) << "\t";
    }
  }
};
struct Patches
{
  unsigned long int Nsub,N;
  vector_real phi;  // potential
  vector_real x,y,z; 
  vector_real D,Dsub; 
  vector_complex psi,tmp,k1,k2,k3,k4,debug; // 
  double L;
  int level;
  vector_complex coeff;
  vector_real coeff2D_R, coeff2D_I;
  vector_real coeff3D_R, coeff3D_I;
  vector_complex Dc;
  double dt;
  int N_ghost;
  fftw_plan plan2D_R, plan2D_I,plan3D_R, plan3D_I;

  public:
  Patches(unsigned long int _N)
  {

    level = 0;
    N = _N;
    //    Nsub = _Nsub;
    phi.reserve(N+1);
    psi.reserve(N+1);
    k1.reserve(N+1);
    k2.reserve(N+1);
    k3.reserve(N+1);
    k4.reserve(N+1);
    tmp.reserve(N+1);
    debug.reserve(N+1);
    x.reserve(N+1);
    D.reserve((N+1)*(N+1));

    for (size_t i = 0; i <(N+1); ++i)
    {
      x.push_back(0.0);
      y.push_back(0.0);
      z.push_back(0.0);
    }
    for (size_t i = 0; i <(N+1)* (N+1)*(N+1); ++i)
    {
      phi.push_back(0.0);
      psi.push_back(complex(0.0,0.0));            
      tmp.push_back(complex(0.0,0.0));            
      k1.push_back(complex(0.0,0.0));
      k2.push_back(complex(0.0,0.0));
      k3.push_back(complex(0.0,0.0));
      k4.push_back(complex(0.0,0.0));
      debug.push_back(complex(0.0,0.0));            
      coeff.push_back(complex(0.0,0.0));            
      coeff3D_R.push_back(0.0); 
      coeff3D_I.push_back(0.0);
    }

    for (size_t i = 0; i <(N+1)*(N+1); ++i){
      D.push_back(0.0);
      Dc.push_back(complex(0.0,0.0));
      coeff2D_R.push_back(0.0); 
      coeff2D_I.push_back(0.0);
    }
    /* for (size_t i = 0; i < (Nsub+1)*(Nsub+1); ++i) */
    /*   Dsub.push_back(0.0); */
  plan2D_R = fftw_plan_r2r_2d(N+1,N+1, &coeff2D_R[0], &coeff2D_R[0], FFTW_REDFT00, FFTW_REDFT00,FFTW_ESTIMATE);
  plan2D_I = fftw_plan_r2r_2d(N+1,N+1, &coeff2D_I[0], &coeff2D_I[0], FFTW_REDFT00, FFTW_REDFT00,FFTW_ESTIMATE);
  plan3D_R = fftw_plan_r2r_3d(N+1,N+1,N+1, &coeff3D_R[0], &coeff3D_R[0], FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00,FFTW_ESTIMATE);
  plan3D_I = fftw_plan_r2r_3d(N+1,N+1,N+1, &coeff3D_I[0], &coeff3D_I[0], FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00,FFTW_ESTIMATE);
  }
};

