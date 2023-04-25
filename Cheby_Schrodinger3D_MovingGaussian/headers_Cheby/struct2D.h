extern const double mfdm, bigG, ivR, ivL, Mh1, Mh2;
struct Params
{      
  double xmax;
  unsigned int res;
  double dt;
  unsigned int timesteps;
  double dx, dy;
  vector_real x, y;
  double dkx, dky;
  vector_real kx, ky;
  bool im_time;
  int np;
  vector_real px, py;
  vector_real pvx, pvy;
  double rn, rp, rs, del, v0; // Sponge parameter

  // initialization
  Params(double _xmax, unsigned int _res, double _dt, unsigned int _timesteps, int _np, bool im, double offset[2])
  {

    xmax = _xmax;
    res = _res;
    dt = _dt;
    timesteps = _timesteps;
    dx = 2.0 * xmax / res;
    dy = 2.0 * xmax / res;
    x.reserve(res);
    y.reserve(res);
    dkx = M_PI / xmax;
    dky = M_PI / xmax;
    kx.reserve(res);
    ky.reserve(res);
    im_time = im;
    np = _np;
    px.reserve(np*2);
    py.reserve(np*2);
    pvx.reserve(np*2);
    pvy.reserve(np*2);

    // Sponge parameter
    rs = xmax; // center of sponge
    del = xmax*0.2; // width of sponge
    v0 = 1e-2;
    
    for (size_t i = 0; i < res; ++i)
    {
      x.emplace_back(xmax / res - xmax + i * (2.0 * xmax / res));
      y.emplace_back(xmax / res - xmax + i * (2.0 * xmax / res));

      if (i < res / 2)
      {
	  kx.push_back(i * M_PI / xmax);
	  ky.push_back(i * M_PI / xmax);
      }
      else
      {
	  kx.push_back((static_cast<double>(i) - res) * M_PI / xmax);
	  ky.push_back((static_cast<double>(i) - res) * M_PI / xmax);
      }
    }

    // N-body particles x,v
    srand(time(NULL));
    double rr, ang;
    for (int i = 0; i < np; ++i)
    {
//      px.push_back(-xmax + i * (2.0*xmax)/(double)np);
//      rr = sqrt(pow(fRand(0,1.0),-1./8))/0.091/2.0;
      rr = sqrt(fRand(0,1.0))/2.0;  // random points within circle
      ang = fRand(0,2*M_PI);

      px.push_back(rr*cos(ang)+offset[0]);
      py.push_back(rr*sin(ang)+offset[1]);
      px.push_back(rr*cos(ang)-offset[0]);
      py.push_back(rr*sin(ang)-offset[1]);
      pvx.push_back(0.0);
      pvy.push_back(0.0);
      pvx.push_back(0.0);
      pvy.push_back(0.0);

    }
  }

};

double SolitonDensity (double r, double M )
{
  double rc = 1.6 * pow(M/1e9, -1./3) * (1e-22/mfdm); // kpc
  double rho0 = 3.1e15 * pow(2.5e-22/mfdm,2) * pow(rc,-4); // Msolar/Mpc**3
  return rho0 * pow(1+0.091*pow(r/rc, 2), -8); // Msolar/Mpc**3
}

struct Operators
{
  size_t size;
  vector_real phi, oldphi;  // potential
  vector_real ax, oldax;
  vector_real ay, olday;
  vector_complex pe; // opr.r
  vector_complex ke; // opr.k
  vector_complex wfc; // wavefunction

  public:
  Operators(Params &par, double voffset, double wfcoffset[2])
  {
    size = par.res;
    phi.reserve(size);
    oldphi.reserve(size);
    ax.reserve(size);
    ay.reserve(size);
    oldax.reserve(size);
    olday.reserve(size);
    pe.reserve(size);
    ke.reserve(size);
    wfc.reserve(size);
    for (size_t i = 0; i < size; ++i)
    {
      for (size_t j = 0; j < size; ++j)
      {
	phi.push_back(0.0);
	oldphi.push_back(0.0);
	ax.push_back(0.0);
	ay.push_back(0.0);
	oldax.push_back(0.0);
	olday.push_back(0.0);
	pe.push_back(complex(0.0,0.0));
	ke.push_back(complex(0.0,0.0));            
	
	//wfc.push_back(sqrt(SolitonDensity(sqrt(pow(par.x[i] - wfcoffset[0], 2) + pow(par.y[j] - wfcoffset[1], 2)), Mh1)) * exp(ivL*(par.x[i]) * complex(0.0, 1.0)));
     // wfc.push_back(1.0*exp(-pow(par.x[i] - wfcoffset, 2) / 2.0) * exp(ivR*(par.x[i] - wfcoffset)* complex(0.0, 1.0))/1.0 + exp(-pow(par.x[i] + wfcoffset, 2) / 2.0) * exp(ivL*(par.x[i] + wfcoffset)* complex(0.0, 1.0)/1.0)); 
	wfc.push_back(sqrt(SolitonDensity(sqrt(pow(par.x[i] - wfcoffset[0], 2) + pow(par.y[j] - wfcoffset[1], 2)), Mh1)) * exp(ivL*(par.x[i]) * complex(0.0, 1.0)) + sqrt(SolitonDensity(sqrt(pow(par.x[i] + wfcoffset[0], 2) + pow(par.y[j] + wfcoffset[1], 2)), Mh2)) * exp(ivR*(par.x[i]) * complex(0.0, 1.0)));
//      std::cout << pow(abs(wfc[i]),2) << "\t";
      }
    }
  };
};

