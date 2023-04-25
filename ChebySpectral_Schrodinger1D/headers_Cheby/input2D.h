const double m2kpc = 3.2408e19; // m  
const double s2Myr = 60.0*60.0*24.0*365.0*1e6; // s
const double hbar = 6.58211928e-16 / s2Myr; // eV * Myr
const double mfdm = 2.5e-22; // eV
const double bigG = 4.51475758330181e-48 * pow(s2Myr,2); // Mpc**3/Msolar/Myr**2
const double c = 2.99792458e8 / m2kpc * s2Myr; // kpc/Myr

const double ivR = 0.0; 
const double ivL = 0.0;
const double Mh1 = 8e9;
const double Mh2 = 1e9;

double inp_xmax = 20.0;
double inp_res = 256*2;
double inp_dt = 2e0;
double inp_timestep = 2000.0;
int inp_np = 2000;
bool inp_im = false;
double inp_offset[2] = {3.0,0.0};
