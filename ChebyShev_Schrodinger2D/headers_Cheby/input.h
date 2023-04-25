const double m2kpc = 3.2408e19; // m  
const double s2Myr = 60.0*60.0*24.0*365.0*1e6; // s

/* const double hbar = 6.58211928e-16 / s2Myr; // eV * Myr */
/* const double mfdm = 2.5e-22; // eV */
/* const double bigG = 4.51475758330181e-48 * pow(s2Myr,2); // Mpc**3/Msolar/Myr**2 */
/* const double c = 2.99792458e8 / m2kpc * s2Myr; // kpc/Myr */

const double hbar = 1.0;
const double mfdm = 1.0;
const double bigG = 1.0;
const double c = 1.0;

const double ivR = 1.0; 
const double ivL = -1.0;
const double Mh1 = 1e9;
const double Mh2 = 1e9;
const double v0_init = -2e1;  //1e1
const double x0_init = 0.5;

double inp_xmax = 1.0;
double inp_res = 64;
double inp_dt = 1e-5;
double inp_timestep = 5e3;
int inp_np = 10;
bool inp_im = false;
double inp_xoffset = 3.0;
