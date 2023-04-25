import numpy as np
import pyfftw.interfaces.numpy_fft as fft
from scipy.fftpack import dct,dctn

# def GWP_vals(x,y,Time):

  
#   double hbar = 1.;
#   double m =1.;
#   double  sigma_px = 8, sigma_py=8;
#   //  double  vx0 =  v0_init ; //2e1;
#   double  vx0 =  -5e0,  vy0 = -5e0; //2e1;
#   double  x0 = 0.5;
#   double  y0 = 0.5;

#   double sigma_x, phase_x, M_x;
#   double sigma_y, phase_y, M_y;

#   sigma_x = hbar/(2*sigma_px) * pow(1+4*pow(sigma_px,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
#   phase_x = 1./hbar * (m*vx0+sigma_px*sigma_px/(sigma_x*sigma_x)*Time/(2*m)*(x-x0-vx0*Time))*
#     (x-x0-vx0*Time) + vx0*m/(2*hbar) * vx0*Time - atan(2*sigma_px*sigma_px*Time/(hbar*m))/ 2 ;
#   M_x = 1./(pow(2*M_PI,1./4)*pow(sigma_x,1./2)) * exp(-(x-x0-vx0*Time)*(x-x0-vx0*Time)/(4*sigma_x*sigma_x)) ;

#   sigma_y = hbar/(2*sigma_py) * pow(1+4*pow(sigma_py,4)/(hbar*hbar)*Time*Time/(m*m),1./2);
#   phase_y = 1./hbar * (m*vy0+sigma_py*sigma_py/(sigma_y*sigma_y)*Time/(2*m)*(y-y0-vy0*Time))*
#     (y-y0-vy0*Time) + vy0*m/(2*hbar) * vy0*Time - atan(2*sigma_py*sigma_py*Time/(hbar*m))/ 2 ;
#   M_y = 1./(pow(2*M_PI,1./4)*pow(sigma_y,1./2)) * exp(-(y-y0-vy0*Time)*(y-y0-vy0*Time)/(4*sigma_y*sigma_y));

#   return complex(M_x*M_y*cos(phase_x+phase_y),M_x*M_y*sin(phase_x+phase_y));

def cheby_exgrid(N,a,b):
    i=np.linspace(0,N,N+1)
    return (b-a)/2.*(-np.cos(i*np.pi/N))+(b+a)/2.0

def cheby_poly(x,n):
    # clenshaw algorithm
    T0 = 1.0
    T1 = x
    if n == 0:
        return  T0
    elif n == 1:
        return  T1
    elif n > 1:
        for i in range(n-1):
            T = 2*x*T1-T0 
            T0 = T1
            T1 = T
        return T

def cheby_coeff_2d(Z, N):
    # Compute Chebyshev nodes
#    x = -np.cos(np.pi * np.linspace(0,N,N+1) / N )
#    y = -np.cos(np.pi * np.linspace(0,N,N+1) / N )

    # Compute function values at Chebyshev nodes
#    X, Y = np.meshgrid(x, y)
#    Z = f(X, Y)
    print(Z.shape)
    # Compute 2D DCT of function values
    c = dctn(Z, type=1)
    
    # Scale coefficients to obtain Chebyshev polynomial coefficients
    a = np.zeros((N+1, N+1))
    a = c / (N**2)

    a[0,:] *= 0.5
    a[:,0] *= 0.5
    a[N,1:N] *= 0.5
    a[1:N,N] *= 0.5

#    print(a)
    return a

def cheby_interp_2d(x,y,coeff,N,ax,bx,ay,by):

    x_chebgrid = -( 2.0 * x - ax  - bx ) / ( bx-ax);
    y_chebgrid = -( 2.0 * y - ay  - by ) / ( by-ay);

    sum = 0
    for i in range(0,N+1):
        for j in range(0,N+1):
                sum += coeff[j,i]*cheby_poly(y_chebgrid,i)*cheby_poly(x_chebgrid,j)
    return sum

Coeff_data = np.loadtxt('output/P4_1.dat',usecols=[6,])
F = np.loadtxt('output/P4_1.dat',usecols=[4,])
xx = np.loadtxt('output/P4_1.dat',usecols=[0,]) 
yy = np.loadtxt('output/P4_1.dat',usecols=[1,]) 
N = np.int64(np.sqrt(xx.shape[0])-1)

ax = np.min(xx)
bx = np.max(xx)
ay = np.min(yy)
by = np.max(yy)
print(ax,ay,bx,by)
F = F.reshape((N+1,N+1))
xx = xx.reshape(N+1,N+1)
yy = yy.reshape(N+1,N+1)
Coeff_data = Coeff_data.reshape(N+1,N+1)

coeff = cheby_coeff_2d(F, N)

x_new = np.linspace(ax,bx,50)
y_new = np.linspace(bx,by,50)
xxnew,yynew = np.meshgrid(x_new,y_new)
Z_interp = cheby_interp_2d(xxnew, yynew ,coeff,N,ax,bx,ay,by)  # answer


# Plot the original function and the interpolated function
import matplotlib.pyplot as plt


fig, axs = plt.subplots(1, 2, figsize=(10, 4))

#print(Coeff_data)
#print()
#print(coeff)
#axs[0].contourf(X_new, Y_new, f(X_new,Y_new))
#axs[0].set_title("Original Function")
#print(Z_interp-F)
#print(F[:,0])
#axs[0].contourf(xx, yy, F)
#plt.contourf(xx, yy, np.log10(np.abs(Z_interp-F)))
#plt.imshow(np.log10(np.abs(Z_interp)))
#plt.imshow(Coeff_data-coeff)
#plt.title("Interpolated Function")
#plt.colorbar()
#plt.plot(yy,F,'bx')
#plt.plot(xx,(F[:,:]),'b-')
#plt.plot(xx,(np.log10(np.abs(Z_interp-F))),'r.')
#plt.xlim()
#plt.show()
#print(x_new)
#print(Z_interp)
#print(yy)
#print(Z_interp)
#print(F[:,0])
#print(f(X_new,Y_new))

#print(np.max(np.abs((Z_interp-f(X_new,Y_new)))))
