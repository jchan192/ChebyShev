import time
import numpy as np
import pyfftw.interfaces.numpy_fft as fft
from scipy.fftpack import dct,dctn

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

s = time.time()
print(cheby_poly(-0.7,100))
print(time.time()-s)

def cheby_coeff_2d(Z, N):
    # Compute Chebyshev nodes
#    x = -np.cos(np.pi * np.linspace(0,N,N+1) / N )
#    y = -np.cos(np.pi * np.linspace(0,N,N+1) / N )

    # Compute function values at Chebyshev nodes
#    X, Y = np.meshgrid(x, y)
#    Z = f(X, Y)
#    print(Z.shape)
    # Compute 2D DCT of function values
    c = dctn(Z, type=1)
    
    # Scale coefficients to obtain Chebyshev polynomial coefficients
    a = np.zeros((N+1, N+1))
    a = c / (N**2)

    a[0,:] *= 0.5
    a[:,0] *= 0.5
    a[N,:] *= 0.5
    a[:,N] *= 0.5

#    print(a)
    return a

def cheby_interp_2d(x,y,coeff,N,ax,bx,ay,by):

    x_chebgrid = -( 2.0 * x - ax  - bx ) / ( bx-ax);
    y_chebgrid = -( 2.0 * y - ay  - by ) / ( by-ay);

    sum = 0
    for i in range(0,N+1):
        for j in range(0,N+1):
                sum += coeff[i,j]*cheby_poly(x_chebgrid,j)*cheby_poly(y_chebgrid,i)
    return sum

#f = lambda x, y: np.exp(x + y)
#f = lambda x, y: (x**4 + y**4)
f = lambda x, y: np.sin(5*x) + np.sin(5*y)
#f = lambda x, y: np.exp(-x**2-y**2)
ax = 1.5
bx = 3.5
ay = 0
by = 1
N = 20

x = cheby_exgrid(N,ax,bx)
y = cheby_exgrid(N,ay,by)
xx, yy = np.meshgrid(x,y)
F = f(xx,yy)
coeff = cheby_coeff_2d(F, N)

x_new = np.linspace(ax, bx, 200)
y_new = np.linspace(ay, by, 200)
X_new, Y_new = np.meshgrid(x_new, y_new)
Z_interp = cheby_interp_2d(X_new, Y_new,coeff,N,ax,bx,ay,by)  # answer

#print(cheby_interp_2d(x[1],x[1],coeff,N,ax,bx,ay,by))
# Plot the original function and the interpolated function
import matplotlib.pyplot as plt


#fig, axs = plt.subplots(1, 2, figsize=(10, 4))

#axs[0].contourf(X_new, Y_new, f(X_new,Y_new))
#axs[0].set_title("Original Function")

print(np.max(np.abs((Z_interp-f(X_new,Y_new)))))
#print(np.min(f(X_new,Y_new)))
#plt.contourf(X_new, Y_new, np.log10(np.abs(Z_interp-f(X_new,Y_new))))
#plt.imshow(np.log10(np.abs(Z_interp-f(X_new,Y_new))))
#plt.title("Interpolated Function")
#plt.colorbar()
#plt.plot(Y_new,f(X_new,Y_new),'bx')
#plt.plot((Y_new),Z_interp,'.')
#plt.plot((Y_new),np.abs(Z_interp-f(X_new,Y_new)),'.')
#plt.show()

#print(Z_interp)
#print(f(X_new,Y_new))

