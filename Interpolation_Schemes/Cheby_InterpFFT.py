import numpy as np
import pyfftw.interfaces.numpy_fft as fft
from scipy.fftpack import dct

def cheby_exgrid(N,a,b):
    i=np.linspace(0,N,N+1)
    return (b-a)/2.*(-np.cos(i*np.pi/N))+(b+a)/2.0

def cheb_interp_1d(y, N):
    # Compute Chebyshev nodes
    # x = -np.cos(np.pi * np.arange(N) / (N - 1))

    # # Compute function values at Chebyshev nodes
    # y = f(x)

    # Compute DCT of function values
    c = dct(y, type=1)

    # Scale coefficients to obtain Chebyshev polynomial coefficients
    a = np.zeros(N+1)
    a = c/(N)

    return a

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


def cheby_interp(x,coeff,N,a,b):

    x_chebgrid = -( 2.0 * x - a  - b ) / ( b-a);
    sum  = coeff[0]*cheby_poly(x_chebgrid,0) * 0.5    
    for j in range(1,N):
        sum += coeff[j]*cheby_poly(x_chebgrid,j)
    sum += coeff[N]*cheby_poly(x_chebgrid,N) * 0.5
    return sum

N = 16
a = -1
b = 1
# Example usage
#f = lambda x: np.exp(x)
f = lambda x: np.sin(3.*x)
#x = -np.cos(np.pi * np.arange(N) / (N - 1))
x = cheby_exgrid(N,a,b)
y = f(x)
coeff = cheb_interp_1d(y, N)
#print(a)

x = np.linspace(a, b, 40)

interp =  cheby_interp(x,coeff,N-1,a,b)

print(f(x))
print(interp)
import matplotlib.pyplot as plt
# #plt.plot(x, y, label='Chebyshev interpolant')
#plt.plot((x), f(x), '--', label='True function')
#plt.plot((x), interp, '.', label='my interp')
plt.plot((x), np.abs(f(x)-interp), '-', label='my interp')

# plt.legend()
plt.show()
