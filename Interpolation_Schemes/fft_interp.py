import numpy as np
import matplotlib.pyplot as plt
# Define the periodic function (in this case, a sine wave)
def f(x):
    return np.sin(2*np.pi*x)

def fft_interp(coeffs, x):
    # 1D FFT interpolate for arbitrary x
    # assume x is periodic on [0,1] interval
    size = len(coeffs)
    kn = np.fft.fftfreq(size)
    eikx = np.exp(2.j*np.pi*x*size*kn)
    return np.dot(coeffs, eikx) / size

def fft_interp_vec(coeffs, xv):
    size = len(coeffs)
    kn = np.fft.fftfreq(size)
    eikx = np.exp( 2.j*np.pi*size*np.outer(xv, kn) )
    return np.einsum('ab,b->a', eikx, pix_fft) / size

# Define the range of x values over which the function is defined
xmin = 0
xmax = 1.

# Define the number of points used to sample the function
n_points = 16

# Generate the x values at which to sample the function
x = np.linspace(xmin, xmax, n_points, endpoint=False)
dx = x[1]-x[0]

# Evaluate the function at the sample points
y = f(x)

# Compute the FFT of the function
fft_y = np.fft.fft(y)


xx = np.linspace(0,1,100)
interp = []
for xxx in xx:
    interp.append(np.real(fft_interp(fft_y,xxx)))

# Print the interpolated value
#print("Interpolated value at x=1.5:", np.real(interp_val))

# Print the interpolated value
#plt.plot(x,f(x),'-')
#plt.plot(x[0]+dx*interp_index,interp_val,'.')
plt.plot(xx,np.log10(np.abs(interp-f(xx))),'.')
plt.show()
