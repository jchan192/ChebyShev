import numpy as np
from scipy.fftpack import dctn
import matplotlib.pyplot as plt

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

def cheby_coeff_3d(F, N):
    # Compute Chebyshev nodes
    # x = -np.cos(np.pi * np.arange(N) / (N - 1))
    # y = -np.cos(np.pi * np.arange(N) / (N - 1))
    # z = -np.cos(np.pi * np.arange(N) / (N - 1))

    # # Compute function values at Chebyshev nodes
    # X, Y, Z = np.meshgrid(x, y, z)
    # W = f(X, Y, Z)

    # Compute 3D DCT of function values
    c = dctn(F, type=1)

    # Scale coefficients to obtain Chebyshev polynomial coefficients
    a = np.zeros((N, N, N))
    a = c / ((N)**3)

    # still need to figure out the scaling .....
    a[0,:,:] *= 0.5
    a[:,0,:] *= 0.5
    a[:,:,0] *= 0.5
    a[N,:,:] *= 0.5
    a[:,N,:] *= 0.5
    a[:,:,N] *= 0.5
    #    a[:,N,1:N+1] *= 0.5
#    a[:,1:N+1,N] *= 0.5

    return a

def cheby_interp_3d(x,y,z,coeff,N,ax,bx,ay,by,az,bz):
    x_chebgrid = -( 2.0 * x - ax  - bx ) / ( bx-ax);
    y_chebgrid = -( 2.0 * y - ay  - by ) / ( by-ay);
    z_chebgrid = -( 2.0 * z - az  - bz ) / ( bz-az);

    sum = 0
    for i in range(0,N+1):
        for j in range(0,N+1):
            for k in range(0,N+1):
                sum += coeff[j,i,k]*cheby_poly(x_chebgrid,i)*cheby_poly(y_chebgrid,j)*cheby_poly(z_chebgrid,k)

    return sum

N = 16
ax = -2
bx = 1
ay = -1
by = 2
az = 0
bz = 1
x = cheby_exgrid(N,ax,bx)
y = cheby_exgrid(N,ay,by)
z = cheby_exgrid(N,az,bz)
xxx, yyy,zzz = np.meshgrid(x,y,z)

# example usage
f = lambda x, y, z: np.exp(x + y**2 + z**4)
f = lambda x, y, z: np.sin(1*x) + np.sin(2*y) + np.sin(4.*z) 
a = cheby_coeff_3d(f(xxx,yyy,zzz), N)
x_new = np.linspace(ax, bx, 10)
y_new = np.linspace(ay, by, 10)
z_new = np.linspace(az, bz, 10)
X_new, Y_new, Z_new = np.meshgrid(x_new, y_new, z_new)
F= f(X_new,Y_new,Z_new)
W_interp = cheby_interp_3d(X_new, Y_new, Z_new, a, N, ax,bx,ay,by,az,bz)


#print(W_interp-f(x_new,y_new,z_new))
fig, axs = plt.subplots(2, 3, figsize=(10, 4))

# ploting 
print(np.max(np.abs(W_interp[:,:,:]-F[:,:,:])))

#contourf(X_new[:,:,5],Y_new[:,:,5],F[:,:,5])
axs[0,0].imshow(F[5,:,:])
axs[0,1].imshow(F[:,5,:])
axs[0,2].imshow(F[:,:,5])
axs[1,0].imshow(W_interp[5,:,:])
axs[1,1].imshow(W_interp[:,5,:])
axs[1,2].imshow(W_interp[:,:,5])

#axs[1].imshow(W_interp[5,:,:])
#plt.imshow(np.log10(np.abs(W_interp[:,:,0]-F[:,:,0])))
#plt.imshow(W_interp[:,:,0]-F[:,:,0])
#plt.colorbar()
plt.show()
