import numpy as np
import matplotlib.pyplot as plt

def cheby_exgrid(N,a,b):
    i=np.linspace(0,N,N+1)
    return (b-a)/2.*(-np.cos(i*np.pi/N))+(b+a)/2.0

def cheby_polylist(x,n):
    T = []
    for i in n:
        T.append(cheby_poly(x,i))

    T = np.array(T)
    return T

def cheby_poly(x,n):
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

def cheby_coeff(x,y,N):
    a = []
    y[0] *= 0.5
    y[N] *= 0.5

    x_chebgrid = cheby_exgrid(N,-1.0,1.0)
#    print(x,x_chebgrid)

    for j in range(N+1):
        sum = 2.0 / N  * np.sum(y*cheby_poly(x_chebgrid,j))
        a.append(sum)

    y[0] /= 0.5
    y[N] /= 0.5
    return a

def cheby_interp(x,coeff,N,a,b):

    x_chebgrid = ( 2.0 * x - a  - b ) / ( b-a);
    sum  = coeff[0]*cheby_poly(x_chebgrid,0) * 0.5    
    for j in range(1,N-1):
        sum += coeff[j]*cheby_poly(x_chebgrid,j)
    sum += coeff[N]*cheby_poly(x_chebgrid,N) * 0.5
    return sum

def cheby_coeff2D(x,y,f,N):
    a = []
    f[:,0] *= 0.5
    f[:,N] *= 0.5
    f[0,:] *= 0.5
    f[N,:] *= 0.5

    x_chebgrid = cheby_exgrid(N,-1.0,1.0)
    xx,yy = np.meshgrid(x_chebgrid,x_chebgrid)
#    print(x,x_chebgrid)

    for i in range(N+1):
        for j in range(N+1):
            sum = 2.0 / N  * np.sum(f*cheby_poly(xx,i)*cheby_poly(yy,j))
            a.append(sum)

    f[:,0] /= 0.5
    f[:,N] /= 0.5
    f[0,:] /= 0.5
    f[N,:] /= 0.5

    a = np.array(a).reshape(N+1,N+1)
    return a

def cheby_interp2D(x,y,coeff,N,a,b):

    # problem !!!! Not finished. how to do chebgridlist

    x_chebgrid = ( 2.0 * x - a  - b ) / ( b-a);
    y_chebgrid = ( 2.0 * y - a  - b ) / ( b-a);

    print(cheby_polylist(y_chebgrid,np.int64(np.arange(N+1))).shape)

    sum = 0
    for j in range(0,N+1):
            sum += coeff[0,j]*cheby_poly(x_chebgrid,0)*cheby_poly(y_chebgrid,j) * 0.5

    for i in range(0,N+1):
            sum += coeff[i,0]*cheby_poly(x_chebgrid,i)*cheby_poly(y_chebgrid,0) * 0.5

    for i in range(1,N-1):
        for j in range(1,N-1):
            sum += coeff[i,j]*cheby_poly(x_chebgrid,i)*cheby_poly(y_chebgrid,j)

    for j in range(0,N+1):
            sum += coeff[N,j]*cheby_poly(x_chebgrid,N)*cheby_poly(y_chebgrid,j) * 0.5

    for i in range(0,N+1):
            sum += coeff[i,N]*cheby_poly(x_chebgrid,i)*cheby_poly(y_chebgrid,N) * 0.5
    return sum

def Make_PolyMatrix(x,y,N,a,b):

    x_chebgrid = ( 2.0 * x - a  - b ) / ( b-a);
    y_chebgrid = ( 2.0 * y - a  - b ) / ( b-a);

    M = np.zeros([N+1,N+1])

    for i in range(N+1):
        for j in range(N+1):
            M[i,j] = cheby_poly(x_chebgrid,i)*cheby_poly(y_chebgrid,j)

    return M

# test it with sin(x)
N = 5
a = -1
b = 1
x = cheby_exgrid(N,a,b)
xx,yy = np.meshgrid(x,x)

print(Make_PolyMatrix(0.8,0.8,N,a,b))
print(Make_PolyMatrix(-0.8,-0.8,N,a,b))
print(Make_PolyMatrix(-0.8,0.8,N,a,b))
print(Make_PolyMatrix(0.8,-0.8,N,a,b))

exit(0);
f = np.sin(5.*xx) + np.sin(8.*yy)
xint1D = np.linspace(a,b,24+1)

coeff = cheby_coeff(yy[:,np.int64(N/2)],f[:,np.int64(N/2)],N)
interp1D =  cheby_interp(xint1D,coeff,N,a,b)

print(xint1D.shape,yy[:,np.int64(N/2)].shape)
coeff2D = cheby_coeff2D(xx,yy,f,N)
interp2D =  cheby_interp2D(xint1D[:],yy[:,np.int64(N/2)],coeff2D,N,a,b)

#plt.plot(x,f,'-',xx,interp,'.')
#plt.plot(x,f,'-',xx,interp,'.')
plt.plot(yy[:,np.int64(N/2)],f[:,np.int64(N/2)],'-',
         yy[:,12],interp2D,'.')
plt.plot()
plt.show()

