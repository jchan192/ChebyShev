import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.interpolate import BSpline
from scipy import interpolate

def MycubicInterpolate(p,xp, x) :
    ax = xp[0]
    bx = xp[-1]
    dx = xp[1]-xp[0]
    idx = np.int64(np.floor(x/dx))
    x = (x-ax)/(bx-ax)

    print(x)
    print(idx)
    size = np.shape(x)[0]
    interp = []
    for i in range(size):
        interp.append( p[idx[i]] + 0.5 * x[i]*(p[idx[i]+1] - p[idx[i]-1] 
                                 + x[i]*(2.0*p[idx[i]-1] - 5.0*p[idx[i]] + 4.0*p[idx[i]+1] - p[idx[i]+2] 
                                 + x[i]*(3.0*(p[idx[i]] - p[idx[i]+1]) + p[idx[i]+2] - p[idx[i]-1]))))

    return interp

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

def function(x):

    return np.sin(2.*x)
    #return np.exp(x)
#f = np.loadtxt('test2')[:,[0,1]]
# N = np.shape(f)[0]-1

# print(f)
# xx     = cheby_exgrid(40,f[0,0],f[-1,0])
# coeff  = cheby_coeff(f[:,0],f[:,1],N)
# interp = cheby_interp(xx,coeff,N,f[0,0],f[-1,0])

# print( cheby_interp(0.247551,coeff,N,f[0,0],f[-1,0]))
# plt.plot(f[:,0],f[:,1],'>',xx,interp,'.')
# plt.show()

# test it with sin(x)
N = 16
a = 0
b = 2
#x = cheby_exgrid(N,a,b)
x = np.linspace(a,b,N)
f = function(x)
xx = np.linspace(x[1],x[-3],10)
finex = np.linspace(a,b,200)
finef = function(finex)

#coeff = cheby_coeff(x,f,N)
#print(coeff)

#cheb_interp =  cheby_interp(xx,coeff,N,a,b)
lin_interp = np.interp(xx,x,f)

spl = CubicSpline(x,f)
cspline_interp = spl(xx)

tck = interpolate.splrep(x, f, s=0, k=3) 
BSpline_interp = interpolate.BSpline(*tck)(xx)

mycubicspline_interp = MycubicInterpolate(f,x,xx)
print(xx)
#plt.plot(finex,finef,'-',xx,cheb_interp,'.',xx,lin_interp,'.',xx,cspline_interp,'.')
#plt.plot(xx,np.log10(np.abs(cheb_interp-function(xx))),'.',label="cheby")
plt.plot(xx,np.log10(np.abs(lin_interp-function(xx))),'.',label="linear")
plt.plot(xx,np.log10(np.abs(cspline_interp-function(xx))),'.',label="cubic spline")
plt.plot(xx,np.log10(np.abs(BSpline_interp-function(xx))),'x',label="B spline")
plt.plot(xx,np.log10(np.abs(mycubicspline_interp-function(xx))),'o',label="my cubicspline")

# plt.plot(x,f,'k-')
# plt.plot(x,f,'k.')
# plt.plot(xx,lin_interp,'o',label="lin")
# plt.plot(xx,cspline_interp,'x',label="cspline")
# plt.plot(xx,BSpline_interp,'.',label="Bspline")
# plt.plot(xx,mycubicspline_interp,'.',label="mycubicspline")

plt.plot()
plt.legend()
#plt.plot(x,f,'-',xx,interp,'.')
plt.show()

