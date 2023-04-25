import numpy as np
import matplotlib.pyplot as plt

# make multidomain with equal length and N

import math

def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy

def makeCheb_oni(a,b,N,ith):
    N=np.int64(N)
    i=np.linspace(0,N,N+1)

    scale = (a-b) / (np.cos((N-ith)*np.pi/N) - np.cos(ith*np.pi/N))
    mid   = b - scale * (-np.cos((N-ith)*np.pi/N)) 
    print("scale=",scale)
#    return (b-a)/2.*(-np.cos(i*np.pi/N))+(b-a)/2.0
#    return (b-a+scale)/2.*(-np.cos(i*np.pi/N))+(b-a)/2.0
    return scale *(-np.cos(i*np.pi/N)) + mid

def makeX(N,L):
    N=np.int64(N)
    i=np.linspace(0,N,N+1)
    return L/2.*(-np.cos(i*np.pi/N))+L/2.0

def MakeL(a,b,N,NsubD):
    return (b-a)/(NsubD+(1-NsubD)*(1.0-np.cos(np.pi/N))/2.)

a = 0
b = 1
L = b-a
N = 16
N2 = 20
Ng = 8
NsubD = 4
NN = N + (N-1)*(NsubD-1)

#L = MakeL(a,b,N,NsubD)
x  = makeCheb_oni(0,a+L/2.,N+2*Ng,Ng)
x2 = makeCheb_oni(a+L/2.,a+L,N+2*Ng,Ng)
#x3 = makeCheb_oni(1,1.5,N,4)
#x  = makeX(N,L)
#x4 = x + 2 * (L-1.*(x[1]-x[0]))

xx, yy = np.meshgrid(x,x)
xx2, yy2 = np.meshgrid(x2,x2)

# for h-adaptivity 
# #L2 = MakeL(x[-2],x4[1],N2,2)
# x2 = makeX(N2,L2) + 1 * (L-1.*(x[1]-x[0]))
# x3 = makeX(N2,L2) + x2[-2]

# print(L,NN)
# print(x)
# print(x2)
# #print(x3)
# #print(x4)

# plt.plot(x,np.ones(N+1),'.',
#          x2,np.ones(N2+1),'x',
#          x3,np.ones(N2+1),'.',
#          x4,np.ones(N+1),'x')
# plt.show()

# for p-adaptivity 
# L2 = L
# x2 = makeX(N2,L2) + 1 * (L-1.*(x[1]-x[0]))
# x3 = makeX(N2,L2) + x2[-2]

#print(L,NN)
#print(x)
#print(x2)
#print(x3)
#print(x4)

#plt.plot(x,np.ones(N+1),'.',
#         x2,np.ones(N+1),'x')
#         x4,np.ones(N+1),'x')
plt.plot(xx,yy,'r.',
#         xx+1.0,yy,'bx'
#         xx,yy+1.0,'cx')
         xx2,yy2,'gx',)
#plt.plot(x,x,'kx')
x2 =  makeCheb_oni(0,1,N,8)

# xa = np.abs(x2[0])
# x2[:] += a
# x2[:] *= np.sqrt(2)
# x2[:] -= a

# print(x2[0],x[0])
# print(x2[-1])
# print(np.sqrt((x[-1]-x[0])**2+(x[-1]-x[0])**2))

# x3=[]
# x3.append(x[0])
# for i in range(1,x.shape[0]):
#     x3.append(x3[i-1] + np.sqrt(2*(x[i]-x[i-1])**2))
# x3 = np.array(x3)
# print(x)
# print(x2)
# print(x3)
# print(x3[-1]-x3[0])
# #plt.plot(x2,yy[0],'k.')
# plt.plot(x3,np.ones(x.shape[0])*-2,'k.')

#plt.plot(x2,np.ones(x.shape[0])*-2,'k.')
#plt.plot(x2,np.ones(x.shape[0]),'kx')
#print(rotate((x3[0],x[0]),(x3[-1],x[0]),math.radians(np.pi/4)))
plt.show()
