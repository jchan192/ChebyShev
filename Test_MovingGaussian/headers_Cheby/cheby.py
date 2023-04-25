import numpy as np

'''Chebushev polynomial differentiation matrix.
Ref.: Trefethen's 'Spectral Methods in MATLAB' book.
'''
def cheb(N):
    N = 2
    x = np.cos(np.pi*np.linspace(N,0,N+1)/N)
    #    x[N/2]=0.0 # only when N is even!
    c=np.zeros(N+1)
    c[0]=2.
    c[1:N]=1.
    c[N]=2.
    c = c * (-1)**np.linspace(0,N,N+1)
    X = np.tile(x, (N+1,1))
    #print(X)
    dX = X - X.T
    #print(dX)
    D = np.dot(c.reshape(N+1,1),(1./c).reshape(1,N+1))
    #print(D)
    #print()
    #print(dX+np.eye(N+1))
    D = D / (dX+np.eye(N+1))    
    D = D - np.diag( D.T.sum(axis=0) )
    return D
