import numpy as np

def RMS(y,y_true):
    N = y.shape[0]
    return np.sqrt(np.sum((y-y_true)*(y-y_true))/N)

data = np.loadtxt('Gaussian_Dens_000001_N128',usecols=[0,3,4])
print(np.max(np.abs(data[:,2]-data[:,1])))
print(RMS((data[:,1]),(data[:,2])))
#print(RMS(np.sqrt(data[:,1]),np.sqrt(data[:,2]) ) )
