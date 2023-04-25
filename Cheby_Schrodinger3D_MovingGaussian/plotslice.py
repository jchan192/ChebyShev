import numpy as np
import matplotlib.pyplot as plt

z = 16

data = np.loadtxt("output/XYZ600.dat")
slice = data[z*17**2:(z+1)*17**2,4].reshape(17,17)
plt.imshow(slice)
plt.show()
