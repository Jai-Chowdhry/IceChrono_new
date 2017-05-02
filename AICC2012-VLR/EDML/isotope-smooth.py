import sys
import numpy as np
import os
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as mpl

filename = sys.argv[1]
smooth = sys.argv[2] #Smoothing parameter
isoarray = np.loadtxt(os.getcwd()+'/'+filename)
spliso = UnivariateSpline(isoarray[:,0],isoarray[:,1],s=smooth)
splsigma = UnivariateSpline(isoarray[:,0],isoarray[:,2])

saveiso = np.vstack((np.vstack((np.array(isoarray[:,0]),np.array(spliso(isoarray[:,0])))),np.array(splsigma(isoarray[:,0]))))
np.savetxt(os.getcwd()+'/smooth_'+filename,np.transpose(saveiso))

mpl.plot(isoarray[:,0],isoarray[:,1])
mpl.plot(isoarray[:,0],spliso(isoarray[:,0]))
mpl.show()
