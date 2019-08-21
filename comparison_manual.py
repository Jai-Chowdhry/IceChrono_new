#!/usr/bin/env python
# -*- coding: utf-8 -*-

#This program compares manual and automated synchronization

import sys
import time
import math as m
import numpy as np
import matplotlib.pyplot as mpl
import warnings
import os
import scipy.linalg #
from scipy.optimize import leastsq, minimize, basinhopping
import emcee #
import scipy.sparse.linalg
import gc
import more_itertools as mit
import copy
import shutil
import h5py #
from sksparse.cholmod import cholesky as cholesky_sparse #
import pickle
from matplotlib.backends.backend_pdf import PdfPages
from scipy.linalg import cholesky
import matplotlib.patches as mpatches

#!/usr/bin/env python
# -*- coding: utf-8 -*-

#This program accesses the posterior result of a dating experiment saved in the emcee-HDF5 framework.

import sys
import time
import math as m
import numpy as np
import matplotlib.pyplot as mpl
import warnings
import os
import scipy.linalg #
from scipy.optimize import leastsq, minimize, basinhopping
import emcee #
import scipy.sparse.linalg
import gc
import more_itertools as mit
import copy
import shutil
import h5py #
from sksparse.cholmod import cholesky as cholesky_sparse #
import pickle
from matplotlib.backends.backend_pdf import PdfPages
from scipy.linalg import cholesky
import matplotlib.patches as mpatches

###Reading parameters directory
datadir=sys.argv[1]
if datadir[-1]!='/':
    datadir=datadir+'/'
print 'Parameters directory is: ',datadir
#os.chdir(datadir)

##Parameters
scale_ageci=10.
show_figures=False
show_airlayerthick=False
execfile(datadir+'/parameters.py')
execfile('DEMove_modified.py') #Fixme maybe we can import this class

##Global
variables=np.array([])
sigmas=np.array([])
D={}
DC={}

##Functions and Classes


execfile('IceChronoClasses.py')

reader = emcee.backends.HDFBackend(datadir+"saved_iterations.hdf5")
flatchain = reader.get_chain(flat=True)
maxindex = np.argmin(reader.get_log_prob(flat=True))
#x = reader.get_autocorr_time(c=5,tol=50,thin=1,quiet=True)
#print x


chains = flatchain #Eventually, we loop over variables, or use a class or function. flatchain should give us one variable set at a time

for dlabel in list_drillings: #This is the initialization loop
    D[dlabel] = Record(dlabel)
    D[dlabel].init() #When we do this we gain access to the prior as well...
    D[dlabel].model(D[dlabel].variables)

#TODO: Now, we want to load manual tie points and uncertainties, and compare with the stats for the matching depth values of the MCMC iterations.

manual_points = np.loadtxt(os.getcwd()+'/'+datadir+'EDC-TALDICE/air_depth-loulergue.txt') #TODO Change according to experiment. Could this be an input??

def depth_equivalent(series,core1,core2,list):
    depth_1 = list[core1].depth
    airage_1 = list[core1].airage
    depth_2 = list[core2].depth
    airage_2 = list[core2].airage

    newage = np.interp(series,depth_1,airage_1)
    newdepth = np.interp(newage,airage_2,depth_2)

    return newdepth

equivalents = []

for chain in chains:
    index = 0
    for dlabel in list_drillings:  # This is the extraction loop
        D[dlabel].variables = chain[index:index + np.size(D[dlabel].variables)]
        index = index + np.size(D[dlabel].variables)
        D[dlabel].model(D[dlabel].variables)
    equivalents.append(depth_equivalent(manual_points[:,0],'EDC','TALDICE',D))

#Get stats from equivalents
equivalent_means = np.mean(equivalents, axis=0)
equivalent_std = np.std(equivalents, axis=0)

print(equivalent_means)
print(np.shape(manual_points[:,0]),np.shape(equivalent_means),np.shape(equivalent_std))


mpl.errorbar(manual_points[:,0],manual_points[:,1],xerr=manual_points[:,2],yerr=manual_points[:,2],color='red',label='Manual Tie Points (Loulergue et al., 2008)') #This is
mpl.errorbar(manual_points[:,0],equivalent_means,xerr=equivalent_std,yerr=equivalent_std,color='blue',label='Automated tie points (this study)')
mpl.xlabel('EDC depth (m)')
mpl.ylabel('TALDICE depth (m)')
mpl.legend()

new_points = np.vstack((manual_points.transpose(),equivalent_means))
new_points = np.vstack((new_points,equivalent_std))
np.savetxt(datadir+'/comparison-synchro.txt',new_points.transpose(),fmt='%0.2f')

#Optional: plot best synchro. Toggle if to True/False

'''best = flatchain[maxindex]
index = 0
for dlabel in list_drillings:  # This is the extraction loop
    D[dlabel].variables = chain[index:index + np.size(D[dlabel].variables)]
    index = index + np.size(D[dlabel].variables)
    D[dlabel].model(D[dlabel].variables)
    x = D[dlabel].residuals(D[dlabel].variables)

if True:
    for dlabel in list_drillings:
        for proxy, tag in D[dlabel].tuning_dict.items():
            fig = mpl.figure(D[dlabel].label + ' ' + proxy)
            ax = mpl.subplot(111)
            # Age at final depth, proxy value minus sigma, proxy value plus sigma
            if D[dlabel].tuning_multi:
                tuning_label = ''
                for drilling in list_drillings:
                    if D[drilling].label != D[dlabel].label:
                        tuning_label += D[drilling].label + ' '

                ax.plot(D[dlabel].tuning_age[proxy], D[dlabel].tuning_target[proxy], color="#3F5D7D", linewidth=0.3, label=tuning_label)

            else:
                ax.plot(D[dlabel].tuning_age[proxy], D[dlabel].tuning_target[proxy], color="#3F5D7D", linewidth=0.3,
                        label='Target')
            if proxy in D[dlabel].tuning_target_sigma:
                ax.fill_between(D[dlabel].tuning_age[proxy],
                                D[dlabel].tuning_target[proxy] - 2 * D[dlabel].tuning_target_sigma[proxy],
                                D[dlabel].tuning_target[proxy] + 2 * D[dlabel].tuning_target_sigma[proxy], color="#3F5D7D", alpha=0.5,
                                linewidth=0.0)
            leg = mpl.legend(frameon=False, loc="best")
            if tag == 'ice':

                ax.plot(D[dlabel].fct_age(D[dlabel].tuning_depth[proxy]), D[dlabel].tuning_proxy[proxy], linewidth=0.3,
                        color=(20 / 255., 20 / 255., 20 / 255.), label=D[dlabel].label)
                ax.fill_between(D[dlabel].fct_age(D[dlabel].tuning_depth[proxy]),
                                D[dlabel].tuning_proxy[proxy] - 2 * D[dlabel].tuning_proxy_sigma[proxy],
                                D[dlabel].tuning_proxy[proxy] + 2 * D[dlabel].tuning_proxy_sigma[proxy],
                                facecolor=(140 / 255., 86 / 255., 75 / 255.), linewidth=0.0, alpha=0.5)
            elif tag == 'air':

                ax.plot(D[dlabel].fct_airage(D[dlabel].tuning_depth[proxy]), D[dlabel].tuning_proxy[proxy], linewidth=0.3,
                        color=(20 / 255., 20 / 255., 20 / 255.), label=D[dlabel].label)
                ax.fill_between(D[dlabel].fct_airage(D[dlabel].tuning_depth[proxy]),
                                D[dlabel].tuning_proxy[proxy] - 2 * D[dlabel].tuning_proxy_sigma[proxy],
                                D[dlabel].tuning_proxy[proxy] + 2 * D[dlabel].tuning_proxy_sigma[proxy],
                                facecolor=(140 / 255., 86 / 255., 75 / 255.), linewidth=0.0, alpha=0.5)
            leg = mpl.legend(frameon=False, loc="best")

            ax.tick_params(axis='x', top='off')
            ax.tick_params(axis='y', right='off')
            mpl.xlabel('Age (yr BP)')
            if hasattr(D[dlabel], 'tuning_units'):
                mpl.ylabel(proxy + ' (' + D[dlabel].tuning_units[proxy] + ')')
            else:
                mpl.ylabel(proxy)

            pp = PdfPages(datadir + D[dlabel].label + '/tuning' + proxy + '.pdf')
            pp.savefig(mpl.figure(D[dlabel].label + ' ' + proxy), bbox_inches='tight')
            pp.close()
            if not show_figures:
                mpl.close()'''

mpl.show()
