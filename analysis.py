#!/usr/bin/env python
# -*- coding: utf-8 -*-

#This program accesses the posterior result of a dating experiment saved in the emcee-HDF5 framework.

import sys
#import time
#import math as m
import numpy as np
import matplotlib.pyplot as mpl
#import warnings
import os
#import scipy.linalg #
#from scipy.optimize import leastsq, minimize, basinhopping
import emcee #
#import scipy.sparse.linalg
#import gc
#import more_itertools as mit
#import copy
#import shutil
#import h5py #
#from sksparse.cholmod import cholesky as cholesky_sparse #
#import pickle
#from matplotlib.backends.backend_pdf import PdfPages
#from scipy.linalg import cholesky
import matplotlib.patches as mpatches
import shutil

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
showtrue = False

lw = 0.1 #0.04 for EDC-Taldice
alpha = 0.1 #0.04 for EDC-Taldice

##Global
variables=np.array([])
sigmas=np.array([])
D={}
DC={}

##Functions and Classes


execfile('IceChronoClasses.py')

reader = emcee.backends.HDFBackend(datadir+"saved_iterations.hdf5")
flatchain = reader.get_chain(flat=True,thin=10)
autocorr_time = reader.get_autocorr_time(discard=0,quiet=True)
print np.mean(autocorr_time)


if False:
    chains = flatchain #Eventually, we loop over variables, or use a class or function. flatchain should give us one variable set at a time

    for dlabel in list_drillings: #This is the initialization loop
        D[dlabel] = Record(dlabel)
        D[dlabel].init() #When we do this we gain access to the prior as well...
        D[dlabel].model(D[dlabel].variables)

    figdict = {}
    for dlabel in list_drillings: #Makes a figure for each core, plots the priors
        for name in ['air_age','ice_age','accu','LID','thinning']:
            figdict[dlabel+name] = mpl.subplots()

    colors=["darkslateblue","saddlebrown"]


    def air_age():
        for chain in chains:
            index = 0
            colorindex = 0
            for dlabel in list_drillings: #This is the extraction loop
                D[dlabel].variables=chain[index:index+np.size(D[dlabel].variables)]
                index=index+np.size(D[dlabel].variables)
                D[dlabel].model(D[dlabel].variables)
                figdict[dlabel+'air_age'][1].plot(D[dlabel].airage, D[dlabel].depth,linewidth=lw,alpha=alpha,color=colors[colorindex])
                figdict[dlabel+'air_age'][1].set_ylabel('Depth')
                figdict[dlabel+'air_age'][1].set_xlabel('Air Age (yr BP)')
                color_patch = mpatches.Patch(color=colors[colorindex],label=dlabel)
                color_patch2 = mpatches.Patch(color='gray', label=dlabel + ' prior')
                figdict[dlabel+'air_age'][1].legend(handles=[color_patch, color_patch2])
                colorindex += 1
                x1, x2, y1, y2 = figdict[dlabel+'air_age'][1].axis()
                figdict[dlabel+'air_age'][1].axis((x1, x2, 1100, 0))
                #mpl.savefig(os.getcwd() +'/'+datadir+ dlabel+'_air_age.png')

    def ice_age():
        for chain in chains:
            index = 0
            colorindex = 0
            for dlabel in list_drillings:  # This is the extraction loop
                D[dlabel].variables=chain[index:index+np.size(D[dlabel].variables)]
                index=index+np.size(D[dlabel].variables)
                D[dlabel].model(D[dlabel].variables)
                figdict[dlabel+'ice_age'][1].plot(D[dlabel].age, D[dlabel].depth,linewidth=lw,alpha=alpha,color=colors[colorindex])
                figdict[dlabel+'ice_age'][1].set_ylabel('Depth')
                figdict[dlabel+'ice_age'][1].set_xlabel('Ice Age (yr BP)')
                if (np.size(D[dlabel].icemarkers_depth) > 0):
                    figdict[dlabel+'ice_age'][1].errorbar(D[dlabel].icemarkers_age, D[dlabel].icemarkers_depth, color='red', xerr=D[dlabel].icemarkers_sigma,
                                 linestyle='', marker='o', markersize=2)
                colorindex += 1

        colorindex = 0
        for dlabel in list_drillings:#This loop makes the legends, axes etc.
            color_patch = mpatches.Patch(color=colors[colorindex],label=dlabel)
            color_patch2 = mpatches.Patch(color='gray', label=dlabel + ' prior')
            color_patch3 = mpatches.Patch(color='red', label="Dated horizons")
            figdict[dlabel+'ice_age'][1].legend(handles=[color_patch, color_patch2, color_patch3])
            colorindex += 1
            x1, x2, y1, y2 = figdict[dlabel+'ice_age'][1].axis()
            figdict[dlabel+'ice_age'][1].axis((D[dlabel].age_top, x2, D[dlabel].depth[-1], D[dlabel].depth[0]))
            #mpl.savefig(os.getcwd() +'/'+datadir+ dlabel + '_ice_age.png')

    def accu():
        for chain in chains:
            index = 0
            colorindex = 0
            for dlabel in list_drillings: #This is the extraction loop
                D[dlabel].variables=chain[index:index+np.size(D[dlabel].variables)]
                index=index+np.size(D[dlabel].variables)
                D[dlabel].model(D[dlabel].variables)
                figdict[dlabel+'accu'][1].step(D[dlabel].age, np.concatenate((D[dlabel].a, np.array([D[dlabel].a[-1]]))),linewidth=lw,alpha=alpha,color=colors[colorindex])
                figdict[dlabel+'accu'][1].set_xlabel('Age (yr BP)')
                figdict[dlabel+'accu'][1].set_ylabel('Accumulation (m/yr)')
                color_patch = mpatches.Patch(color=colors[colorindex],label=dlabel)
                color_patch2 = mpatches.Patch(color='gray', label=dlabel + ' prior')
                figdict[dlabel+'accu'][1].legend(handles=[color_patch, color_patch2])
                colorindex += 1
                #mpl.savefig(os.getcwd() +'/'+datadir+ dlabel + 'accu.png')

    def LID():
        for chain in chains:
            index = 0
            colorindex = 0
            for dlabel in list_drillings: #This is the extraction loop
                D[dlabel].variables=chain[index:index+np.size(D[dlabel].variables)]
                index=index+np.size(D[dlabel].variables)
                D[dlabel].model(D[dlabel].variables)
                figdict[dlabel+'LID'][1].step(D[dlabel].age, D[dlabel].LID,linewidth=lw,alpha=alpha,color=colors[colorindex])
                figdict[dlabel+'LID'][1].set_ylabel('LID (m)')
                figdict[dlabel+'LID'][1].set_xlabel('Age (yr BP)')
                color_patch = mpatches.Patch(color=colors[colorindex],label=dlabel)
                color_patch2 = mpatches.Patch(color='gray', label=dlabel + ' prior')
                figdict[dlabel+'LID'][1].legend(handles=[color_patch, color_patch2])
                colorindex += 1
                #mpl.savefig(os.getcwd() +'/'+datadir+ dlabel + 'LID.png')

    def thinning():
        for chain in chains:
            index = 0
            colorindex = 0
            for dlabel in list_drillings: #This is the extraction loop
                D[dlabel].variables=chain[index:index+np.size(D[dlabel].variables)]
                index=index+np.size(D[dlabel].variables)
                D[dlabel].model(D[dlabel].variables)
                figdict[dlabel+'thinning'][1].step(D[dlabel].tau,D[dlabel].depth_mid,linewidth=lw,alpha=alpha,color=colors[colorindex])
                figdict[dlabel+'thinning'][1].set_xlabel('Depth')
                figdict[dlabel+'thinning'][1].set_ylabel(r'Thinning Function $\tau$')
                color_patch = mpatches.Patch(color=colors[colorindex],label=dlabel)
                color_patch2 = mpatches.Patch(color='gray', label=dlabel+' prior')
                figdict[dlabel+'thinning'][1].legend(handles=[color_patch,color_patch2])
                colorindex += 1
                #mpl.savefig(os.getcwd() +'/'+datadir+ dlabel + 'thinning.png')


    for dlabel in list_drillings: #Plot the priors
        D[dlabel].init()  # Reinitialize
        D[dlabel].model(D[dlabel].variables)
        figdict[dlabel + 'air_age'][1].plot(D[dlabel].airage_model, D[dlabel].depth, color='gray')
        figdict[dlabel + 'ice_age'][1].plot(D[dlabel].age_model, D[dlabel].depth, color='gray')
        figdict[dlabel + 'accu'][1].step(D[dlabel].age,
                                             np.concatenate((D[dlabel].a_model, np.array([D[dlabel].a_model[-1]]))),
                                             color='gray')
        figdict[dlabel + 'LID'][1].step(D[dlabel].age, D[dlabel].LID_model, color='gray')
        figdict[dlabel + 'thinning'][1].step(D[dlabel].tau_model, D[dlabel].depth_mid, color='gray')
        if showtrue:
            shutil.move(os.getcwd()+'/'+datadir+'VK/thinning-prior.txt',os.getcwd()+'/'+datadir+'VK/thinning-prior-save.txt')
            shutil.move(os.getcwd()+'/'+datadir+'VK/thinning-prior-true.txt',os.getcwd()+'/'+datadir+'VK/thinning-prior.txt')
            D[dlabel].init()  # Reinitialize
            D[dlabel].model(D[dlabel].variables)
            figdict['VK' + 'ice_age'][1].plot(D[dlabel].age_model, D[dlabel].depth, color='black')
            shutil.move(os.getcwd()+'/'+datadir+'VK/thinning-prior.txt',os.getcwd()+'/'+datadir+'VK/thinning-prior-true.txt')
            shutil.move(os.getcwd()+'/'+datadir+'VK/thinning-prior-save.txt',os.getcwd()+'/'+datadir+'VK/thinning-prior.txt')

    air_age()
    #ice_age()
    #accu()
    #LID()
    #thinning()

    mpl.show()