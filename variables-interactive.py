from bokeh.io import output_file, show, save
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure
import os
import numpy as np
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
showtrue = True

lw = 0.1 #0.04 for EDC-Taldice
alpha = 0.1 #0.04 for EDC-Taldice

##Global
variables=np.array([])
sigmas=np.array([])
D={}
DC={}


execfile('IceChronoClasses.py')

reader = emcee.backends.HDFBackend(datadir+"saved_iterations.hdf5")
flatchain = reader.get_chain(flat=True)

chains = flatchain[1::100000] #Eventually, we loop over variables, or use a class or function. flatchain should give us one variable set at a time

for dlabel in list_drillings: #This is the initialization loop
    D[dlabel] = Record(dlabel)
    D[dlabel].init() #When we do this we gain access to the prior as well...
    D[dlabel].model(D[dlabel].variables)

air_age = {}
depth = {}
ice_age = {}
methaneage = {}
methane = {}

for dlabel in list_drillings: #Makes a figure for each core, plots the priors
    air_age[dlabel] = [D[dlabel].airage]
    ice_age[dlabel] = [D[dlabel].age]
    depth[dlabel] = [D[dlabel].depth]
    methaneage[dlabel] = [D[dlabel].fct_age(D[dlabel].tuning_depth['CH4'])]
    methane[dlabel] = [D[dlabel].tuning_proxy['CH4']]

colors=["darkslateblue","saddlebrown"]


def fct_age(newdepth, depth, age):
    return np.interp(newdepth, depth, age)

for chain in chains:
    index = 0
    colorindex = 0
    for dlabel in list_drillings: #This is the extraction loop
        D[dlabel].variables=chain[index:index+np.size(D[dlabel].variables)]
        index=index+np.size(D[dlabel].variables)
        D[dlabel].model(D[dlabel].variables)
        depth[dlabel].append(D[dlabel].depth)
        air_age[dlabel].append(D[dlabel].airage)
        ice_age[dlabel].append(D[dlabel].age)
        methaneage[dlabel].append(fct_age(D[dlabel].tuning_depth['CH4'],D[dlabel].depth,D[dlabel].airage))
        methane[dlabel].append(D[dlabel].tuning_proxy['CH4'])

xpts = np.array([-.09, -.12, .0, .12,  .09])
ypts = np.array([-.1,   .02, .1, .02, -.1])
ypts2 = np.array([-3,-2,-1,0,1])

source = ColumnDataSource(dict(
        airageEDC=air_age['EDC'],
        iceageEDC=ice_age['EDC'],
        airageTALOS=air_age['TALDICE'],
        iceageTALOS=ice_age['TALDICE'],
        depthEDC=depth['EDC'],
        depthTALOS=depth['TALDICE'],
        methEDC=methane['EDC'],
        methTALOS=methane['TALDICE'],
        methageEDC = methaneage['EDC'],
        methageTALOS = methaneage['TALDICE'],
    )
)

TOOLS = "tap,box_zoom,pan,reset"

# create a new plot and add a renderer
iceEDC = figure(tools=TOOLS, plot_width=300, plot_height=300, title='EDC', x_axis_label='Ice Age (yr BP)', y_axis_label='Depth (m)')
iceEDC.multi_line(xs = 'iceageEDC', ys = 'depthEDC', source=source, line_color="darkslateblue", line_width=0.2)

# create another new plot and add a renderer
iceTALOS = figure(tools=TOOLS, plot_width=300, plot_height=300, title='TALDICE', x_axis_label='Ice Age (yr BP)', y_axis_label='Depth (m)')
iceTALOS.multi_line(xs = 'iceageTALOS', ys = 'depthTALOS', source=source, line_color="saddlebrown", line_width=0.2)

# create a new plot and add a renderer
airEDC = figure(tools=TOOLS, plot_width=300, plot_height=300, title='EDC', x_axis_label='Air Age (yr BP)', y_axis_label='Depth (m)')
airEDC.multi_line(xs = 'airageEDC', ys = 'depthEDC', source=source, line_color="darkslateblue", line_width=0.2)

# create another new plot and add a renderer
airTALOS = figure(tools=TOOLS, plot_width=300, plot_height=300, title='TALDICE', x_axis_label='Air Age (yr BP)', y_axis_label='Depth (m)')
airTALOS.multi_line(xs = 'airageTALOS', ys = 'depthTALOS', source=source, line_color="saddlebrown", line_width=0.2)

methane = figure(tools=TOOLS, plot_width=300, plot_height=300, title='CH4 Synchronization', x_axis_label='Air Age (yr BP)', y_axis_label='CH4 (ppbv)')
methane.multi_line(xs = 'methageEDC', ys = 'methEDC', source=source, line_color="darkslateblue", line_width=0.2)
methane.multi_line(xs = 'methageTALOS', ys = 'methTALOS', source=source, line_color="saddlebrown", line_width=0.2)


p = gridplot([[iceEDC, iceTALOS],[airEDC, airTALOS],[methane]])

output_file(os.getcwd() + '/bokehtest.html')
save(p)