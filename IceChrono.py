#TODO: what about symbolic links in github?

import sys
import time
import math as m
import numpy as np
import matplotlib.pyplot as mpl
import multiprocessing
import warnings
import os
import scipy.linalg
from scipy.interpolate import interp1d
from scipy.optimize import leastsq, minimize, basinhopping
from matplotlib.backends.backend_pdf import PdfPages
from scipy.linalg import cholesky
from scipy.misc import derivative


###Registration of start time
start_time = time.clock()      #Use time.clock() for processor time

###Reading parameters directory
datadir=sys.argv[1]
if datadir[-1]!='/':
    datadir=datadir+'/'
print 'Parameters directory is: ',datadir
#os.chdir(datadir)

###Opening of output.txt file
output_file = open(datadir+'output.txt','a')

##Parameters
scale_ageci=10.
show_figures=False
show_airlayerthick=False
execfile(datadir+'/parameters.py')


##Global
variables=np.array([])
sigmas=np.array([])
D={}
DC={}

##Functions and Classes


execfile('IceChronoClasses.py')


def residuals(var):
    """Calculate the residuals."""
    resi=np.array([])
    index=0
    for i,dlabel in enumerate(list_drillings):
        D[dlabel].variables=var[index:index+np.size(D[dlabel].variables)]
        index=index+np.size(D[dlabel].variables)
        resi=np.concatenate((resi,D[dlabel].residuals(D[dlabel].variables)))
        for j,dlabel2 in enumerate(list_drillings):
            if j<i:
                resi=np.concatenate((resi,DC[dlabel2+'-'+dlabel].residuals()))
    return resi

def cost_function(var):
    cost=np.dot(residuals(var),np.transpose(residuals(var)))
    return cost


def Dres(var):
    """Calculate derivatives for each parameter using pool."""
    zeropred = residuals(var)
    derivparams = []
    results=[]
    delta = m.sqrt(np.finfo(float).eps) #Stolen from the leastsq code
    for i in range(len(var)): #fixme: This loop is probably sub-optimal. Have a look at what does leastsq to improve this.
        copy = np.array(var)
        copy[i] += delta
        derivparams.append(copy)
#        results.append(residuals(derivparams))
    if __name__ == "__main__":
        pool = multiprocessing.Pool(nb_nodes)
    results = pool.map(residuals, derivparams)
    pool.close() # Patch to avoid memory allocation errors
    derivs = [ (r - zeropred)/delta for r in results ]
    return derivs

##MAIN


##Initialisation
for i,dlabel in enumerate(list_drillings):

    print 'Initialization of drilling '+dlabel
        
    D[dlabel]=Drilling(dlabel)
    D[dlabel].init()
    D[dlabel].model(D[dlabel].variables)
#    D[dlabel].a_init=D[dlabel].a
#    D[dlabel].LID_init=D[dlabel].LID
    D[dlabel].write_init()
#    D[dlabel].display_init()
    variables=np.concatenate((variables,D[dlabel].variables))
    sigmas=np.concatenate((sigmas,D[dlabel].sigmap_all))

for i,dlabel in enumerate(list_drillings):
    for j,dlabel2 in enumerate(list_drillings):
        if j<i:
            print 'Initialization of drilling pair '+dlabel2+'-'+dlabel
            DC[dlabel2+'-'+dlabel]=DrillingPair(D[dlabel2],D[dlabel])
            DC[dlabel2+'-'+dlabel].init()
#            DC[dlabel2+'-'+dlabel].display_init()

##Optimization
start_time_opt = time.time()
print 'cost function: ',cost_function(variables)
if opt_method=='leastsq':
    print 'Optimization by leastsq'
    variables,hess,infodict,mesg,ier=leastsq(residuals, variables, full_output=1)
elif opt_method=='leastsq-parallel':
    print 'Optimization by leastsq-parallel'
    variables,hess,infodict,mesg,ier=leastsq(residuals, variables, Dfun=Dres, col_deriv=1, full_output=1)
elif opt_method=="L-BFGS-B":
    print 'Optimization by L-BFGS-B'
    res=minimize(cost_function, variables, method='L-BFGS-B', jac=False)
    variables=res.x
    print 'number of iterations: ',res.nit
    hess=np.zeros((np.size(variables),np.size(variables)))
    print 'Message: ',res.message
#    cost=cost_function(variables)
elif opt_method=="dogleg":
    print 'Optimization by dogleg trust regions'
    res=minimize(cost_function,variables, method='dogleg', jac=Dres)
    variables = res.x
    print 'number of iterations: ', res.nit
    hess = np.zeros((np.size(variables), np.size(variables)))
    print 'Message: ', res.message

elif opt_method=='basinhopping':
    print 'Optimization by basin hopping'
    leak_f = []
    leak_hess = [np.zeros((np.size(variables),np.size(variables)))]

    class StepTaker(object):
        def __init__(self,  stepsize = 0.01): #0.01
            self.stepsize=stepsize

        def __call__(self, stratopause):
            s = self.stepsize  # Allows the distribution to be widened or narrowed.
            uni = np.random.multivariate_normal(np.zeros(np.size(variables)),np.diag(np.ones(np.size(variables))))
            stratopause += uni * s
            return stratopause

    # class CallbackFun(object):
    #     def __init__(self):
    #         self.callbackiter = 0
    #         self.fvalue = cost_function(variables)
    #         print ('Self fvalue is ' + str(self.fvalue))
    #         if hasattr(D[dlabel], 'dict'):  # Tuning files
    #           for proxy, tag in D[dlabel].dict.items():
    #               D[dlabel].tuning_uncertainty.update({proxy: np.mean(D[dlabel].partsix)})
    #               print proxy, D[dlabel].tuning_uncertainty[proxy]
    #
    #     def __call__(self, x, f, accept):
    #
    #             if f < self.fvalue:
    #                 self.fvalue = f
    #                 print ('Self fvalue is ' + str(self.fvalue))
    #                 for i, dlabel in enumerate(list_drillings):
    #                     if hasattr(D[dlabel], 'dict'):  # Tuning files
    #                         for proxy, tag in D[dlabel].dict.items():
    #                             D[dlabel].tuning_uncertainty.update({proxy: np.mean(D[dlabel].partsix)})
    #                             print proxy, D[dlabel].tuning_uncertainty[proxy]
    #             self.callbackiter += 1

    def method_leastsq(*args,**kwargs):
        opt_result = scipy.optimize.OptimizeResult()
        variables, hess, infodict, mesg, ier = leastsq(residuals,args[1], Dfun=Dres, col_deriv=1, full_output=1)
        opt_result.x = variables
        opt_result.success = True
        opt_result.status = ier
        opt_result.message = mesg
        opt_result.jac = Dres(variables)
        opt_result.hess = hess
        opt_result.nfev = infodict['nfev']
        opt_result.nit = infodict['nfev']
        opt_result.maxcv = 0
        opt_result.fun = cost_function(variables)
        leak_f.append(opt_result.fun)
        if min(leak_f) == leak_f[-1]:
            leak_hess[0] = opt_result.hess
        if np.isnan(opt_result.fun):
            opt_result.fun = 10.**(100.)
        return opt_result

    steptaker = StepTaker()
    #callbackfun = CallbackFun()
    res=basinhopping(cost_function,variables, niter=basin_iter, T=temperature, take_step=steptaker, disp=True, minimizer_kwargs={'method': method_leastsq, 'hess':True, 'jac':True})
    variables = res.x
    hess = leak_hess[0]
    print 'number of iterations: ', res.nit
    print 'Message: ', res.message
elif opt_method=='none':
    print 'No optimization'
    hess=np.zeros((np.size(variables),np.size(variables)))
else:
    print opt_method,': Optimization method not recognized.'
    quit()
print 'Optimization execution time: ', time.time() - start_time_opt, 'seconds'
#print 'solution: ',variables
print 'cost function: ',cost_function(variables)
if hess==None:
    print 'singular matrix encountered (flat curvature in some direction)'
    quit()
print 'Calculation of confidence intervals'
index=0
for dlabel in list_drillings:
    D[dlabel].variables=variables[index:index+np.size(D[dlabel].variables)]
    D[dlabel].hess=hess[index:index+np.size(D[dlabel].variables),index:index+np.size(D[dlabel].variables)]
    index=index+np.size(D[dlabel].variables)
    D[dlabel].sigma()

###Final display and output
print 'Display of results'
for i,dlabel in enumerate(list_drillings):
#    print dlabel+'\n'
    D[dlabel].save()
    D[dlabel].figures()
    for j,dlabel2 in enumerate(list_drillings):
        if j<i:
#            print dlabel2+'-'+dlabel+'\n'
            DC[dlabel2+'-'+dlabel].figures()
            
###Program execution time
message='Program execution time: '+str(time.clock()-start_time)+' seconds.' 
print  message
output_file.write(message)

if show_figures:
    mpl.show()
