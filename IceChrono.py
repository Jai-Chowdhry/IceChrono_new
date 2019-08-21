#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
print_synchro=False
sys.path.extend(['','/home/jcb232/nix-tests/python/venv/lib/python2.7/site-packages', '/nix/store/9w2zsb4skzs1j516wl12akywn5f086b6-python2.7-setuptools-19.4/lib/python2.7/site-packages/setuptools-19.4-py2.7.egg', '/nix/store/9w2zsb4skzs1j516wl12akywn5f086b6-python2.7-setuptools-19.4/lib/python2.7/site-packages/setuptools-19.4-py2.7.egg', '/nix/store/rv593nah61ym1yklsgm0hawxn07f6r85-python-2.7.11/lib/python2.7/site-packages', '/nix/store/40k7ssmlsih2bhzdja6xadgwv1jcyvv8-python2.7-numpy-1.10.4/lib/python2.7/site-packages', '/nix/store/d945ibfx9x185xf04b890y4f9g3cbb63-python-2.7.11/lib/python2.7/site-packages', '/nix/store/9w2zsb4skzs1j516wl12akywn5f086b6-python2.7-setuptools-19.4/lib/python2.7/site-packages', '/nix/store/5f2050cfj9r09sadkadq1rxpcqylrvmp-python2.7-scipy-0.17.0/lib/python2.7/site-packages', '/nix/store/s1zhj6w7lgk5y5i6131dcvmpdji0vd56-python2.7-pip-8.0.2/lib/python2.7/site-packages', '/nix/store/ra025z4fajhnnsz631d9k4kx64wn552n-python2.7-matplotlib-1.5.0/lib/python2.7/site-packages', '/nix/store/w7v6i3ckznl4kffdiyf0a3fza3pjik4s-python2.7-cycler-0.9.0/lib/python2.7/site-packages', '/nix/store/n086vp0yfjb9lnlzcsf479s2l934zvyf-python2.7-six-1.10.0/lib/python2.7/site-packages', '/nix/store/imrn17l9hmz29z9c9z66p09swgighb6q-python2.7-dateutil-2.4.2/lib/python2.7/site-packages', '/nix/store/bmkfsfwch6jj7il5nz53wxy7wdvbsj33-python2.7-nose-1.3.7/lib/python2.7/site-packages', '/nix/store/69z9cq3wp14cm1fig61x1k3vxs2nfhky-python2.7-coverage-4.0.1/lib/python2.7/site-packages', '/nix/store/s2wrbk2vndy44wdkn35ffc8zg3q63ykm-python2.7-pyparsing-2.0.1/lib/python2.7/site-packages', '/nix/store/flw17nylvz02iif8i5m5a1p5194qhc3c-python2.7-tornado-4.2.1/lib/python2.7/site-packages', '/nix/store/b0pf625rinj4b9z34rlmzcb6b1s15j8i-python2.7-backports.ssl_match_hostname-3.4.0.2/lib/python2.7/site-packages', '/nix/store/qmjg50bz6i7wh0nfxkj1p78164253ii6-python2.7-certifi-2015.9.6.2/lib/python2.7/site-packages', '/nix/store/4zdaqg59n7c9yg4ibnw3xq3gxv7mfahl-python2.7-pkgconfig-1.1.0/lib/python2.7/site-packages', '/nix/store/mw82kzsx06ry9a45vnn92xs86d8kq4r9-python2.7-mock-1.3.0/lib/python2.7/site-packages', '/nix/store/g9z0rnr68nb0asznw1x0y2g4jplnn5bi-python2.7-funcsigs-0.4/lib/python2.7/site-packages', '/nix/store/2whif46cpv8jv178bl23vyq8d0vm4rr0-python2.7-pbr-1.8.1/lib/python2.7/site-packages', '/nix/store/szwk9q45am71mkwlllc7r447rpdg5hs1-python2.7-pytz-2015.7/lib/python2.7/site-packages', '/nix/store/vzsfdbc3dnzrdqywz0mvvi06f7f1i92q-python2.7-virtualenv-13.1.2/lib/python2.7/site-packages', '/nix/store/9d8lbpa8lvla8s4q40pijg1ny7gya3qm-python-readline-2.7.11/lib/python2.7/site-packages', '/nix/store/a01d3bvj4nj8pkjdr6jq7wvpqzyi4qz6-python-sqlite3-2.7.11/lib/python2.7/site-packages', '/nix/store/pbsbjh1ma1gp6icbzv9y4vhls1kkgbif-python-curses-2.7.11/lib/python2.7/site-packages', '/nix/store/rv593nah61ym1yklsgm0hawxn07f6r85-python-2.7.11/lib/python27.zip', '/nix/store/rv593nah61ym1yklsgm0hawxn07f6r85-python-2.7.11/lib/python2.7', '/nix/store/rv593nah61ym1yklsgm0hawxn07f6r85-python-2.7.11/lib/python2.7/plat-linux2', '/nix/store/rv593nah61ym1yklsgm0hawxn07f6r85-python-2.7.11/lib/python2.7/lib-tk', '/nix/store/rv593nah61ym1yklsgm0hawxn07f6r85-python-2.7.11/lib/python2.7/lib-old', '/nix/store/rv593nah61ym1yklsgm0hawxn07f6r85-python-2.7.11/lib/python2.7/lib-dynload', '/home/jcb232/.local/lib/python2.7/site-packages'])


#TODO: what about symbolic links in github?

import sys
import time
import math as m
import numpy as np
import matplotlib.pyplot as mpl
import multiprocessing
import tqdm
import warnings
import os
import scipy.linalg #
from scipy.optimize import leastsq, minimize, basinhopping
import emcee #
from schwimmbad import MPIPool, MultiPool
import scipy.sparse.linalg
import gc
import more_itertools as mit
import copy
import shutil
from sksparse.cholmod import cholesky as cholesky_sparse #
from matplotlib.backends.backend_pdf import PdfPages
from scipy.linalg import cholesky

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
execfile('DEMove_modified.py') #Fixme maybe we can import this class

##Global
variables=np.array([])
sigmas=np.array([])
D={}
DC={}

##Functions and Classes


execfile('IceChronoClasses.py')

def residuals(var): #FIXME something gets unintentionally updated here with the wrong timing...
    """Calculate the residuals."""
    resi=np.array([])
    index=0
    for i,dlabel in enumerate(list_drillings):
        D[dlabel].variables=var[index:index+np.size(D[dlabel].variables)]
        index=index+np.size(D[dlabel].variables)
    for i, dlabel in enumerate(list_drillings):
        resi=np.concatenate((resi,D[dlabel].residuals(D[dlabel].variables)))
        for j,dlabel2 in enumerate(list_drillings):
            if j<i:
                resi=np.concatenate((resi,DC[dlabel2+'-'+dlabel].residuals()))
    return resi

def residuals_obs(var):
    """Calculate the residuals."""
    resi=np.array([])
    index=0
    for i,dlabel in enumerate(list_drillings):
        D[dlabel].variables=var[index:index+np.size(D[dlabel].variables)]
        index=index+np.size(D[dlabel].variables)
    for i, dlabel in enumerate(list_drillings):
        resi=np.concatenate((resi,D[dlabel].residuals_separated(D[dlabel].variables)[1]))
        for j,dlabel2 in enumerate(list_drillings):
            if j<i:
                resi=np.concatenate((resi,DC[dlabel2+'-'+dlabel].residuals()))
    return resi

def residuals_prior(var):
    """Calculate the residuals."""
    resi=np.array([])
    index=0
    for i,dlabel in enumerate(list_drillings):
        D[dlabel].variables=var[index:index+np.size(D[dlabel].variables)]
        index=index+np.size(D[dlabel].variables)
    for i, dlabel in enumerate(list_drillings):
        resi=np.concatenate((resi,D[dlabel].residuals_separated(D[dlabel].variables)[0]))
        #for j,dlabel2 in enumerate(list_drillings):
        #    if j<i:
        #        resi=np.concatenate((resi,DC[dlabel2+'-'+dlabel].residuals()))
    return resi

blob2 = {}
def cost_function(var):
    #gc.collect()
    cost=np.dot(residuals(var),np.transpose(residuals(var)))
    #cost=np.dot(residuals(var),np.transpose(residuals(var)))
    #gc.collect()
    return cost #, blob2

#blob = {}
def cost_function_negative(var):
    #depthlist = []
    #for dlabel in list_drillings:
    #    bloblist.append([D[dlabel].depth, D[dlabel].age, D[dlabel].airage])
    cost=-np.dot(residuals(var),np.transpose(residuals(var)))
    #cost=-np.dot(residuals(var),np.transpose(residuals(var)))
    #blob = np.array(bloblist)
    return cost

def cost_function_observations(var): #Fixme how to add blobs to parallel temp?
    #for dlabel in list_drillings:
    #    blob.update({dlabel : [D[dlabel].depth, D[dlabel].age, D[dlabel].airage]})
    cost=-np.dot(residuals_obs(var),np.transpose(residuals_obs(var)))
    return cost

def cost_function_prior(var): #Fixme how to add blobs to parallel temp?
    #for dlabel in list_drillings:
    #    blob.update({dlabel : [D[dlabel].depth, D[dlabel].age, D[dlabel].airage]})
    cost=-np.dot(residuals_prior(var),np.transpose(residuals_prior(var)))
    return cost

def prior_dummy(var):
    return 1.0

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
        
    D[dlabel]=Record(dlabel)
    D[dlabel].init()
    D[dlabel].model(D[dlabel].variables)
#    D[dlabel].a_init=D[dlabel].a
#    D[dlabel].LID_init=D[dlabel].LID
    D[dlabel].write_init()
#    D[dlabel].display_init()
    variables=np.concatenate((variables,D[dlabel].variables))
    D[dlabel].synchro_covar()

for i,dlabel in enumerate(list_drillings):
    for j,dlabel2 in enumerate(list_drillings):
        if j<i:
            print 'Initialization of drilling pair '+dlabel2+'-'+dlabel
            DC[dlabel2+'-'+dlabel]=RecordPair(D[dlabel2],D[dlabel])
            DC[dlabel2+'-'+dlabel].init()
#            DC[dlabel2+'-'+dlabel].display_init()

def complist(counter): #Hack for parallelization
    return variables

def corrmat(A): #Hack, thanks to https://github.com/CamDavidsonPilon/Python-Numerics/blob/master/utils/cov2corr.py
    d = np.sqrt(A.diagonal())
    A = ((A.T/d).T)/d
    return A

def emcee_init_all(variables): #Hack for parallelization
    fit = []
    for dlabel in list_drillings:
        if hasattr(D[dlabel],'reloadcovariance') and D[dlabel].reloadcovariance:
            filename=datadir+dlabel+'/covariance.npz'
            cmatrix = np.load(filename)
            c_tau=corrmat(cmatrix['tau'])
            c_a=corrmat(cmatrix['a'])
            c_LID=corrmat(cmatrix['LID'])
        else:
            c_tau = D[dlabel].correlation_corr_tau
            c_a = D[dlabel].correlation_corr_a
            c_LID = D[dlabel].correlation_corr_LID
        #print(np.shape(D[dlabel].corr_tau) , np.shape(c_tau))
        random_corr_tau = D[dlabel].corr_tau + np.random.multivariate_normal(D[dlabel].corr_tau, c_tau) *D[dlabel].MC_init_width #+ np.random.uniform(-D[dlabel].MC_init_width,D[dlabel].MC_init_width,np.shape(D[dlabel].corr_tau))
        random_corr_a = D[dlabel].corr_a + np.random.multivariate_normal(D[dlabel].corr_a, c_a)  *D[dlabel].MC_init_width #+ np.random.uniform(-D[dlabel].MC_init_width,D[dlabel].MC_init_width,np.shape(D[dlabel].corr_a))
        random_corr_LID = D[dlabel].corr_LID + np.random.multivariate_normal(D[dlabel].corr_LID, c_LID)  *D[dlabel].MC_init_width #+ np.random.uniform(-D[dlabel].MC_init_width,D[dlabel].MC_init_width,np.shape(D[dlabel].corr_LID))
        if D[dlabel].corr_tau_depth[0] == 0.0: #Fixme this initializes all of the surface values of thinning to 0. Implement everywhere?
            random_corr_tau[0] = 0.0 #Fixme doesn't work because of random component of proposals.
        fit = np.concatenate((fit,random_corr_tau,random_corr_a,random_corr_LID))
    return fit

##Optimization
start_time_opt = time.time()
print 'cost function: ',cost_function(variables)
print 'Variable size ', np.size(variables)
#print 'cost function negative: ',cost_function_negative(variables)
#print 'cost function observation: ',cost_function_observations(variables)
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

elif opt_method == 'MC':

    mc_backend = emcee.backends.HDFBackend(datadir+'saved_iterations.hdf5',name='mcmc',read_only=False) #See emcee backends documentation for how to load

    print 'Optimization by [parallel] affine invariant MCMC'

    pos0 = []

    #print'Initializing walkers'

    if not mpi_activate:
        pool=MultiPool(processes=nb_nodes)
        #if not pool.is_master():
        #    pool.wait()
        #    sys.exit(0)
        variables_list = map(complist, range(
            MC_walkers))  # FIXME reset to pool.map after tests
        #print("Variable list size ", np.shape(variables_list))
        pos0 = map(emcee_init_all, variables_list)
        pos0[0] = variables #Hack to ensure the non-perturbed initial position is included #FIXME why is this not working???
        steps = emcee.EnsembleSampler(MC_walkers, np.size(variables), cost_function_negative, a=MC_step,
                                      backend=mc_backend, moves=emcee.moves.StretchMove(a=MC_step, live_dangerously=True)) #Fixme add pool=pool after profiling

        print 'Walkers initialized. Starting optimization'

        steps.run_mcmc(pos0, nsteps=MC_iter, thin=MC_thin, progress=True)

        pool.close()

    if mpi_activate: #We wait until here to send a copy of the initialised variables, classes and functions to each node.
        with MPIPool() as pool:
            if not pool.is_master():
                pool.wait()
                sys.exit(0)
            if MC_restart_backend:
                print('Restarting from backend')
                variables_list = mc_backend.get_chain()[0]
                pos0 = variables_list
                print np.shape(variables_list)
            else:
                mc_backend.reset(MC_walkers,len(variables))
                variables_list = pool.map(complist,range(MC_walkers)) #FIXME why not include variables as an argument? Would allow for earlier parallelization
                #print("Variable list size ", np.shape(variables_list))
                pos0 = pool.map(emcee_init_all, variables_list)
            pos0[0] = variables  # Hack to ensure the non-perturbed initial position is included. #FIXME why is this not working???
            #steps = emcee.EnsembleSampler(MC_walkers, np.size(variables), cost_function_negative, pool=pool, moves=emcee.moves.StretchMove(a=MC_step,live_dangerously=True))

            moves = [
                (emcee.moves.DESnookerMove(gammas=(2.38 / np.sqrt(2 * np.size(variables))),live_dangerously=True), 0.4),
                (emcee.moves.DEMove(gamma0=(2.38 / np.sqrt(2 * np.size(variables))), live_dangerously=True), 0.5),
                (emcee.moves.DEMove(gamma0=0.98,live_dangerously=True), 0.1)]

            if MC_adaptive:
                moves_adaptive = [
                    (TerBraak2008_snookermove(backend=mc_backend,live_dangerously=True), 0.1),
                    (TerBraak2008_DEMove(backend=mc_backend,live_dangerously=True), 0.8),
                    (TerBraak2008_DEMove(gamma0=0.98,backend=mc_backend,live_dangerously=True), 0.1),]

            if not MC_adaptive and not MC_kombine:
                steps = emcee.EnsembleSampler(MC_walkers, np.size(variables), cost_function_negative, pool=pool, backend=mc_backend, moves=moves)
            #steps = emcee.EnsembleSampler(MC_walkers, np.size(variables), cost_function_negative, pool=pool, moves=emcee.moves.DESnookerMove(live_dangerously=True))
            elif MC_adaptive:
                steps = emcee.EnsembleSampler(MC_walkers, np.size(variables), cost_function_negative, pool=pool,
                                          backend=mc_backend, moves=moves_adaptive)
            elif MC_kombine:
                import kombine
                steps = kombine.Sampler(MC_walkers, np.size(variables), cost_function_negative, pool=pool)
                steps.burnin(p0=pos0,max_steps=100,verbose=True) #testing emcee for burnin
                for result in tqdm.tqdm (steps.sample(p0=pos0,iterations=MC_iter),total=MC_iter):
                    pass
                #steps.run_mcmc(p0=pos0, N=MC_iter)

            print 'Walkers initialized. Starting optimization'
            if not MC_kombine:
                steps.run_mcmc(pos0, nsteps=MC_iter,  progress=True, tune=True, thin=MC_thin) #Tune should only work with moves_adaptive
            #Fixme, with adaptive should we then stop, clear and rerun without adaptation???
    if not MC_kombine:
        maxindex = np.argmax(steps.get_log_prob(flat=True))
        #print maxindex
        #print np.shape(steps.get_log_prob(flat=True))
        a=mpl.figure()
        mpl.plot(steps.get_log_prob(flat=True))
        mpl.plot(maxindex,steps.get_log_prob(flat=True)[maxindex],'r*')
        print(steps.get_log_prob(flat=True)[maxindex])
        a.savefig('cost_evolution.png')
        maxfit = steps.get_chain(flat=True)[maxindex]
        #print np.shape(steps.get_chain(flat=True))
        print steps.acceptance_fraction
        variables = res = maxfit
        hess = np.zeros((np.size(variables), np.size(variables)))
        np.savetxt('bestresi.txt', residuals(variables))
        #b = mpl.figure()
        autocorr_time = steps.get_autocorr_time(c=5, tol=50, quiet=True) #Will send a warning if autocorrelation time is too large/Not enough iterations
        #b.savefig('autocorrelationtime_chain0.png')
    elif MC_kombine:
        maxindex = np.argmax(steps.lnpost.flatten())
        a = mpl.figure()
        mpl.plot(steps.lnpost.flatten())
        mpl.plot(maxindex, steps.lnpost.flatten()[maxindex], 'r*')
        print(steps.lnpost.flatten()[maxindex])
        a.savefig('cost_evolution.png')
        maxindex = np.argmax(steps.lnpost)
        maxfit = steps.get_samples()[maxindex]
        print steps.acceptance_fraction
        variables = res = maxfit
        hess = np.zeros((np.size(variables), np.size(variables)))
        np.savetxt('bestresi.txt', residuals(variables))
        b = mpl.figure()
        mpl.hist(steps.autocorrelation_times)
        b.savefig('autocorrelationtime_chain0.png')

    #maxindex = np.unravel_index(np.argmax(steps.lnprobability), np.shape(steps.lnprobability))
    #maxfit = steps.chain[maxindex]
    #variables = res = maxfit
    hess = np.zeros((np.size(variables), np.size(variables)))

elif opt_method == 'PT':
    import ptemcee
    print 'Optimization by parallel tempering ensemble MCMC'

    pos0 = []

    print'Initializing walkers'

    if not mpi_activate:
        variables_list = pool.map(complist, range(
            MC_walkers))  # FIXME why not include variables as an argument? Would allow for earlier parallelization
        print("Variable list size ", np.shape(variables_list))
        pos0 = pool.map(emcee_init_all, variables_list)
        pos0 = np.array(pos0).reshape((ntemps, MC_walkers, len(variables)))

    if mpi_activate: #We wait until here to send a copy of the initialised variables, classes and functions to each node.

            pool = MPIPool()
            if not pool.is_master():
                pool.wait()
                sys.exit(0)

            variables_list = pool.map(complist,range(MC_walkers*ntemps))
            print("Variable list size ", np.shape(variables_list))
            pos0 = pool.map(emcee_init_all, variables_list)
            pos0 = np.array(pos0).reshape((ntemps, MC_walkers, len(variables)))
            steps = ptemcee.Sampler(MC_walkers, dim=np.size(variables), logl=cost_function_observations,
                                    logp=cost_function_prior, ntemps=ntemps, adaptation_lag=1000,adaptation_time=10, pool=pool)
            print 'Walkers initialized. Starting optimization'
            for result in tqdm.tqdm(steps.sample(pos0, iterations=MC_iter),total=MC_iter):  #,thin=MC_thin)):
                pass
                #if (i) % 1000 == 0:
                #if (i) % 1 == 0:
                #    print str(i) + " iterations completed"
                #    #steps.clear_blobs()
            pool.close()

    else:
        print np.size(variables)
        steps = emcee.PTSampler(ntemps, MC_walkers, np.size(variables), cost_function_observations, cost_function_prior, threads=10)

        print 'Walkers initialized. Starting optimization'

        i = 0
        print np.shape(pos0)
        for prob, lnprob, lnlike in steps.sample(pos0, iterations=MC_iter ,thin=MC_thin):
            #if (i) % 1000 == 0:
            i += 1
            if (i) % 10 == 0:
                print str(i) + " iterations completed"

    #maxindex = np.argmax(steps.lnlikelihood)
    #print maxindex
    #print np.shape(steps.get_log_prob(flat=True))
    #a = mpl.figure()
    #mpl.plot(steps.lnlikelihood)
    #mpl.plot(maxindex, steps.get_log_prob(flat=True)[maxindex], 'r*')
    #print(steps.get_log_prob(flat=True)[maxindex])
    #a.savefig('cost_evolution.png')
    b = mpl.figure()
    #print steps.get_autocorr_time()
    mpl.hist(steps.get_autocorr_time())
    b.savefig('autocorrelationtime_chain0.png')
    #maxfit = steps.get_chain(flat=True)[maxindex]
    #print np.shape(steps.get_chain(flat=True))
    print steps.acceptance_fraction
    #variables = res = maxfit
    hess = np.zeros((np.size(variables), np.size(variables)))
    #np.savetxt('bestresi.txt', residuals(variables))
elif opt_method=='basinhopping':
    print 'Optimization by basin hopping'
    leak_f = []
    leak_hess = [np.zeros((np.size(variables),np.size(variables)))]

    class StepTaker(object):
        def __init__(self,  stepsize = 0.05): #0.01
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
        #if np.isnan(opt_result.fun):
        #    opt_result.fun = 10.**(100.)
        return opt_result

    steptaker = StepTaker()
    #callbackfun = CallbackFun()
    res=basinhopping(cost_function,variables, niter=MC_iter, T=temperature, take_step=steptaker, disp=True, minimizer_kwargs={'method': method_leastsq, 'hess':True, 'jac':True})
    variables = res.x
    hess = leak_hess[0]
    print 'number of iterations: ', res.nit
    print 'Message: ', res.message

elif opt_method=='none':
    print 'No optimization'
#    hess=np.zeros((np.size(variables),np.size(variables)))
else:
    print opt_method,': Optimization method not recognized.'
    quit()

print 'Optimization execution time: ', time.time() - start_time_opt, 'seconds'
#print 'solution: ',variables
time.sleep(0.001)
#print 'cost function: ',cost_function(variables)
#np.savetxt('residuals.txt',residuals(variables))

if opt_method!='none' and hess.__str__()=='None':
    print 'singular matrix encountered (flat curvature in some direction)'
    quit()
print 'Calculation of confidence intervals'
index=0
for dlabel in list_drillings:
    if opt_method=='none': #Fixme
        D[dlabel].sigma_zero()
    else:
        D[dlabel].variables=variables[index:index+np.size(D[dlabel].variables)]
        D[dlabel].hess=hess[index:index+np.size(D[dlabel].variables),index:index+np.size(D[dlabel].variables)]
        index=index+np.size(D[dlabel].variables)
        D[dlabel].sigma()

print 'cost function: ',cost_function(variables)
print_synchro=False
# 'cost function negative: ',cost_function_negative(variables)
#print 'cost function prior', cost_function_prior(variables)
#print 'cost function observation: ',cost_function_observations(variables)
#print 'cost function: ',cost_function(variables)

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
message=                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        'Program execution time: '+str(time.clock()-start_time)+' seconds.'
print  message
output_file.write(message)

if show_figures:
    mpl.show()
