list_drillings=['EDC','TALDICE']
opt_method='MC'  #leastsq, leastsq-parallel, MC, PT, basinhopping none
mpi_activate=True #Use schwimmbad mpi implementation for parallel computing? Boolean
nb_nodes=10         #Number of nodes for the leastsq-parallel, MC and basinhopping modes
temperature=500      #Temperature for basinhopping algorithm
ntemps = 4 #Number of temperatures for parallel tempering
MC_iter=200000 #Number of iterations for basinhopping, PT or MC algorithms
MC_walkers = 128
MC_step = 0.2
MC_thin = 0.01
MC_write=False #Fixme we can get rid of this parameter because backend automatically writes...
MC_adaptive=False
MC_kombine=False
MC_restart_backend=False

#Defines the colors for the graphs
color_obs='r'       #color for the observations
color_opt='k'       #color for the posterior scenario
color_mod='b'       #color for the prior scenario
color_ci='0.8'      #color for the confidence intervals
color_sigma='m'     #color for the uncertainty
color_di='g'        #color for the dated intervals
show_initial=False  #always put to False for now
color_init='c'      #always put to 'c' for now
scale_ageci=10.     #scaling of the confidence interval in the ice and air age figures
show_figures=False  #whether to show or not the figures at the end of the run
show_airlayerthick=False #whether to show the air layer thickness figure (buggy on anaconda)
