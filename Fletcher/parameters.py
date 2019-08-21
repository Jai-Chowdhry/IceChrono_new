list_drillings=['FP']
opt_method='MC'  #leastsq, leastsq-parallel, basinhopping, MC, none
mpi_activate=False #Use schwimmbad mpi implementation for parallel computing? Boolean
nb_nodes=1         #Number of nodes for the leastsq-parallel and basinhopping modes
temperature=20      #Temperature for basinhopping algorithm
MC_iter=100   #Number of iterations for basinhopping or MC algorithms
MC_thin=1 #Number of iterations between true proposals to guarantee independence (MC algorithm only, must be less than MC_iter)
MC_walkers = 128
MC_step = 0.2
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

