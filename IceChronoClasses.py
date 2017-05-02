#TODO: extend the chronology down to the bedrock by extrapolating the accumulation
#TODO: optinally use a restart file to have a bootstrap method
#TODO: is there an elegant way to unpack the variables vector in the model function?
#TODO: allow to save the correction vector to be able to restart while changing the resolution
#TODO: include some checks for when dDdepth/dz>1
#TODO: Delta-depth observations should be lognormal?
#TODO: we should superpose two charts for ice and air ages, one for the age and one for the uncertainty, since the min age is not always near 0.
#TODO: also compute the prior uncertainties and show them in the figures.
#TODO: the reading of observations does not work if there is only one observation (since the read matrix is 1D in this case).
#TODO: is there really a computation gain with the change of variable for the correction functions? Avoiding this change of variables would make the code easier to understand. I think there is no gain since solving A^-1 b when we have the LU factorisation of A does not cost more than computing A^-1 * b when we have computed A^-1.

#TODO: Comment all lines for implementation of orbital tuning.

def interp1d_extrap(x,y):
    def f(xp):
        g=interp1d(x,y, bounds_error=False)
        return np.where(xp<x[0],y[0],np.where(xp>x[-1],y[-1],g(xp)))    
    return f

def interp1d_extrap_linear(x,y):
    def f(xp):
        g=interp1d(x,y, bounds_error=False)
        return np.where(xp<x[0],(xp-x[0])*(y[1]-y[0])/(x[1]-x[0])+x[0],np.where(xp>x[-1],(xp-x[-1])*(y[-1]-y[-2])/(x[-1]-x[-2])+x[-1],g(xp)))
    return f

def interp1d_lin_aver_extrap(x, y):
    def f(xp):
        yp=np.nan*np.zeros(np.size(xp)-1)
        if xp[0]<min(x):
            xmod=np.concatenate((np.array([xp[0]]),x))
            ymod=np.concatenate((np.array([y[0]]),y))
        else:
            xmod=x+0
            ymod=y+0
        if xp[-1]>max(x):
            xmod=np.concatenate((xmod,np.array([xp[-1]])))
            ymod=np.concatenate((ymod,np.array([y[-1]])))
        for i in range(np.size(xp)-1):
            xx=xmod[np.where(np.logical_and(xmod>xp[i],xmod<xp[i+1]))]
            xx=np.concatenate((np.array([xp[i]]),xx,np.array([xp[i+1]])))
            f=interp1d(xmod,ymod)
            yy=f(xx)
            yp[i]=np.sum((yy[1:]+yy[:-1])/2*(xx[1:]-xx[:-1]))/(xp[i+1]-xp[i])
        return yp

    return f

def interp1d_stair_aver_extrap(x, y):
    def f(xp):
        xmod=x+0
        ymod=y+0
        if xp[0]<x[0]:
            xmod=np.concatenate((np.array([xp[0]]),xmod))
            ymod=np.concatenate((np.array([y[0]]),ymod))
        if xp[-1]>x[-1]:
            xmod=np.concatenate((xmod,np.array([xp[-1]])))
            ymod=np.concatenate((ymod,np.array([y[-1]])))
        yint=np.cumsum(np.concatenate((np.array([0]),ymod[:-1]*(xmod[1:]-xmod[:-1]))))
        g=interp1d(xmod,yint)
        yp=(g(xp[1:])-g(xp[:-1]))/(xp[1:]-xp[:-1])     #Maybe this is suboptimal since we compute twice g(xp[i])
        return yp

    return f

def gaussian(x):
    return np.exp(-x**2/2)

class Drilling:

    def __init__(self, dlabel):
        self.label=dlabel

    def init(self):

#        print 'Initialization of drilling '+self.label

        self.accu_prior_rep='staircase'

        execfile(datadir+'/parameters-AllDrillings.py')
        execfile(datadir+self.label+'/parameters.py')

        self.depth_mid=(self.depth[1:]+self.depth[:-1])/2
        self.depth_inter=(self.depth[1:]-self.depth[:-1])

## We set up the raw model

        if self.calc_a:
            readarray=np.loadtxt(datadir+self.label+'/isotopes.txt')
            if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
            self.iso_depth=readarray[:,0]
            if self.calc_a_method=='fullcorr':
                self.iso_d18Oice=readarray[:,1]
                f=interp1d_stair_aver_extrap(self.iso_depth, self.iso_d18Oice)
                self.d18Oice=f(self.depth)
                self.iso_deutice=readarray[:,2]
                f=interp1d_stair_aver_extrap(self.iso_depth, self.iso_deutice)
                self.deutice=f(self.depth)
                self.iso_d18Osw=readarray[:,3]
                f=interp1d_stair_aver_extrap(self.iso_depth, self.iso_d18Osw)
                self.d18Osw=f(self.depth)
                self.excess=self.deutice-8*self.d18Oice   # dans Uemura : d=excess
                self.a=np.empty_like(self.deutice)
                self.d18Oice_corr=self.d18Oice-self.d18Osw*(1+self.d18Oice/1000)/(1+self.d18Osw/1000)	#Uemura (1)
                self.deutice_corr=self.deutice-8*self.d18Osw*(1+self.deutice/1000)/(1+8*self.d18Osw/1000)	#Uemura et al. (CP, 2012) (2) 
                self.excess_corr=self.deutice_corr-8*self.d18Oice_corr
                self.deutice_fullcorr=self.deutice_corr+self.gamma_source/self.beta_source*self.excess_corr
            elif self.calc_a_method=='deut':
                self.iso_deutice=readarray[:,1]
                f=interp1d_stair_aver_extrap(self.iso_depth, self.iso_deutice)
                self.deutice_fullcorr=f(self.depth)
            elif selc.calc_a_method=='d18O':
                self.d18Oice=readarray[:,1]
                f=interp1d_stair_aver_extrap(self.iso_depth, self.iso_d18Oice)
                self.deutice_fullcorr=8*f(depth)
            else:
                print 'Accumulation method not recognized'
                quit()
        else:
            readarray=np.loadtxt(datadir+self.label+'/accu-prior.txt')
            if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
            self.a_depth=readarray[:,0]
            self.a_a=readarray[:,1]
            if readarray.shape[1]>=3:
                self.a_sigma=readarray[:,2]
            if self.accu_prior_rep=='staircase':
                f=interp1d_stair_aver_extrap(self.a_depth, self.a_a)
            elif self.accu_prior_rep=='linear':
                f=interp1d_lin_aver_extrap(self.a_depth,self.a_a)    #FIXME: We should implement a interp1d_lin_aver_extrap function
            else:
                print 'Representation of prior accu scenario not recognized'
            self.a_model=f(self.depth)
            self.a=self.a_model


        
        self.age=np.empty_like(self.depth)
        self.airage=np.empty_like(self.depth)
        

        readarray=np.loadtxt(datadir+self.label+'/density-prior.txt')
#        self.density_depth=readarray[:,0]
        if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
        self.D_depth=readarray[:,0]
        self.D_D=readarray[:,1]
        f=interp1d_extrap(self.D_depth, self.D_D)
        self.D=f(self.depth_mid)

        self.iedepth=np.cumsum(np.concatenate((np.array([0]), self.D*self.depth_inter)))
        self.iedepth_mid=(self.iedepth[1:]+self.iedepth[:-1])/2
        if self.calc_tau:
            self.thickness_ie=self.thickness-self.depth[-1]+self.iedepth[-1]
        
        if self.calc_LID:
            if self.depth[0]<self.LID_value:
                self.LID_depth=np.array([self.depth[0], self.LID_value, self.depth[-1]])
                self.LID_LID=np.array([self.depth[0], self.LID_value, self.LID_value])
            else:
                self.LID_depth=np.array([self.depth[0], self.depth[-1]])
                self.LID_LID=np.array([self.LID_value, self.LID_value])
        else:
#            self.LID_model=np.loadtxt(datadir+self.label+'/LID-prior.txt')
            readarray=np.loadtxt(datadir+self.label+'/LID-prior.txt')
            if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
            self.LID_depth=readarray[:,0]
            self.LID_LID=readarray[:,1]
            if readarray.shape[1]>=3:
                self.LID_sigma=readarray[:,2]
        f=interp1d_extrap(self.LID_depth, self.LID_LID)
        self.LID_model=f(self.depth)




        self.Ddepth=np.empty_like(self.depth)
        self.udepth=np.empty_like(self.depth)

#        print 'depth_mid ', np.size(self.depth_mid)
#        print 'zeta ', np.size(self.zeta)
        if self.calc_tau:
            self.thicknessie=self.thickness-self.depth[-1]+self.iedepth[-1]
            self.zeta=(self.thicknessie-self.iedepth_mid)/self.thicknessie  #FIXME: maybe we should use iedepth and thickness_ie here?
            self.tau=np.empty_like(self.depth_mid)
        else:
            readarray=np.loadtxt(datadir+self.label+'/thinning-prior.txt')
            if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
            self.tau_depth=readarray[:,0]
            self.tau_tau=readarray[:,1]
            if readarray.shape[1]>=3:
                self.tau_sigma=readarray[:,2]
            f=interp1d_extrap(self.tau_depth, self.tau_tau)
            self.tau_model=f(self.depth_mid)
            self.tau=self.tau_model

        self.raw_model()

## Now we set up the correction functions

        if self.start=='restart':
            vartemp = np.loadtxt(datadir+self.label+'/restart.txt')
            xtemp = np.loadtxt(datadir+self.label+'/restart_x.txt')
            corr_tau_length, corr_a_length, corr_LID_length = vartemp[0], vartemp[1], vartemp[2]
            index = 3
            corr_tau_temp = vartemp[index:index + corr_tau_length]
            corr_a_temp = vartemp[index + corr_tau_length: index + corr_tau_length + corr_a_length]
            corr_LID_temp = vartemp[index + corr_tau_length + corr_a_length: index + corr_tau_length + corr_a_length + corr_LID_length]
            corr_tau_x = xtemp[index:index + corr_tau_length]
            corr_a_x = xtemp[index + corr_tau_length: index + corr_tau_length + corr_a_length]
            corr_LID_x = xtemp[index + corr_tau_length + corr_a_length: index + corr_tau_length + corr_a_length + corr_LID_length]
            interptau = interp1d_extrap(corr_tau_x, corr_tau_temp)
            interpa = interp1d_extrap(corr_a_x,corr_a_temp)
            interpLID = interp1d_extrap(corr_LID_x,corr_LID_temp)
            self.corr_a = interpa(self.corr_a_age)
            self.corr_tau = interptau(self.corr_tau_depth)
            self.corr_LID = interpLID(self.corr_LID_age)

        elif self.start=='default':
            self.corr_a=np.zeros(np.size(self.corr_a_age))
            self.corr_LID=np.zeros(np.size(self.corr_LID_age))
            self.corr_tau=np.zeros(np.size(self.corr_tau_depth))
        elif self.start=='random':
            self.corr_a=np.random.normal(loc=0., scale=1., size=np.size(self.corr_a_age))
            self.corr_LID=np.random.normal(loc=0., scale=1., size=np.size(self.corr_LID_age))
            self.corr_tau=np.random.normal(loc=0., scale=1., size=np.size(self.corr_tau_depth))
        else:
            print 'Start option not recognized. Try restart, default or random.'

## Now we set up the correlation matrices

        self.correlation_corr_a=np.diag(np.ones(np.size(self.corr_a)))
        self.correlation_corr_LID=np.diag(np.ones(np.size(self.corr_LID)))
        self.correlation_corr_tau=np.diag(np.ones(np.size(self.corr_tau)))

        self.chol_a=np.diag(np.ones(np.size(self.corr_a)))
        self.chol_LID=np.diag(np.ones(np.size(self.corr_LID)))
        self.chol_tau=np.diag(np.ones(np.size(self.corr_tau)))



## Definition of the covariance matrix of the background

        try:
            xx=np.where(self.a_depth<=self.depth[-1],self.fct_age_model(np.minimum(self.a_depth,self.depth[-1])),np.nan)
            yy=np.where(self.a_depth<=self.depth[-1],self.a_sigma,np.nan)
            f=interp1d_extrap(xx,yy) #, bounds_error=False, fill_value=self.a_sigma[-1])
            self.sigmap_corr_a=f(self.corr_a_age)           #FIXME: we should average here since it would be more representative
        except AttributeError:
            print 'Sigma on prior accu scenario not defined in the accu-prior.txt file'

        try:
            xx=np.where(self.LID_depth<=self.depth[-1],self.fct_airage_model(np.minimum(self.LID_depth,self.depth[-1])),np.nan)
            xx=np.concatenate((np.array([self.age_top]),xx))
            yy=np.where(self.LID_depth<=self.depth[-1],self.LID_sigma,np.nan)
            yy=np.concatenate((np.array([self.LID_sigma[0]]),yy))
            f=interp1d_extrap(xx,yy)  #, bounds_error=False, fill_value=self.LID_sigma[-1])
            self.sigmap_corr_LID=f(self.corr_LID_age)           #FIXME: we should average here since it would be more representative
        except AttributeError:
            print 'Sigma on prior LID scenario not defined in the LID-prior.txt file'

        try:
            f=interp1d_extrap(self.tau_depth,self.tau_sigma)  #, bounds_error=False, fill_value=self.tau_sigma[-1])
            self.sigmap_corr_tau=f(self.corr_tau_depth)           #FIXME: we should average here since it would be more representative
        except AttributeError:
            print 'Sigma on prior thinning scenario not defined in the thinning-prior.txt file'

        self.correlation_corr_a_before=self.correlation_corr_a+0
        self.correlation_corr_LID_before=self.correlation_corr_LID+0
        self.correlation_corr_tau_before=self.correlation_corr_tau+0

        filename=datadir+'/parameters-CovariancePrior-AllDrillings-init.py'
        if os.path.isfile(filename):
            execfile(filename)
        filename=datadir+self.label+'/parameters-CovariancePrior-init.py'
        if os.path.isfile(filename):
            execfile(filename)


        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_a=cholesky(self.correlation_corr_a)
        if (self.correlation_corr_LID_before!=self.correlation_corr_LID).any():
            self.chol_LID=cholesky(self.correlation_corr_LID)
        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_tau=cholesky(self.correlation_corr_tau)

        self.variables=np.array([])
        #        if self.calc_a==True:
        #            self.variables=np.concatenate((self.variables, np.array([self.A0]), np.array([self.beta])))
        #        if self.calc_tau==True:
        #            self.variables=np.concatenate((self.variables, np.array([self.pprime]), np.array([self.muprime])))
        self.variables = np.concatenate((self.variables, self.corr_tau, self.corr_a, self.corr_LID))
        self.sigmap_all=np.array([])
        self.sigmap_all=np.concatenate((self.sigmap_all, self.sigmap_corr_tau, self.sigmap_corr_a, self.sigmap_corr_LID))

#Reading of observations

        filename=datadir+self.label+'/ice_age.txt'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename))>0:
                readarray=np.loadtxt(filename)
                if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                self.icemarkers_depth=readarray[:,0]
                self.icemarkers_age=readarray[:,1]
                self.icemarkers_sigma=readarray[:,2]
            else:
                self.icemarkers_depth=np.array([])
                self.icemarkers_age=np.array([])
                self.icemarkers_sigma=np.array([])

        filename=datadir+self.label+'/air_age.txt'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename))>0:
                readarray=np.loadtxt(filename)
                if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                self.airmarkers_depth=readarray[:,0]
                self.airmarkers_age=readarray[:,1]
                self.airmarkers_sigma=readarray[:,2]
            else:
                self.airmarkers_depth=np.array([])
                self.airmarkers_age=np.array([])
                self.airmarkers_sigma=np.array([])

        filename=datadir+self.label+'/ice_age_intervals.txt'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename))>0:
                readarray=np.loadtxt(filename)
                if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                self.iceintervals_depthtop=readarray[:,0]
                self.iceintervals_depthbot=readarray[:,1]
                self.iceintervals_duration=readarray[:,2]
                self.iceintervals_sigma=readarray[:,3]
            else:
                self.iceintervals_depthtop=np.array([])
                self.iceintervals_depthbot=np.array([])
                self.iceintervals_duration=np.array([])
                self.iceintervals_sigma=np.array([])

        filename=datadir+self.label+'/air_age_intervals.txt'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename))>0:
                readarray=np.loadtxt(filename)
                if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                self.airintervals_depthtop=readarray[:,0]
                self.airintervals_depthbot=readarray[:,1]
                self.airintervals_duration=readarray[:,2]
                self.airintervals_sigma=readarray[:,3]
            else:
                self.airintervals_depthtop=np.array([])
                self.airintervals_depthbot=np.array([])
                self.airintervals_duration=np.array([])
                self.airintervals_sigma=np.array([])

        filename=datadir+self.label+'/Ddepth.txt'
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename))>0:
                readarray=np.loadtxt(filename)
                if (np.size(readarray)==np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                self.Ddepth_depth=readarray[:,0]
                self.Ddepth_Ddepth=readarray[:,1]
                self.Ddepth_sigma=readarray[:,2]
            else:
                self.Ddepth_depth=np.array([])
                self.Ddepth_Ddepth=np.array([])
                self.Ddepth_sigma=np.array([])

        self.tuning_depth = {}
        self.tuning_proxy = {}
        self.tuning_proxy_sigma = {}

        self.tuning_age = {}
        self.tuning_target = {}
        self.tuning_target_sigma = {}
        if not hasattr(self,'tuning_uncertainty'):
            self.tuning_uncertainty={}

        if hasattr(self,'dict'):  # Tuning files
            for proxy, tag in self.dict.items():
                filename = datadir+self.label+'/' + proxy + '.txt'   # Tuning proxies
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename)) > 0:
                        readarray = np.loadtxt(filename)
                        if (np.size(readarray) == np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                        self.tuning_depth.update({proxy: readarray[:, 0]})
                        self.tuning_proxy.update({proxy: readarray[:, 1]})
                        self.tuning_proxy_sigma.update({proxy: readarray[:, 2]})

                    else:
                        raise ValueError('Tuning proxy file ' + element + '.txt not found in ' + datadir+self.label+'/')
                if not hasattr(self, 'tuning_uncertainty'):
                    self.tuning_uncertainty.update({proxy:0})

                filename = datadir + self.label + '/' + proxy + '_target.txt'  # Tuning targets
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename)) > 0:
                        readarray = np.loadtxt(filename)
                        if (np.size(readarray) == np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                        self.tuning_age.update({proxy: readarray[:, 0]})
                        self.tuning_target.update({proxy: readarray[:, 1]})
                        self.tuning_target_sigma.update({proxy: readarray[:, 2]})

                    else:
                        raise ValueError('Tuning target file ' + element + '_target.txt not found in ' + datadir+self.label+'/')
        else:
            self.dict = {}

        self.icemarkers_correlation=np.diag(np.ones(np.size(self.icemarkers_depth)))
        self.airmarkers_correlation=np.diag(np.ones(np.size(self.airmarkers_depth)))
        self.iceintervals_correlation=np.diag(np.ones(np.size(self.iceintervals_depthtop)))
        self.airintervals_correlation=np.diag(np.ones(np.size(self.airintervals_depthtop)))
        self.Ddepth_correlation=np.diag(np.ones(np.size(self.Ddepth_depth)))

        self.tuning_correlation = {}
        for proxy, tag in self.dict.items():
                self.tuning_correlation.update({proxy: np.diag(np.ones(np.size(self.tuning_proxy[proxy])))})

        filename=datadir+'/parameters-CovarianceObservations-AllDrillings.py'
        if os.path.isfile(filename):
            execfile(filename)
        filename=datadir+self.label+'/parameters-CovarianceObservations.py'
        if os.path.isfile(filename):
            execfile(filename)
        if np.size(self.icemarkers_depth)>0:
            self.icemarkers_chol=cholesky(self.icemarkers_correlation)
            self.icemarkers_lu_piv=scipy.linalg.lu_factor(np.transpose(self.icemarkers_chol))  #FIXME: we LU factor a triangular matrix. This is suboptimal. We should set lu_piv directly instead.
        if np.size(self.airmarkers_depth)>0:
            self.airmarkers_chol=cholesky(self.airmarkers_correlation)
            self.airmarkers_lu_piv=scipy.linalg.lu_factor(np.transpose(self.airmarkers_chol))
        if np.size(self.iceintervals_depthtop)>0:
            self.iceintervals_chol=cholesky(self.iceintervals_correlation)
            self.iceintervals_lu_piv=scipy.linalg.lu_factor(np.transpose(self.iceintervals_chol))
        if np.size(self.airintervals_depthtop)>0:
            self.airintervals_chol=cholesky(self.airintervals_correlation)
            self.airintervals_lu_piv=scipy.linalg.lu_factor(np.transpose(self.airintervals_chol))
        if np.size(self.Ddepth_depth)>0:
            self.Ddepth_chol=cholesky(self.Ddepth_correlation)
            self.Ddepth_lu_piv=scipy.linalg.lu_factor(np.transpose(self.Ddepth_chol))
        if np.size(self.tuning_proxy)>0:  # Tuning correlation matrix
            self.tuning_chol={}
            self.tuning_lu_piv={}
            for proxy, tag in self.dict.items():
                temp_chol = cholesky(self.tuning_correlation[proxy])
                self.tuning_chol.update({proxy: temp_chol})
                self.tuning_lu_piv.update({proxy: scipy.linalg.lu_factor(np.transpose(temp_chol))})


    def raw_model(self):



        #Accumulation
        if self.calc_a:
            self.a_model=self.A0*np.exp(self.beta*(self.deutice_fullcorr-self.deutice_fullcorr[0])) #Parrenin et al. (CP, 2007a) 2.3 (6)

        #Thinning
        if self.calc_tau:
            self.p=-1+m.exp(self.pprime)
            self.mu=m.exp(self.muprime)
#            self.s=m.tanh(self.sprime)
            omega_D=1-(self.p+2)/(self.p+1)*(1-self.zeta)+1/(self.p+1)*(1-self.zeta)**(self.p+2)	#Parrenin et al. (CP, 2007a) 2.2 (3)
            omega=self.s*self.zeta+(1-self.s)*omega_D   #Parrenin et al. (CP, 2007a) 2.2 (2)
            self.tau_model=(1-self.mu)*omega+self.mu 

        #udepth
        self.udepth_model=self.udepth_top+np.cumsum(np.concatenate((np.array([0]), self.D/self.tau_model*self.depth_inter)))
        
        g_model=interp1d(self.iedepth, self.udepth_model)
        self.LIDIE_model=self.LID_model*self.Dfirn
        self.ULIDIE_model=g_model(self.LIDIE_model)
        i_model=interp1d(self.udepth_model, self.depth)

        #Ice age
        self.icelayerthick_model=self.tau_model*self.a_model/self.D
        self.age_model=self.age_top+np.cumsum(np.concatenate((np.array([0]), self.D/self.tau_model/self.a_model*self.depth_inter)))
            
        f_model=interp1d(self.depth, self.age_model, bounds_error=False, fill_value=np.nan)

        #air age
        self.ice_equiv_depth_model=i_model(np.where(self.udepth_model-self.ULIDIE_model>self.udepth_top, self.udepth_model-self.ULIDIE_model, np.nan))  
        self.Ddepth_model=self.depth-self.ice_equiv_depth_model
        self.airage_model=f_model(self.ice_equiv_depth_model)
        self.airlayerthick_model=1/np.diff(self.airage_model)

    def corrected_model(self):

        self.correlation_corr_a_before=self.correlation_corr_a+0
        self.correlation_corr_LID_before=self.correlation_corr_LID+0
        self.correlation_corr_tau_before=self.correlation_corr_tau+0

        filename=datadir+'/parameters-CovariancePrior-AllDrillings.py'
        if os.path.isfile(filename):
            execfile(filename)
        filename=datadir+self.label+'/parameters-CovariancePrior.py'
        if os.path.isfile(filename):
            execfile(filename)

        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_a=cholesky(self.correlation_corr_a)
        if (self.correlation_corr_LID_before!=self.correlation_corr_LID).any():
            self.chol_LID=cholesky(self.correlation_corr_LID)
        if (self.correlation_corr_a_before!=self.correlation_corr_a).any():
            self.chol_tau=cholesky(self.correlation_corr_tau)


        #Accu
        #print np.shape(self.chol_a)
        #print np.shape(self.corr_a)
        corr=np.dot(self.chol_a,self.corr_a)*self.sigmap_corr_a
        j=interp1d_extrap(self.corr_a_age, corr)
        self.a=self.a_model*np.exp(j(self.age_model[:-1])) #FIXME: we should use mid-age and not age

        #Thinning
        h=interp1d(self.corr_tau_depth, np.dot(self.chol_tau,self.corr_tau)*self.sigmap_corr_tau)
        self.tau=self.tau_model*np.exp(h(self.depth_mid))
        self.udepth=self.udepth_top+np.cumsum(np.concatenate((np.array([0]), self.D/self.tau*self.depth_inter)))
        g=interp1d_extrap(self.iedepth, self.udepth)
        corr=np.dot(self.chol_LID,self.corr_LID)*self.sigmap_corr_LID
        j=interp1d_extrap(self.corr_LID_age, corr)
        self.LID=self.LID_model*np.exp(j(self.age_model))
        self.LIDIE=self.LID*self.Dfirn
        self.ULIDIE=g(self.LIDIE)
        i=interp1d(self.udepth, self.depth)

        #Ice age
        self.icelayerthick=self.tau*self.a/self.D
        self.age=self.age_top+np.cumsum(np.concatenate((np.array([0]), self.D/self.tau/self.a*self.depth_inter)))
        f=interp1d(self.depth,self.age, bounds_error=False, fill_value=np.nan)

        self.ice_equiv_depth=i(np.where(self.udepth-self.LIDIE>self.udepth_top, self.udepth-self.LIDIE, np.nan))
        self.Ddepth=self.depth-self.ice_equiv_depth
        self.airage=f(self.ice_equiv_depth)
        self.airlayerthick=1/np.diff(self.airage)


    def model(self, variables):
        index=0
#        if self.calc_a==True:
#            self.A0=variables[index]
#            self.beta=variables[index+1]
#            index=index+2
#        if self.calc_tau==True:
##            self.p=-1+m.exp(variables[index])
##            self.s=variables[index+1]
##            self.mu=variables[index+2]
##            index=index+3
#            self.pprime=variables[index]
#            self.muprime=variables[index+1]
#            index=index+2
        self.corr_tau=variables[index:index+np.size(self.corr_tau)]
        self.corr_a=variables[index+np.size(self.corr_tau):index+np.size(self.corr_tau)+np.size(self.corr_a)]
        self.corr_LID=variables[index+np.size(self.corr_tau)+np.size(self.corr_a): index+np.size(self.corr_tau)+np.size(self.corr_a)+np.size(self.corr_LID)]

        ##Raw model

        self.raw_model()

        ##Corrected model

        self.corrected_model()

        return np.concatenate((self.age,self.airage,self.Ddepth,self.a,self.tau,self.LID,self.icelayerthick,self.airlayerthick)) 


    def write_init(self):
        self.a_init=self.a
        self.LID_init=self.LID
        self.tau_init=self.tau
        self.icelayerthick_init=self.icelayerthick
        self.airlayerthick_init=self.airlayerthick
        self.age_init=self.age
        self.airage_init=self.airage
        self.Ddepth_init=self.Ddepth

    def fct_age(self, depth):
        f=interp1d(self.depth,self.age)
        return f(depth)

    def fct_age_init(self, depth):
        f=interp1d(self.depth,self.age_init)
        return f(depth)
   
    def fct_age_model(self, depth):
        f=interp1d(self.depth,self.age_model)
        return f(depth)
   
    def fct_airage(self, depth): #Fixme: should this be interp1d_extrap?
        f=interp1d(self.depth,self.airage)
        return f(depth)

    def fct_airage_init(self, depth): #Fixme: should this be interp1d_extrap? or interp1d?
        f=interp1d(self.depth,self.airage_init)
        return f(depth)

    def fct_airage_model(self, depth):
        f=interp1d(self.depth,self.airage_model)
        return f(depth)

    def fct_tuning_target(self,age,proxy):  # New for tuning
        f = interp1d_extrap(self.tuning_age[proxy], self.tuning_target[proxy])
        return f(age)

    def fct_tuning_sigma(self,age,proxy):
        f = interp1d_extrap(self.tuning_age[proxy], self.tuning_target_sigma[proxy])
        return f(age)

    def fct_Ddepth(self, depth):
        f=interp1d(self.depth,self.Ddepth)
        return f(depth)

    def residuals(self, variables):
        self.model(variables)

        if hasattr(self, 'tuning_correlation'):
            self.tuning_correlation_before=self.tuning_correlation.copy()

        filename=datadir+self.label+'/parameters-CovarianceTuning.py'
        if os.path.isfile(filename):
            execfile(filename)

        for proxy in self.dict.keys():  # Synchronization correlation matrix
            if not np.array_equal(self.tuning_correlation_before[proxy],self.tuning_correlation[proxy]):
                self.tuning_chol.update({proxy:np.array([])})
                self.tuning_lu_piv.update({proxy:np.array([])})
                if np.size(self.tuning_correlation[proxy]) > 0:
                    self.tuning_chol.update({proxy: cholesky(self.tuning_correlation[proxy])})
                    self.tuning_lu_piv.update({proxy: scipy.linalg.lu_factor(self.tuning_chol[proxy])})

        resi_corr_a=self.corr_a
        resi_corr_LID=self.corr_LID
        resi_corr_tau=self.corr_tau
        resi_age=(self.fct_age(self.icemarkers_depth)-self.icemarkers_age)/self.icemarkers_sigma
        if np.size(self.icemarkers_depth)>0:
            resi_age=scipy.linalg.lu_solve(self.icemarkers_lu_piv,resi_age)
        resi_airage=(self.fct_airage(self.airmarkers_depth)-self.airmarkers_age)/self.airmarkers_sigma
        if np.size(self.airmarkers_depth)>0:
            resi_airage=scipy.linalg.lu_solve(self.airmarkers_lu_piv,resi_airage)
        resi_iceint=(self.fct_age(self.iceintervals_depthbot)-self.fct_age(self.iceintervals_depthtop)-self.iceintervals_duration)/self.iceintervals_sigma
        if np.size(self.iceintervals_depthtop)>0:
            resi_iceint=scipy.linalg.lu_solve(self.iceintervals_lu_piv,resi_iceint)
        resi_airint=(self.fct_airage(self.airintervals_depthbot)-self.fct_airage(self.airintervals_depthtop)-self.airintervals_duration)/self.airintervals_sigma
        if np.size(self.airintervals_depthtop)>0:
            resi_airint=scipy.linalg.lu_solve(self.airintervals_lu_piv,resi_airint)
        resi_Ddepth=(self.fct_Ddepth(self.Ddepth_depth)-self.Ddepth_Ddepth)/self.Ddepth_sigma
        if np.size(self.Ddepth_depth)>0:
            resi_Ddepth=scipy.linalg.lu_solve(self.Ddepth_lu_piv,resi_Ddepth)

        resi_tuning = np.array([])
        if np.size(self.tuning_proxy)>0:
            for proxy, tag in self.dict.items():
                if tag == 'ice':
                    try:
                        resi_tuning_temp = (self.fct_tuning_target(self.fct_age(self.tuning_depth[proxy]), proxy)
                                        - self.tuning_proxy[proxy]) / np.sqrt(self.tuning_proxy_sigma[proxy]**2
                                                             + self.fct_tuning_sigma(self.fct_age(self.tuning_depth[proxy]),proxy)**2)
                        # correlation = np.corrcoef(self.fct_tuning_target(self.fct_age(self.tuning_depth[proxy]), proxy),self.tuning_proxy[proxy])[0,1]
                        # resi_tuning_temp *= (1-correlation)
                    except ValueError as e:
                        e.args += ('Orbital tuning proxy may be out of provided target age range. Try giving older target ages.',)
                        raise
                elif tag == 'air':
                    try:
                        partone = self.fct_airage(self.tuning_depth[proxy])
                        parttwo = self.fct_tuning_target(partone,proxy)
                        partfour = self.fct_tuning_sigma(partone,proxy)
                        self.partfive = parttwo-self.tuning_proxy[proxy]
                        self.partsix = np.sqrt(self.tuning_proxy_sigma[proxy]**2 + partfour**2 + self.tuning_uncertainty[proxy]**2)
                        resi_tuning_temp = self.partfive/self.partsix
                        # resi_tuning_temp = (self.fct_tuning_target(self.fct_airage(self.tuning_depth[proxy]), proxy)
                        #                 - self.tuning_proxy[proxy]) / np.sqrt(self.tuning_proxy_sigma[proxy]**2
                        #                                      + self.fct_tuning_sigma(self.fct_airage(self.tuning_depth[proxy]),proxy)**2)
                        # correlation = np.corrcoef(self.fct_tuning_target(self.fct_airage(self.tuning_depth[proxy]), proxy),self.tuning_proxy[proxy])[0,1]
                        # resi_tuning_temp *= (1-correlation)
                    except ValueError as e:
                        e.args += ('Orbital tuning proxy may be out of provided target age range. Try giving older target ages.',)
                        raise
                else: raise ValueError("'"+tag+
                                       "' age scale not recognized. Try 'ice' or 'air' in "+self.label+"/parameters.py")

                resi_tuning = np.append(resi_tuning,(scipy.linalg.lu_solve(self.tuning_lu_piv[proxy],resi_tuning_temp)))

        return np.concatenate((resi_corr_a, resi_corr_LID, resi_corr_tau, resi_age,resi_airage, resi_iceint, resi_airint, resi_Ddepth, resi_tuning))


    def cost_function(self):
        cost=np.dot(self.residuals,np.transpose(self.residuals))
        return cost

    def jacobian(self):
        epsilon=np.sqrt(np.diag(self.hess))/100000000.
        model0=self.model(self.variables)
        jacob=np.empty((np.size(model0), np.size(self.variables)))
        for i in np.arange(np.size(self.variables)):
            var=self.variables+0
            var[i]=var[i]+epsilon[i]
            model1=self.model(var)
            jacob[:,i]=(model1-model0)/epsilon[i]
        model0=self.model(self.variables)

        return jacob
    
    
    def optimisation(self) : 
        self.variables,self.hess,self.infodict,mesg,ier=leastsq(self.residuals, self.variables, full_output=1)
        print self.variables
        print self.hess
        return self.variables, self.hess
      
        
    def sigma(self):
        jacob=self.jacobian()

        index=0
        c_model=np.dot(jacob[index:index+np.size(self.age),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.age),:])))
        self.sigma_age=np.sqrt(np.diag(c_model))
        index=index+np.size(self.age)
        c_model=np.dot(jacob[index:index+np.size(self.airage),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.airage),:])))
        self.sigma_airage=np.sqrt(np.diag(c_model))
        index=index+np.size(self.airage)
        c_model=np.dot(jacob[index:index+np.size(self.Ddepth),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.Ddepth),:])))
        self.sigma_Ddepth=np.sqrt(np.diag(c_model))
        index=index+np.size(self.Ddepth)
        c_model=np.dot(jacob[index:index+np.size(self.a),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.a),:])))
        self.sigma_a=np.sqrt(np.diag(c_model))
        index=index+np.size(self.a)
        c_model=np.dot(jacob[index:index+np.size(self.tau),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.tau),:])))
        self.sigma_tau=np.sqrt(np.diag(c_model))
        index=index+np.size(self.tau)
        c_model=np.dot(jacob[index:index+np.size(self.LID),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.LID),:])))
        self.sigma_LID=np.sqrt(np.diag(c_model))
        index=index+np.size(self.LID)
        c_model=np.dot(jacob[index:index+np.size(self.icelayerthick),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.icelayerthick),:])))
        self.sigma_icelayerthick=np.sqrt(np.diag(c_model))
        index=index+np.size(self.icelayerthick)
        c_model=np.dot(jacob[index:index+np.size(self.airlayerthick),:],np.dot(self.hess,np.transpose(jacob[index:index+np.size(self.airlayerthick),:])))
        self.sigma_airlayerthick=np.sqrt(np.diag(c_model))

        f=interp1d_extrap(self.corr_a_age, self.sigmap_corr_a)
        self.sigma_a_model=f((self.age_model[1:]+self.age_model[:-1])/2)
        f=interp1d_extrap(self.corr_LID_age, self.sigmap_corr_LID)
        self.sigma_LID_model=f(self.age_model)
        f=interp1d_extrap(self.corr_tau_depth, self.sigmap_corr_tau)
        self.sigma_tau_model=f(self.depth_mid)

        
    

        
    def figures(self):

        mpl.figure(self.label+' thinning')
        mpl.title(self.label+' thinning')
        mpl.xlabel('Thinning')
        mpl.ylabel('Depth')
        if show_initial:
            mpl.plot(self.tau_init, self.depth_mid, color=color_init, label='Initial')
        mpl.plot(self.tau_model, self.depth_mid, color=color_mod, label='Prior')
        mpl.plot(self.tau, self.depth_mid, color=color_opt, label='Posterior +/-$\sigma$')
        mpl.fill_betweenx(self.depth_mid, self.tau-self.sigma_tau, self.tau+self.sigma_tau, color=color_ci)
#        mpl.plot(self.tau+self.sigma_tau, self.depth_mid, color='k', linestyle='-', label='+/- 1 sigma')
#        mpl.plot(self.tau-self.sigma_tau, self.depth_mid, color='k', linestyle='-')
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((x1,x2,self.depth[-1],self.depth[0]))
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/thinning.pdf')
        pp.savefig(mpl.figure(self.label+' thinning'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' ice layer thickness')
        mpl.title(self.label+' ice layer thickness')
        mpl.xlabel('thickness of annual layers (m/yr)')
        mpl.ylabel('Depth')
        if show_initial:
            mpl.plot(self.icelayerthick_init, self.depth_mid, color=color_init, label='Initial')
#        for i in range(np.size(self.iceintervals_duration)):
#            y1=self.iceintervals_depthtop[i]
#            y2=self.iceintervals_depthbot[i]
#            x1=(y2-y1)/(self.iceintervals_duration[i]+self.iceintervals_sigma[i])
#            x2=(y2-y1)/(self.iceintervals_duration[i]-self.iceintervals_sigma[i])
#            yserie=np.array([y1,y1,y2,y2,y1])
#            xserie=np.array([x1,x2,x2,x1,x1])
#            if i==0:
#                mpl.plot(xserie,yserie, color=color_obs, label="observations")
#            else:
#                mpl.plot(xserie,yserie, color=color_obs)
        mpl.plot(self.icelayerthick_model, self.depth_mid, color=color_mod, label='Prior')
        mpl.plot(self.icelayerthick, self.depth_mid, color=color_opt, label='Posterior +/-$\sigma$')
        mpl.fill_betweenx(self.depth_mid, self.icelayerthick-self.sigma_icelayerthick, self.icelayerthick+self.sigma_icelayerthick, color=color_ci)
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((0,x2,self.depth[-1],self.depth[0]))
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/icelayerthick.pdf')
        pp.savefig(mpl.figure(self.label+' ice layer thickness'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' air layer thickness')
        mpl.title(self.label+' air layer thickness')
        mpl.xlabel('thickness of annual layers (m/yr)')
        mpl.ylabel('Depth')
        if show_initial:
            mpl.plot(self.airlayerthick_init, self.depth_mid, color=color_init, label='Initial')
#        for i in range(np.size(self.airintervals_duration)):
#            y1=self.airintervals_depthtop[i]
#            y2=self.airintervals_depthbot[i]
#            x1=(y2-y1)/(self.airintervals_duration[i]+self.airintervals_sigma[i])
#            x2=(y2-y1)/(self.airintervals_duration[i]-self.airintervals_sigma[i])
#            yserie=np.array([y1,y1,y2,y2,y1])
#            xserie=np.array([x1,x2,x2,x1,x1])
#            if i==0:
#                mpl.plot(xserie,yserie, color=color_obs, label='observations')
#            else:
#                mpl.plot(xserie,yserie, color=color_obs)
        mpl.plot(self.airlayerthick_model, self.depth_mid, color=color_mod, label='Prior')
        mpl.plot(self.airlayerthick, self.depth_mid, color=color_opt, label='Posterior +/-$\sigma$')
        mpl.fill_betweenx(self.depth_mid, self.airlayerthick-self.sigma_airlayerthick, self.airlayerthick+self.sigma_airlayerthick, color=color_ci)
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((0, 2*max(self.icelayerthick),self.depth[-1],self.depth[0]))
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/airlayerthick.pdf')
        if show_airlayerthick:
            pp.savefig(mpl.figure(self.label+' air layer thickness'))  #Fixme: buggy line on anaconda
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' accumulation')
        mpl.title(self.label+' accumulation')
        mpl.xlabel('Optimized age (yr)')
        mpl.ylabel('Accumulation (m/yr)')
        if show_initial:
            mpl.step(self.age, np.concatenate((self.a_init, np.array([self.a_init[-1]]))), color=color_init, where='post', label='Initial')
        mpl.step(self.age, np.concatenate((self.a_model, np.array([self.a_model[-1]]))), color=color_mod, where='post', label='Prior')
        mpl.step(self.age, np.concatenate((self.a, np.array([self.a[-1]]))), color=color_opt, where='post', label='Posterior +/-$\sigma$')
        mpl.fill_between(self.age[:-1], self.a-self.sigma_a, self.a+self.sigma_a, color=color_ci)
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((self.age_top,x2,y1,y2))
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/accumulation.pdf')
        pp.savefig(mpl.figure(self.label+' accumulation'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' LID')
        mpl.title(self.label+' LID')
        mpl.xlabel('Optimized age (yr)')
        mpl.ylabel('LID')
        if show_initial:
            mpl.plot(self.age, self.LID_init, color=color_init, label='Initial')
        mpl.plot(self.age, self.LID_model, color=color_mod, label='Prior')
        mpl.plot(self.age, self.LID, color=color_opt, label='Posterior +/-$\sigma$')
        mpl.fill_between(self.age, self.LID-self.sigma_LID, self.LID+self.sigma_LID, color=color_ci)
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((self.age_top,x2,y1,y2))
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/LID.pdf')
        pp.savefig(mpl.figure(self.label+' LID'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' ice age')
        mpl.title(self.label+' ice age')
        mpl.xlabel('age (yr b1950)')
        mpl.ylabel('depth (m)')
        if show_initial:
            mpl.plot(self.age_init, self.depth, color=color_init, label='Initial')
        if (np.size(self.icemarkers_depth)>0):
            mpl.errorbar(self.icemarkers_age, self.icemarkers_depth, color=color_obs, xerr=self.icemarkers_sigma, linestyle='', marker='o', markersize=2, label="dated horizons")
#        mpl.ylim(mpl.ylim()[::-1])
        for i in range(np.size(self.iceintervals_duration)):
            y1=self.iceintervals_depthtop[i]
            y2=self.iceintervals_depthbot[i]
            x1=self.fct_age(y1)  #(y2-y1)/(self.iceintervals_duration[i]+self.iceintervals_sigma[i])
            x2=x1+self.iceintervals_duration[i]  #(y2-y1)/(self.iceintervals_duration[i]-self.iceintervals_sigma[i])
            xseries=np.array([x1,x2,x2,x1,x1])
            yseries=np.array([y1,y1,y2,y2,y1])
            if i==0:
                mpl.plot(xseries, yseries, color=color_di, label="dated intervals")
                mpl.errorbar(x2, y2, color=color_di, xerr=self.iceintervals_sigma[i], capsize=1)
            else:
                mpl.plot(xseries, yseries, color=color_di)
                mpl.errorbar(x2, y2, color=color_di, xerr=self.iceintervals_sigma[i], capsize=1)
#            mpl.arrow(x1, y1, x2-x1, y2-y1, fc=color_di, ec=color_di, head_width=0.02, head_length=0.05)        
#        if (np.size(self.iceintervals_depthtop)>0):
#            mpl.errorbar(self.fct_age(self.iceintervals_depthtop)+self.iceintervals_duration, self.iceintervals_depthbot, color=color_di, xerr=self.iceintervals_sigma, linestyle='', marker='o', markersize='2', label="dated intervals")
        mpl.plot(self.age_model, self.depth, color=color_mod, label='Prior')
        mpl.plot(self.age, self.depth, color=color_opt, label='Posterior +/-$\sigma$')
        mpl.fill_betweenx(self.depth, self.age-self.sigma_age, self.age+self.sigma_age , color=color_ci)
#        mpl.plot(self.age-self.sigma_age, self.depth, color='k', linestyle='-')
        mpl.plot(self.sigma_age*scale_ageci, self.depth, color=color_sigma, label='$\sigma$ x'+str(scale_ageci))   
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((self.age_top,x2,self.depth[-1],self.depth[0]))    
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/ice_age.pdf')
        pp.savefig(mpl.figure(self.label+' ice age'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' air age')
        mpl.title(self.label+' air age')
        mpl.xlabel('age (yr b1950)')
        mpl.ylabel('depth (m)')
        if show_initial:
            mpl.plot(self.airage_init, self.depth, color=color_init, label='Initial')
        if (np.size(self.airmarkers_depth)>0):
            mpl.errorbar(self.airmarkers_age, self.airmarkers_depth, color=color_obs, xerr=self.airmarkers_sigma, linestyle='', marker='o', markersize=2, label="observations")
#        mpl.ylim(mpl.ylim()[::-1])
        for i in range(np.size(self.airintervals_duration)):
            y1=self.airintervals_depthtop[i]
            y2=self.airintervals_depthbot[i]
            x1=self.fct_airage(y1)  #(y2-y1)/(self.iceintervals_duration[i]+self.iceintervals_sigma[i])
            x2=x1+self.airintervals_duration[i]  #(y2-y1)/(self.iceintervals_duration[i]-self.iceintervals_sigma[i])
            xseries=np.array([x1,x2,x2,x1,x1])
            yseries=np.array([y1,y1,y2,y2,y1])
            if i==0:
                mpl.plot(xseries, yseries, color=color_di, label="dated intervals")
                mpl.errorbar(x2, y2, color=color_di, xerr=self.airintervals_sigma[i], capsize=1)
            else:
                mpl.plot(xseries, yseries, color=color_di)
                mpl.errorbar(x2, y2, color=color_di, xerr=self.airintervals_sigma[i], capsize=1)
        mpl.plot(self.airage_model, self.depth, color=color_mod, label='Prior')
        mpl.fill_betweenx(self.depth, self.airage-self.sigma_airage, self.airage+self.sigma_airage , color=color_ci)
        mpl.plot(self.airage, self.depth, color=color_opt, label='Posterior +/-$\sigma$')
#        mpl.plot(self.airage+self.sigma_airage, self.depth, color='k', linestyle='-', label='+/- 1 sigma')
#        mpl.plot(self.airage-self.sigma_airage, self.depth, color='k', linestyle='-')
        mpl.plot(self.sigma_airage*scale_ageci, self.depth, color=color_sigma, label='$\sigma$ x'+str(scale_ageci))  
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((self.age_top,x2,self.depth[-1],self.depth[0]))    
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/air_age.pdf')
        pp.savefig(mpl.figure(self.label+' air age'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' Ddepth')
        mpl.title(self.label+' $\Delta$depth')
        mpl.xlabel('$\Delta$depth (m)')
        mpl.ylabel('Air depth (m)')
        if show_initial:
            mpl.plot(self.Ddepth_init, self.depth, color=color_init, label='Initial')
        if (np.size(self.Ddepth_depth)>0):
            mpl.errorbar(self.Ddepth_Ddepth, self.Ddepth_depth, color=color_obs, xerr=self.Ddepth_sigma, linestyle='', marker='o', markersize=2, label="observations")
#        mpl.ylim(mpl.ylim()[::-1])
        mpl.plot(self.Ddepth_model, self.depth, color=color_mod, label='Prior')
        mpl.plot(self.Ddepth, self.depth, color=color_opt, label='Posterior +/-$\sigma$')
        mpl.fill_betweenx(self.depth, self.Ddepth-self.sigma_Ddepth, self.Ddepth+self.sigma_Ddepth, color=color_ci)
#        mpl.plot(self.Ddepth+self.sigma_Ddepth, self.depth, color='k', linestyle='-', label='+/- 1 sigma')
#        mpl.plot(self.Ddepth-self.sigma_Ddepth, self.depth, color='k', linestyle='-')
        x1,x2,y1,y2 = mpl.axis()
        mpl.axis((x1,x2,self.depth[-1],self.depth[0]))
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/Ddepth.pdf')
        pp.savefig(mpl.figure(self.label+' Ddepth'))
        pp.close()
        if not show_figures:
            mpl.close()

        for proxy, tag in self.dict.items():
            fig = mpl.figure(self.label + ' ' + proxy+ 'black')
            ax = mpl.subplot(111)
            ax.patch.set_facecolor((0.02,0.02,0.02))
            #mpl.xlabel(tag + ' age (yr BP)')
            #mpl.ylabel(proxy)
            # Age at final depth, proxy value minus sigma, proxy value plus sigma
            ax.plot(self.tuning_age[proxy], self.tuning_target[proxy], color='white', label='Target')
            if tag == 'ice':
                ax.plot(self.fct_age_init(self.tuning_depth[proxy]),self.tuning_proxy[proxy],'-',color='grey',label='Prior')
                #ax.plot(self.fct_age(self.tuning_depth[proxy]), self.tuning_proxy[proxy], color='red', label='Posterior')
            elif tag == 'air':
                ax.plot(self.fct_airage_init(self.tuning_depth[proxy]), self.tuning_proxy[proxy],'-',color='grey', label='Prior')
                #ax.plot(self.fct_airage(self.tuning_depth[proxy]), self.tuning_proxy[proxy], color='red', label='Posterior')

            leg = mpl.legend(frameon=False,loc="best")
            for text in leg.get_texts():
                text.set_color("white")
            ax.tick_params(axis='x', colors='white', top='off')
            ax.tick_params(axis='y', colors='white', right='off')
            fig.patch.set_facecolor((0.02, 0.02, 0.02))
            pp = PdfPages(datadir + self.label + '/tuning' + proxy + 'black.pdf')
            pp.savefig(mpl.figure(self.label + ' ' + proxy + 'black'),facecolor=fig.get_facecolor(), edgecolor='none')
            pp.close()
            if not show_figures:
                mpl.close()

        for proxy, tag in self.dict.items():
            fig = mpl.figure(self.label + ' ' + proxy)
            ax = mpl.subplot(111)
            #mpl.xlabel(tag + ' age (yr BP)')
            #mpl.ylabel(proxy)
            # Age at final depth, proxy value minus sigma, proxy value plus sigma
            ax.plot(self.tuning_age[proxy], self.tuning_target[proxy], color='black', label='Target')
            if tag == 'ice':
                ax.plot(self.fct_age_init(self.tuning_depth[proxy]),self.tuning_proxy[proxy],'-',color='grey',label='Prior')
                #ax.plot(self.fct_age(self.tuning_depth[proxy]), self.tuning_proxy[proxy], color='red', label='Posterior')
            elif tag == 'air':
                ax.plot(self.fct_airage_init(self.tuning_depth[proxy]), self.tuning_proxy[proxy],'-', color='grey', label='Prior')
                #ax.plot(self.fct_airage(self.tuning_depth[proxy]), self.tuning_proxy[proxy], color='red', label='Posterior')

            leg = mpl.legend(frameon=False,loc="best")

            pp = PdfPages(datadir + self.label + '/tuning' + proxy + '.pdf')
            pp.savefig(mpl.figure(self.label + ' ' + proxy),facecolor=fig.get_facecolor(), edgecolor='none')
            pp.close()
            if not show_figures:
                mpl.close()


    def save(self):
        output=np.vstack((self.depth,self.age,self.sigma_age,self.airage,self.sigma_airage,np.append(self.a,self.a[-1]),np.append(self.sigma_a,self.sigma_a[-1]),np.append(self.tau,self.tau[-1]),np.append(self.sigma_tau,self.sigma_tau[-1]),self.LID,self.sigma_LID, self.Ddepth,self.sigma_Ddepth,np.append(self.a_model,self.a_model[-1]),np.append(self.sigma_a_model,self.sigma_a_model[-1]),np.append(self.tau_model,self.tau_model[-1]),np.append(self.sigma_tau_model,self.sigma_tau_model[-1]),self.LID_model,self.sigma_LID_model,np.append(self.icelayerthick,self.icelayerthick[-1]),np.append(self.sigma_icelayerthick,self.sigma_icelayerthick[-1]),np.append(self.airlayerthick,self.airlayerthick[-1]),np.append(self.sigma_airlayerthick,self.sigma_airlayerthick[-1])))
        with open(datadir+self.label+'/output.txt','w') as f:
            f.write('#depth\tage\tsigma_age\tair_age\tsigma_air_age\taccu\tsigma_accu\tthinning\tsigma_thinning\tLID\tsigma_LID\tDdepth\tsigma_Ddepth\taccu_model\tsigma_accu_model\tthinning_model\tsigma_thinning_model\tLID_model\tsigma_LID_model\ticelayerthick\tsigma_icelayerthick\tairlayerthick\tsigma_airlayerthick\n')
            np.savetxt(f,np.transpose(output), delimiter='\t')
        self.variables_x = []
        self.variables_x = np.concatenate((self.variables_x, self.corr_tau_depth, self.corr_a_age, self.corr_LID_age))
        np.savetxt(datadir+self.label+'/restart.txt',np.column_stack((len(self.corr_tau),len(self.corr_a),len(self.corr_LID), [np.transpose(self.variables)])))
        np.savetxt(datadir+self.label+'/restart_x.txt',np.column_stack((len(self.corr_tau),len(self.corr_a),len(self.corr_LID), [np.transpose(self.variables_x)])))
    
#    def udepth_save(self):
#        np.savetxt(datadir+self.label+'/udepth.txt',self.udepth)


class DrillingPair:

    def __init__(self, D1, D2):
        self.D1=D1
        self.D2=D2


    def init(self):
        self.label=self.D1.label+'-'+self.D2.label
        try:
            execfile(datadir+self.D1.label+'-'+self.D2.label+'/parameters.py')
            print('Automated synchronization for drilling pair '+self.D1.label+'-'+self.D2.label)
        except IOError: pass
#        print 'Initialization of drilling pair ',self.label

        if os.path.isdir(datadir+self.D1.label+'-'+self.D2.label): filename=datadir+self.D1.label+'-'+self.D2.label+'/ice_depth.txt'
        else: filename=datadir+self.D2.label+'-'+self.D1.label+'/ice_depth.txt'

        if os.path.isfile(filename) and open(filename).read():
            readarray=np.loadtxt(filename)
            self.iceicemarkers_depth1=readarray[:,0]
            self.iceicemarkers_depth2=readarray[:,1]
            self.iceicemarkers_sigma=readarray[:,2]
        else:
            self.iceicemarkers_depth1=np.array([])
            self.iceicemarkers_depth2=np.array([])
            self.iceicemarkers_sigma=np.array([])

        filename=datadir+self.D1.label+'-'+self.D2.label+'/air_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray=np.loadtxt(filename)
            self.airairmarkers_depth1=readarray[:,0]
            self.airairmarkers_depth2=readarray[:,1]
            self.airairmarkers_sigma=readarray[:,2]
        else:
            self.airairmarkers_depth1=np.array([])
            self.airairmarkers_depth2=np.array([])
            self.airairmarkers_sigma=np.array([])

        filename=datadir+self.D1.label+'-'+self.D2.label+'/iceair_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray=np.loadtxt(filename)
            self.iceairmarkers_depth1=readarray[:,0]
            self.iceairmarkers_depth2=readarray[:,1]
            self.iceairmarkers_sigma=readarray[:,2]
        else:
            self.iceairmarkers_depth1=np.array([])
            self.iceairmarkers_depth2=np.array([])
            self.iceairmarkers_sigma=np.array([])

        filename=datadir+self.D1.label+'-'+self.D2.label+'/airice_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray=np.loadtxt(filename)
            self.airicemarkers_depth1=readarray[:,0]
            self.airicemarkers_depth2=readarray[:,1]
            self.airicemarkers_sigma=readarray[:,2]
        else:
            self.airicemarkers_depth1=np.array([])
            self.airicemarkers_depth2=np.array([])
            self.airicemarkers_sigma=np.array([])

        if hasattr(self,'dict'):  # Synchronization files
            for record in self.D1, self.D2:
                record.synchro_depth = {}
                record.synchro_proxy = {}
                record.synchro_sigma = {}
                record.synchro_correlation = {}
                record.synchro_std_series = {}
                for proxy in self.dict.keys():
                    filename = datadir+record.label+'/' + proxy + '.txt'   # Synchronization proxies
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        if os.path.isfile(filename) and open(filename).read() and np.size(np.loadtxt(filename)) > 0:
                            readarray = np.loadtxt(filename)
                            if (np.size(readarray) == np.shape(readarray)[0]): readarray.resize(1, np.size(readarray))
                            record.synchro_depth.update({proxy: readarray[:, 0]})
                            record.synchro_proxy.update({proxy: readarray[:,1]})
                            record.synchro_sigma.update({proxy: readarray[:,2]})
                            if hasattr(self,'synchro_uncertainty'):
                                record.synchro_sigma.update({proxy: np.sqrt(record.synchro_sigma[proxy]**2 + self.synchro_uncertainty[proxy]**2)})
                            # record.synchro_std_series.update({proxy: np.std(record.synchro_proxy[proxy])})
                            # synchro_sigma = readarray[:,2]
                            # synchro_std_series = np.std(record.synchro_proxy[proxy])
                            # record.synchro_sigma.update({proxy: np.sqrt(synchro_sigma**2+synchro_std_series**2)})
                            record.synchro_correlation.update({proxy: np.diag(np.ones(np.size(self.D1.synchro_proxy[proxy]) + np.size(self.D1.synchro_proxy[proxy])))})
                            print (record.label+' '+proxy+' record read with '+str(len(record.synchro_depth[proxy]))+ ' values')
                        else:
                            raise ValueError('Synchronization proxy file ' + proxy + '.txt not found in ' + datadir+record.label+'/')
        else:
            self.dict = {}

        # self.ln1 = {}
        # self.ln2 = {}
        # self.axes = {}
        # self.text = {}
        # self.fignum = 1
        # for proxy in self.dict.keys():
        #     mpl.ion()
        #     self.ln1, = mpl.plot(self.D1.fct_airage_init(self.D1.synchro_depth[proxy]),
        #              self.D1.synchro_proxy[proxy], color=color_obs, label=self.D1.label)
        #     self.ln2, = mpl.plot(self.D2.fct_airage_init(self.D2.synchro_depth[proxy]),
        #              self.D2.synchro_proxy[proxy], color=color_mod, label=self.D2.label)
        #     mpl.xlim(max(min(self.D1.fct_airage_init(self.D1.synchro_depth[proxy])),
        #                  min(self.D2.fct_airage_init(self.D2.synchro_depth[proxy]))),
        #              min(max(self.D1.fct_airage_init(self.D1.synchro_depth[proxy])),
        #                  max(self.D2.fct_airage_init(self.D2.synchro_depth[proxy]))))
        #     mpl.legend(loc="best")
        #     self.axes =  mpl.gca()
        #     self.text= mpl.text(60000,800,self.fignum)
        #     self.axes.set_xlabel('Air age (yr BP)')
        #     self.axes.set_ylabel(proxy)
        #     mpl.show()

        self.iceicemarkers_correlation=np.diag(np.ones(np.size(self.iceicemarkers_depth1)))
        self.airairmarkers_correlation=np.diag(np.ones(np.size(self.airairmarkers_depth1)))
        self.iceairmarkers_correlation=np.diag(np.ones(np.size(self.iceairmarkers_depth1)))
        self.airicemarkers_correlation=np.diag(np.ones(np.size(self.airicemarkers_depth1)))

        filename=datadir+'/parameters-CovarianceObservations-AllDrillingPairs.py'
        if os.path.isfile(filename):
            execfile(filename)
        filename=datadir+self.label+'/parameters-CovarianceObservations.py'
        if os.path.isfile(filename):
            execfile(filename)
        if np.size(self.iceicemarkers_depth1)>0:
            self.iceicemarkers_chol=cholesky(self.iceicemarkers_correlation)
            self.iceicemarkers_lu_piv=scipy.linalg.lu_factor(self.iceicemarkers_chol)
        if np.size(self.airairmarkers_depth1)>0:
            self.airairmarkers_chol=cholesky(self.airairmarkers_correlation)
            self.airairmarkers_lu_piv=scipy.linalg.lu_factor(self.airairmarkers_chol)
        if np.size(self.iceairmarkers_depth1)>0:
            self.iceairmarkers_chol=cholesky(self.iceairmarkers_correlation)
            self.iceairmarkers_lu_piv=scipy.linalg.lu_factor(self.iceairmarkers_chol)
        if np.size(self.airicemarkers_depth1)>0:
            self.airicemarkers_chol=cholesky(self.airicemarkers_correlation)
            self.airicemarkers_lu_piv=scipy.linalg.lu_factor(self.airicemarkers_chol)

        for record in self.D1, self.D2:
            record.synchro_chol = {}
            record.synchro_lu_piv = {}
            for proxy in self.dict.keys(): # Synchronization correlation matrix
                if np.size(record.synchro_proxy) > 0:
                    record.synchro_chol.update({proxy: cholesky(record.synchro_correlation[proxy])})
                    record.synchro_lu_piv.update({proxy: scipy.linalg.lu_factor(record.synchro_chol[proxy])})

    def fct_synchro_proxy(self, age, record, proxy):
        f = interp1d_extrap(record.synchro_age[proxy], record.synchro_proxy[proxy])
        return f(age)

    def fct_synchro_sigma(self, age, record, proxy):
        f = interp1d_extrap(record.synchro_age[proxy], record.synchro_sigma[proxy])
        return f(age)

    def residuals(self):
        for record in self.D1, self.D2:
            record.synchro_age = {}

            for proxy, tag in self.dict.items():
                if tag == 'ice':
                    check0 = record.synchro_depth[proxy]
                    check00 = record.fct_age(check0)
                    record.synchro_age.update({proxy: check00})
                elif tag == 'air':
                    check0 = record.synchro_depth[proxy]
                    check00 = record.fct_age(check0)
                    record.synchro_age.update({proxy: check00})

        for record in self.D1, self.D2:
            if hasattr(record, 'synchro_correlation'):
                record.synchro_correlation_before=record.synchro_correlation.copy()
            for proxy in self.dict.keys(): #Fixme: This is probably not needed (number of points evaluated in residual remains the same...)
                record.synchro_correlation.update({proxy: np.diag(np.ones(np.size(record.synchro_age[proxy])))})

        filename=datadir+self.label+'/parameters-CovarianceTuning.py'
        if os.path.isfile(filename):
            execfile(filename)

        for record in self.D1, self.D2:
            for proxy in self.dict.keys():  # Synchronization correlation matrix
                if not np.array_equal(record.synchro_correlation_before[proxy],record.synchro_correlation[proxy]):
                    record.synchro_chol.update({proxy:np.array([])})
                    record.synchro_lu_piv.update({proxy:np.array([])})
                    if np.size(record.synchro_correlation[proxy]) > 0:
                        record.synchro_chol.update({proxy: cholesky(record.synchro_correlation[proxy])})
                        record.synchro_lu_piv.update({proxy: scipy.linalg.lu_factor(record.synchro_chol[proxy])})

        #TODO: Make sure size of correlation matrix and size of synchronization residuals are coherent...

        resi_iceice=(self.D1.fct_age(self.iceicemarkers_depth1)-self.D2.fct_age(self.iceicemarkers_depth2))/self.iceicemarkers_sigma
        if np.size(self.iceicemarkers_depth1)>0:
            resi_iceice=scipy.linalg.lu_solve(self.iceicemarkers_lu_piv,resi_iceice)
        a = self.D1.fct_airage(self.airairmarkers_depth1)
        b = self.D2.fct_airage(self.airairmarkers_depth2)
        resi_airair=(a-b)/self.airairmarkers_sigma
        if np.size(self.airairmarkers_depth1)>0:
            resi_airair=scipy.linalg.lu_solve(self.airairmarkers_lu_piv,resi_airair)
        resi_iceair=(self.D1.fct_age(self.iceairmarkers_depth1)-self.D2.fct_airage(self.iceairmarkers_depth2))/self.iceairmarkers_sigma
        if np.size(self.iceairmarkers_depth1)>0:
            resi_iceair=scipy.linalg.lu_solve(self.iceairmarkers_lu_piv,resi_iceair)
        resi_airice=(self.D1.fct_airage(self.airicemarkers_depth1)-self.D2.fct_age(self.airicemarkers_depth2))/self.airicemarkers_sigma
        if np.size(self.airicemarkers_depth1)>0:
            resi_airice=scipy.linalg.lu_solve(self.airicemarkers_lu_piv,resi_airice)

        resi_synchro = np.array([])

        for record, alternate in (self.D1, self.D2),(self.D2, self.D1):
            for proxy in self.dict.keys():  # Synchronization residuals
                if len(record.synchro_age[proxy])>0:
                    checkone = record.synchro_proxy[proxy]
                    checkone_a=record.synchro_age[proxy]
                    checktwo = self.fct_synchro_proxy(checkone_a,alternate,proxy)
                    checkthree = record.synchro_sigma[proxy]**2
                    checkfour = self.fct_synchro_sigma(checkone_a,alternate,proxy)**2
                    resi_synchro_temp = (checkone-checktwo) /\
                                        np.sqrt(checkthree+checkfour)

                    resi_synchro = np.append(resi_synchro, scipy.linalg.lu_solve(record.synchro_lu_piv[proxy], resi_synchro_temp))

        resi=np.concatenate((resi_iceice,resi_airair,resi_iceair,resi_airice,resi_synchro))

        # self.fignum += 1
        # print self.fignum
        # if opt_method !='leastsq_parallel':
        #     for proxy in self.dict.keys():
        #         self.ln1.set_xdata(self.D1.fct_airage_init(self.D1.synchro_depth[proxy]))
        #         self.ln1.set_ydata(self.D1.synchro_proxy[proxy])
        #         self.ln2.set_xdata(self.D2.fct_airage_init(self.D2.synchro_depth[proxy]))
        #         self.ln2.set_ydata(self.D2.synchro_proxy[proxy])
        #         self.text.set_text(self.fignum)
        #         mpl.pause(0.00000000001)

        # print("\n")
        # print np.dot(resi_synchro, resi_synchro)
        # print np.dot(resi,resi)
        return resi
    

    def figures(self):

        if not os.path.isdir(datadir+self.label):
            os.mkdir(datadir+self.label)


        mpl.figure(self.label+' ice-ice')
        mpl.xlabel(self.D1.label+' ice age (yr b1950)')
        mpl.ylabel(self.D2.label+' ice age (yr b1950)')
        if (np.size(self.iceicemarkers_depth1)>0):
            if show_initial:
                mpl.errorbar(self.D1.fct_age_init(self.iceicemarkers_depth1),self.D2.fct_age_init(self.iceicemarkers_depth2), color=color_init, xerr=self.iceicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Initial")
            mpl.errorbar(self.D1.fct_age_model(self.iceicemarkers_depth1),self.D2.fct_age_model(self.iceicemarkers_depth2), color=color_mod, xerr=self.iceicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Prior")
            mpl.errorbar(self.D1.fct_age(self.iceicemarkers_depth1),self.D2.fct_age(self.iceicemarkers_depth2), color=color_opt, xerr=self.iceicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Posterior")
        x1,x2,y1,y2 = mpl.axis()
        x1=self.D1.age_top
        y1=self.D2.age_top
        mpl.axis((x1,x2,y1,y2))
        range=np.array([max(x1,y1),min(x2,y2)])
        mpl.plot(range,range, color=color_obs, label='perfect agreement')
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/ice-ice.pdf')
        pp.savefig(mpl.figure(self.label+' ice-ice'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' air-air')
        mpl.xlabel(self.D1.label+' air age (yr b1950)')
        mpl.ylabel(self.D2.label+' air age (yr b1950)')
        if (np.size(self.airairmarkers_depth1)>0):
            if show_initial:
                mpl.errorbar(self.D1.fct_airage_init(self.airairmarkers_depth1),self.D2.fct_airage_init(self.airairmarkers_depth2), color=color_init, xerr=self.airairmarkers_sigma, linestyle='', marker='o', markersize=2, label="Initial")
            mpl.errorbar(self.D1.fct_airage_model(self.airairmarkers_depth1),self.D2.fct_airage_model(self.airairmarkers_depth2), color=color_mod, xerr=self.airairmarkers_sigma, linestyle='', marker='o', markersize=2, label="Prior")
            mpl.errorbar(self.D1.fct_airage(self.airairmarkers_depth1),self.D2.fct_airage(self.airairmarkers_depth2), color=color_opt, xerr=self.airairmarkers_sigma, linestyle='', marker='o', markersize=2, label="Posterior")
        x1,x2,y1,y2 = mpl.axis()
        x1=self.D1.age_top
        y1=self.D2.age_top
        mpl.axis((x1,x2,y1,y2))
        range=np.array([max(x1,y1),min(x2,y2)])
        mpl.plot(range,range, color=color_obs, label='perfect agreement')
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/air-air.pdf')
        pp.savefig(mpl.figure(self.label+' air-air'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' ice-air')
        mpl.xlabel(self.D1.label+' ice age (yr b1950)')
        mpl.ylabel(self.D2.label+' air age (yr b1950)')
        if (np.size(self.iceairmarkers_depth1)>0):
            if show_initial:
                mpl.errorbar(self.D1.fct_age_init(self.iceairmarkers_depth1),self.D2.fct_airage_init(self.iceairmarkers_depth2), color=color_init, xerr=self.iceairmarkers_sigma, linestyle='', marker='o', markersize=2, label="Initial")
            mpl.errorbar(self.D1.fct_age_model(self.iceairmarkers_depth1),self.D2.fct_airage_model(self.iceairmarkers_depth2), color=color_mod, xerr=self.iceairmarkers_sigma, linestyle='', marker='o', markersize=2, label="Prior")
            mpl.errorbar(self.D1.fct_age(self.iceairmarkers_depth1),self.D2.fct_airage(self.iceairmarkers_depth2), color=color_opt, xerr=self.iceairmarkers_sigma, linestyle='', marker='o', markersize=2, label="Posterior")
        x1,x2,y1,y2 = mpl.axis()
        x1=self.D1.age_top
        y1=self.D2.age_top
        mpl.axis((x1,x2,y1,y2))
        range=np.array([max(x1,y1),min(x2,y2)])
        mpl.plot(range,range, color=color_obs, label='perfect agreement')
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/ice-air.pdf')
        pp.savefig(mpl.figure(self.label+' ice-air'))
        pp.close()
        if not show_figures:
            mpl.close()

        mpl.figure(self.label+' air-ice')
        mpl.xlabel(self.D1.label+' air age (yr b1950)')
        mpl.ylabel(self.D2.label+' ice age (yr b1950)')
        if (np.size(self.airicemarkers_depth1)>0):
            if show_initial:
                mpl.errorbar(self.D1.fct_airage_init(self.airicemarkers_depth1),self.D2.fct_age_init(self.airicemarkers_depth2), color=color_init, xerr=self.airicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Initial")
            mpl.errorbar(self.D1.fct_airage_model(self.airicemarkers_depth1),self.D2.fct_age_model(self.airicemarkers_depth2), color=color_mod, xerr=self.airicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Prior")
            mpl.errorbar(self.D1.fct_airage(self.airicemarkers_depth1),self.D2.fct_age(self.airicemarkers_depth2), color=color_opt, xerr=self.airicemarkers_sigma, linestyle='', marker='o', markersize=2, label="Posterior")
        x1,x2,y1,y2 = mpl.axis()
        x1=self.D1.age_top
        y1=self.D2.age_top
        mpl.axis((x1,x2,y1,y2))
        range=np.array([max(x1,y1),min(x2,y2)])
        mpl.plot(range,range, color=color_obs, label='perfect agreement')
        mpl.legend(loc="best")
        pp=PdfPages(datadir+self.label+'/air-ice.pdf')
        pp.savefig(mpl.figure(self.label+' air-ice'))
        pp.close()
        if not show_figures:
            mpl.close()

        if hasattr(self,'dict'): #Fixme: needs to be generalized for air and ice proxies
            for proxy, tag in self.dict.items():
                mpl.figure(self.label+' '+proxy)
                mpl.xlabel('Air age (yr BP)')
                mpl.ylabel(proxy)
                #Age at final depth, proxy value minus sigma, proxy value plus sigma
                mpl.fill_between(self.D1.fct_airage(self.D1.synchro_depth[proxy]),
                                 self.D1.synchro_proxy[proxy] - self.D1.synchro_sigma[proxy],
                                 self.D1.synchro_proxy[proxy] + self.D1.synchro_sigma[proxy], color=color_ci)
                mpl.fill_between(self.D2.fct_airage(self.D2.synchro_depth[proxy]),
                                 self.D2.synchro_proxy[proxy] - self.D2.synchro_sigma[proxy],
                                 self.D2.synchro_proxy[proxy] + self.D2.synchro_sigma[proxy], color=color_ci)
                mpl.plot(self.D1.fct_airage(self.D1.synchro_depth[proxy]),
                                 self.D1.synchro_proxy[proxy], color=color_obs, label=self.D1.label)
                mpl.plot(self.D2.fct_airage(self.D2.synchro_depth[proxy]),
                         self.D2.synchro_proxy[proxy], color=color_mod, label=self.D2.label)
                mpl.xlim(max(np.nanmin(self.D1.fct_airage(self.D1.synchro_depth[proxy])),np.nanmin(self.D2.fct_airage(self.D2.synchro_depth[proxy]))),
                         min(np.nanmax(self.D1.fct_airage(self.D1.synchro_depth[proxy])),np.nanmax(self.D2.fct_airage(self.D2.synchro_depth[proxy]))))
                mpl.legend(loc="best")
                pp = PdfPages(datadir + self.label + '/synchro' + proxy + '.pdf')
                pp.savefig(mpl.figure(self.label+' '+proxy))
                pp.close()
                if not show_figures:
                    mpl.close()

                for proxy, tag in self.dict.items():
                    mpl.figure(self.label + ' ' + proxy)
                    mpl.xlabel('Air age (yr BP)')
                    mpl.ylabel(proxy)
                    # Age at final depth, proxy value minus sigma, proxy value plus sigma
                    mpl.fill_between(self.D1.fct_airage_init(self.D1.synchro_depth[proxy]),
                                     self.D1.synchro_proxy[proxy] - self.D1.synchro_sigma[proxy],
                                     self.D1.synchro_proxy[proxy] + self.D1.synchro_sigma[proxy], color=color_ci)
                    mpl.fill_between(self.D2.fct_airage_init(self.D2.synchro_depth[proxy]),
                                     self.D2.synchro_proxy[proxy] - self.D2.synchro_sigma[proxy],
                                     self.D2.synchro_proxy[proxy] + self.D2.synchro_sigma[proxy], color=color_ci)
                    mpl.plot(self.D1.fct_airage_init(self.D1.synchro_depth[proxy]),
                             self.D1.synchro_proxy[proxy], color=color_obs, label=self.D1.label)
                    mpl.plot(self.D2.fct_airage_init(self.D2.synchro_depth[proxy]),
                             self.D2.synchro_proxy[proxy], color=color_mod, label=self.D2.label)
                    mpl.xlim(max(min(self.D1.fct_airage_init(self.D1.synchro_depth[proxy])),
                                 min(self.D2.fct_airage_init(self.D2.synchro_depth[proxy]))),
                             min(max(self.D1.fct_airage_init(self.D1.synchro_depth[proxy])),
                                 max(self.D2.fct_airage_init(self.D2.synchro_depth[proxy]))))
                    mpl.legend(loc="best")
                    pp = PdfPages(datadir + self.label + '/synchro' + proxy + '-prior.pdf')
                    pp.savefig(mpl.figure(self.label + ' ' + proxy))
                    pp.close()
                    if not show_figures:
                        mpl.close()