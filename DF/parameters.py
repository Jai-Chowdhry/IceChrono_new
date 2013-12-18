#Parameters specific to the DF ice core
#TODO: check this file with the Parrenin et al. (CP, 2007) parameters
self.dim=1
self.calc_a=True
self.calc_tau=True
self.age_surf=-50
self.max_depth=3032.
self.step=1.
self.gamma_source=3.2 
self.beta_source=1.6  
self.A0=0.02913
self.lambda_a=4000
self.beta=0.0147
self.sigmap_corr_a=0.2
self.thickness=3032.
self.pprime=1.19
self.s=0.
self.mu=0.00011/self.A0
self.k=0.1
self.lambda_tau=70
self.lambda_LIDIE=4000
self.sigmap_corr_LIDIE=0.4
self.age_max=1000000.
self.corr_a=np.zeros(101)
self.corr_LIDIE=np.zeros(101)
self.LID_depth=np.array([0., self.max_depth])
self.LID_LID=np.array([98., 98.])
self.corr_tau=np.zeros(101)
self.restart=False