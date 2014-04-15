#Parameters specific to the Vostok ice core
self.dim=2
self.udepth_top=0.
self.calc_a=False
self.calc_tau=False
self.calc_LID=False
self.age_top=-47
self.depth=np.arange(0., 3501.+0.01, 1.)
self.thickness=3767.
#self.gamma_source=3.4
#self.beta_source=1.5
self.lambda_a=4000
self.sigmap_corr_a=0.2
self.k=0.1
self.lambda_tau=70
self.lambda_LID=4000
self.sigmap_corr_LID=0.4
self.age_bot=800000.+self.age_top
self.corr_a_age=np.arange(-50, 800000-50+0.01, self.age_step)
self.corr_LID_age=np.arange(-50, 800000-50+0.01, self.age_step)
self.LID_value=98.   #TODO: change this EDC value
self.Dfirn=0.698    #This is the EDC value
self.restart=False


#self.pprime=1.19
#self.s=0.
#self.mu=0.00066/self.A0
#self.beta=0.0157
#self.A0=0.02841
