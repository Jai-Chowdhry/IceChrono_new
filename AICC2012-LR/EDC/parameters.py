#Parameters specific to the EDC ice core
#TODO: Check these parameters with Parrenin et al. (CP, 2007), especially ice thickness
self.dim=1
self.calc_a=False
self.calc_tau=False
self.calc_LID=False
self.calc_udepth_init=False
self.udepth_min=0.
self.age_min=-55.
self.depth_min=0.
self.depth_max=3259.3
self.step=0.55
self.gamma_source=3.4
self.beta_source=1.5
self.A0=3.30e-02
self.lambda_a=4000
self.beta=1.65e-02
self.sigmap_corr_a=0.2
self.thickness=3273.
self.pprime=1.59
self.s=-2.82e-01
self.mu=5.34e-02
self.k=0.1
self.lambda_tau=70
self.lambda_LID=4000
self.sigmap_corr_LID=0.4
self.age_max=1000000.+self.age_min
self.corr_a=np.zeros(101)
self.corr_LID=np.zeros(101)
self.LID_value=98.
self.Dfirn=0.698    #From Parrenin et al. (CP, 2012b)
self.corr_tau=np.zeros(51)
self.restart=False
