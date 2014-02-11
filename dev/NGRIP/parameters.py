#Parameters specific to the NGRIP ice core
self.dim=1
self.calc_a=False
self.calc_tau=False
self.calc_LID=False
self.age_min=-30.
self.depth_min=8.
self.depth_max=3084.
self.step=1
self.lambda_a=4000
self.sigmap_corr_a=0.2
self.thickness=3085. #From NGRIP community members (2004)
self.k=0.1
self.lambda_tau=70
self.lambda_LID=4000
self.sigmap_corr_LID=0.4
self.age_max=150000.+self.age_min
self.corr_a=np.zeros(16)
self.corr_LID=np.zeros(16)
self.Dfirn=0.698    #From Parrenin et al. (CP, 2012b)
self.corr_tau=np.zeros(101)
self.restart=False
