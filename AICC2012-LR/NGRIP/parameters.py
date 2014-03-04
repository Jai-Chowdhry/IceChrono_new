#Parameters specific to the NGRIP ice core
self.dim=1
self.calc_a=False
self.calc_tau=False
self.calc_LID=False
self.calc_udepth=True
self.udepth_min=3.313966267738
self.age_min=-30.
self.depth_min=8.
self.depth_max=3084.
self.step=1
self.thickness=3085. #From NGRIP community members (2004)
self.age_max=150000.+self.age_min
self.corr_a=np.zeros(16)
self.corr_LID=np.zeros(16)
self.Dfirn=0.698    #From Parrenin et al. (CP, 2012b)
self.restart=False
