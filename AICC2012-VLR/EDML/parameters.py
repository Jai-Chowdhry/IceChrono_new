#Parameters specific to the EDML ice core
self.udepth_top=8.54025
self.age_top=70.
self.depth=np.arange(18., 2564.+0.01, 1.)
self.corr_a_age=np.arange(self.age_top, 300000+self.age_top+0.01, self.age_step)
self.corr_LID_age=np.arange(self.age_top, 300000+self.age_top+0.01, self.age_step)
self.corr_tau_depth=np.arange(self.depth[0], self.depth[-1]+0.01, (self.depth[-1]-self.depth[0])/(self.corr_tau_nodes-1))
self.tuning_dict = {'CH4':'air'}  # [Orbital] tuning proxies and respective age scales, 'ice' or 'air'.
self.tuning_uncertainty = {'CH4':5}  # Modeling uncertainty for tuning proxy. Depends on confidence with respect to rest of residuals. Synchronization may bug if too small...
self.tuning_multi = {'CH4':True} #Which proxies to compare with other cores??
self.tuning_matrices = 'sparse'  # Specify if tuning correlation matrices are sparse, greatly improves computational efficiency for large matrices. Options 'sparse', 'dense'
self.calc_corr_tuning = False # Calculate the correlation matrix for tuning ? If false the correlation matrix is assumed identity (can speed up calculation)
self.restart=False


#self.thickness=3000
#self.cT2=0.000078
#self.sigmabA=0.5
#self.cA1=1.
#self.sigmabL=0.6
