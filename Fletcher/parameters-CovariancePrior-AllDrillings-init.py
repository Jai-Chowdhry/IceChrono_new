#Accu correlation matrix
#f=interp1d(np.array([0,self.lambda_a,10000000]),np.array([1, 0, 0]))
self.correlation_corr_a=np.interp(np.abs(np.ones((np.size(self.corr_a_age),np.size(self.corr_a_age)))*self.corr_a_age-np.transpose(np.ones((np.size(self.corr_a_age),np.size(self.corr_a_age)))*self.corr_a_age)),np.array([0,self.lambda_a,10000000]),np.array([1, 0, 0]))
#Gaussian shape, does not work for a too high resolution
#self.correlation_corr_a=gaussian(np.abs(np.ones((np.size(self.corr_a_age),np.size(self.corr_a_age)))*self.corr_a_age-np.transpose(np.ones((np.size(self.corr_a_age),np.size(self.corr_a_age)))*self.corr_a_age))/self.lambda_a)

#LID correlation matrix
#f=interp1d(np.array([0,self.lambda_LID,10000000]),np.array([1, 0, 0]))
self.correlation_corr_LID=np.interp(np.abs(np.ones((np.size(self.corr_LID_age),np.size(self.corr_LID_age)))*self.corr_LID_age-np.transpose(np.ones((np.size(self.corr_LID_age),np.size(self.corr_LID_age)))*self.corr_LID_age)),np.array([0,self.lambda_LID,10000000]),np.array([1, 0, 0]))
#Gaussian shape, does not work for a too high resolution
#self.correlation_corr_LID=gaussian(np.abs(np.ones((np.size(self.corr_LID_age),np.size(self.corr_LID_age)))*self.corr_LID_age-np.transpose(np.ones((np.size(self.corr_LID_age),np.size(self.corr_LID_age)))*self.corr_LID_age))/self.lambda_LID)

#Thinning correlation matrix
#g=interp1d(np.array([0,self.lambda_tau,5000]),np.array([1, 0, 0]))
self.correlation_corr_tau=np.interp(np.abs(np.ones((np.size(self.corr_tau_depth),np.size(self.corr_tau_depth)))*self.corr_tau_depth-np.transpose(np.ones((np.size(self.corr_tau_depth),np.size(self.corr_tau_depth)))*self.corr_tau_depth)),np.array([0,self.lambda_tau,5000]),np.array([1, 0, 0]))

#        f=interpolate.interp1d(self.depth,self.udepth_init)
#self.thickness_ie=self.thickness-self.depth[-1]+self.iedepth[-1]
#f=interp1d(self.depth,self.iedepth, bounds_error=False, fill_value=self.iedepth[-1]) #We should not need the bounds_error option. Check what is the problem.
#self.sigmap_corr_tau=self.k/self.thickness_ie*f(self.corr_tau_depth)

#f=interp1d(self.depth,self.udepth_model, bounds_error=False, fill_value=self.udepth_model[-1]) #We should not need the bounds_error option. Check what is the problem.
self.sigmap_corr_tau=self.k/self.thickness_ie*np.interp(self.corr_tau_depth,self.depth,self.udepth_model)


self.sigmap_corr_a=self.sigmap_corr_a*np.ones(np.size(self.corr_a_age))
self.sigmap_corr_LID=self.sigmap_corr_LID*np.ones(np.size(self.corr_LID_age))
