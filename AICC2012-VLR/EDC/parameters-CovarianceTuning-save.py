#CH4 Tuning correlation matrix
self.tuning_correlation = {}
for proxy in self.tuning_dict:
    self.tuning_correlation.update({proxy:np.interp(np.abs(np.ones((np.size(self.tuning_depth[proxy]),np.size(self.tuning_depth[proxy])))*np.linspace(min(self.fct_airage(self.tuning_depth[proxy])),max(self.fct_airage(self.tuning_depth[proxy])),np.size(self.tuning_depth[proxy]))-np.transpose(np.ones((np.size(self.tuning_depth[proxy]),np.size(self.tuning_depth[proxy])))*np.linspace(min(self.fct_airage(self.tuning_depth[proxy])),max(self.fct_airage(self.tuning_depth[proxy])),np.size(self.tuning_depth[proxy])))), np.array([0,self.lambda_tuning]),np.array([1, 0]))})
