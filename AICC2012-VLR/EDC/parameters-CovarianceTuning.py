# #A covariance matrix for orbital tuning...
# f = interp1d_extrap([0,40000],[1,0])
# if hasattr(self,'tuning_depth'):
#     for proxy in self.dict.keys():
#         dummy_matrix = np.abs(
#             np.ones((np.size(self.tuning_depth[proxy]), np.size(self.tuning_depth[proxy])))
#             * self.fct_age(self.tuning_depth[proxy]) - np.transpose(np.ones((np.size(self.tuning_depth[proxy]), np.size(self.tuning_depth[proxy])))
#             * self.fct_age(self.tuning_depth[proxy])))
#
#         # self.synchro_correlation.update({proxy: np.diag(np.ones(np.size(self.synchro_age_eval[proxy])))})
#         self.tuning_correlation.update({proxy: f(dummy_matrix)})