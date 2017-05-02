#A covariance matrix for tuning...
# f = interp1d_extrap([0,100],[1,0])
# if hasattr(self,'synchro_age_eval'):
#     for proxy in self.dict.keys():
#         dummy_matrix = np.abs(
#             np.ones((np.size(self.synchro_age_eval[proxy]), np.size(self.synchro_age_eval[proxy])))
#             * self.synchro_age_eval[proxy] - np.transpose(np.ones((np.size(self.synchro_age_eval[proxy]), np.size(self.synchro_age_eval[proxy]))))
#             * self.synchro_age_eval[proxy])
#
#         self.synchro_correlation.update({proxy: np.diag(np.ones(np.size(self.synchro_age_eval[proxy])))})