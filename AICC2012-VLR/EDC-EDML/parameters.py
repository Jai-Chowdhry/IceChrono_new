#Parameters specific to the EDC-EDML drilling pair.

self.dict = {'CH4':'air'}  # [Inter-core] tuning proxies and respective age scales, 'ice' or 'air'.
self.synchro_uncertainty = {'CH4':0}  # Modeling uncertainty for tuning proxy. Depends on confidence with respect to rest of residuals. Synchronization may bug if too small...