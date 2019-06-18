from schwimmbad import MPIPool
import numpy as np
import time
import sys

def complist(counter):
    variables = np.random.multivariate_normal(np.zeros(1000), np.diag(np.ones(1000)))
    return variables

def initial(variables):
    fit = variables + np.random.multivariate_normal(np.zeros(np.size(variables)),
                                                    np.diag(np.ones(np.size(variables)))) * 0.02
    print('Initialized once ')
    return fit



pool = MPIPool()

if not pool.is_master():
    pool.wait()
    sys.exit(0)

print('Starting at '+str(time.time()))
variables = pool.map(complist, range(0,100))
values = np.array(pool.map(initial, variables))
print('Finishing at '+str(time.time()))

pool.close()