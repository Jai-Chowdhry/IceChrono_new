
# -*- coding: utf-8 -*-

#Modification of emcee implementation to allow for gamma = 1 every 10 iterations. (Highly multimodal problem)
#Could adjust gamma as done by Nelson et al., would need to pass acceptance fraction for each generation.

from __future__ import division, print_function

import numpy as np
from emcee.moves.red_blue import RedBlueMove
from emcee.backends import HDFBackend

__all__ = ["DEMove_modified"]


class DEMove_modified(RedBlueMove):
    """A proposal using differential evolution.
    This `Differential evolution proposal
    <http://www.stat.columbia.edu/~gelman/stuff_for_blog/cajo.pdf>`_ is
    implemented following `Nelson et al. (2013)
    <http://arxiv.org/abs/1311.5229>`_.
    Args:
        sigma (float): The standard deviation of the Gaussian used to stretch
            the proposal vector.
        gamma0 (Optional[float]): The mean stretch factor for the proposal
            vector. By default, it is :math:`2.38 / \sqrt{2\,\mathrm{ndim}}`
            as recommended by the two references.
    """
    def __init__(self, sigma=1.0e-5, gamma0=None, **kwargs):
        self.sigma = sigma
        self.gamma0 = gamma0
        self.counter = 0
        kwargs["nsplits"] = 3
        super(DEMove_modified, self).__init__(**kwargs)

    def setup(self, coords):
        self.g0 = self.gamma0
        if self.g0 is None:
            # Pure MAGIC:
            ndim = coords.shape[1]
            self.g0 = 2.38 / np.sqrt(2 * ndim)

    def get_proposal(self, s, c, random):
        Ns = len(s)
        Nc = list(map(len, c))
        ndim = s.shape[1]
        q = np.empty((Ns, ndim), dtype=np.float64)
        f = self.sigma * random.randn(Ns)
        if self.counter%25 == 0:
            self.g = 1.0
        else:
            self.g = self.g0
        for i in range(Ns):
            w = np.array([c[j][random.randint(Nc[j])] for j in range(2)])
            random.shuffle(w)
            g = np.diff(w, axis=0) * self.g + f[i]
            q[i] = s[i] + g
        self.counter+=1
        return q, np.zeros(Ns, dtype=np.float64)

    def tune(self,state,accepted):
        fraction = np.count_nonzero(accepted)/np.size(accepted)
        #print("Fraction is", fraction)
        if fraction > 0.31:
            self.g0 = self.g0 * 1.1
        elif fraction < 0.2:
            self.g0 = self.g0 * 0.9
        else:
            self.g0 = self.g0 * np.sqrt(fraction / 0.25)

class TerBraak2008_DEMove(RedBlueMove): #Fixme is this really the snooker move?

    def __init__(self, sigma=1.0e-5, gamma0=None, backend = None, **kwargs):
        self.sigma = sigma
        self.gamma0 = gamma0
        self.counter = 0
        kwargs["nsplits"] = 3
        super(TerBraak2008_DEMove, self).__init__(**kwargs)
        if backend == None:
            self.backend = HDFBackend('saved_iterations.hdf5', name='mcmc', read_only=True)
        else: self.backend = backend

    def setup(self, coords):
        self.g0 = self.gamma0
        if self.g0 is None:
            # Pure MAGIC:
            ndim = coords.shape[1]
            self.g0 = 2.38 / np.sqrt(2 * ndim)

    def get_proposal(self, s, c, random):

        backend = self.backend.get_chain()
        Ns = len(s)
        ndim = s.shape[1]
        q = np.empty((Ns, ndim), dtype=np.float64)
        f = self.sigma * random.randn(Ns)

        for i in range(Ns):
            w = np.array([backend[random.randint(0,np.shape(backend)[0])][random.randint(0,np.shape(backend)[1])] for j in range(2)])
            random.shuffle(w)
            g = np.diff(w, axis=0) * self.g0 + f[i]
            q[i] = s[i] + g
        self.counter+=1
        if True in np.isnan(q[i]):
            print('Warning : the adaptive DE move gave NaN')

        return q, np.zeros(Ns, dtype=np.float64)

class TerBraak2008_snookermove(RedBlueMove):

    def __init__(self, sigma=1.0e-5, gammas=1.7, backend = None, **kwargs):
        self.sigma = sigma
        self.gammas = gammas
        self.counter = 0
        kwargs["nsplits"] = 4
        super(TerBraak2008_snookermove, self).__init__(**kwargs)
        if backend == None:
            self.backend = HDFBackend('saved_iterations.hdf5', name='mcmc', read_only=True)
        else: self.backend = backend

    def get_proposal(self, s, c, random):

        backend = self.backend.get_chain()
        Ns = len(s)
        ndim = s.shape[1]
        q = np.empty((Ns, ndim), dtype=np.float64)
        metropolis = np.empty(Ns, dtype=np.float64)

        for i in range(Ns):
            w = np.array([backend[random.randint(0,np.shape(backend)[0])][random.randint(0,np.shape(backend)[1])] for j in range(3)])
            random.shuffle(w)
            z, z1, z2 = w
            delta = s[i] - z
            norm = np.linalg.norm(delta)
            if np.sqrt(norm) == 0: #Avoiding divide by zero
                u = np.ones(np.shape(delta)) / np.sqrt(2 * ndim) #Change to DE move
            else:
                u = delta/np.sqrt(norm)
            q[i] = s[i] + u*self.gammas*(np.dot(u,z1)-np.dot(u,z2))
            metropolis[i] =np.log(np.linalg.norm(q[i]-z)) -np.log(norm)
            if True in np.isnan(q[i]):
                print('Warning : the adaptive snooker move gave NaN')


        return q, 0.5 * (ndim-1.0) * metropolis