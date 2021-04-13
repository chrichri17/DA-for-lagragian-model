# pylint: disable=import-error, unbalanced-tuple-unpacking

import numpy as np
import numpy.random as random
import numpy.linalg as la
from copy import deepcopy
from Lagrange import Lagrange

class Ensemble:
    
    def __init__(self):
        self.list = []
    
    def __iter__(self):
        return iter(self.list)
    
    def __getitem__(self, key):
        return self.list[key]
    
    def append(self, member):
        self.list.append(member)
    
    def propagate(self, steps):
        for member in self.list:
            member.propagate(steps)

class Member:

    def __init__(self, grid, NSTEPS=400):
        self.model = Lagrange(grid, NSTEPS)
        self.model.start_fire(ITYPESPARK=1)
    
    def propagate(self, steps):
        for istep in steps:
            self.model.propagate(istep)
    
    def get_props(self, key):
        output = np.empty(self.model.NPART)

        for ip in range(self.model.NPART):
            particle   = self.model.particles[ip]
            output[ip] = particle[key]
        
        return output
    
    def update(self, key, values):
        for i, particle in enumerate(self.model.particles):
            particle[key] = values[i]


class EnKF:

    def __init__(self, grid, N=10, k_meas=25, sigma=0.1, NSTEPS=200):
        self.grid   = grid         # model grid 
        self.N      = N            # number of ensemble members
        self.k_meas = k_meas       # number of timesteps between measurements
        self.sigma  = sigma
        self.NSTEPS = NSTEPS       # number of time step for model
        # measurement covariance
        self.Cee = sigma**2 * np.eye(grid.get_npart())
        # measurement matrix
        self.M = np.eye(grid.get_npart())
        # initialize ensemble
        self.ensemble = Ensemble()
        for _ in range(N):
            self.ensemble.append(Member(grid, NSTEPS))
    

    def twin_experiment(self, data=None, key="state"):
        print('Starting the twin experiment')
        np.random.seed(2021)
        if data:
            self.experiment_data = data
            return
        # create model
        model = Lagrange(self.grid, self.NSTEPS)
        model.start_fire(ITYPESPARK=1)
        # twin experiment data
        self.experiment_data = np.empty((self.NSTEPS, self.grid.get_npart()), dtype=object)
        for ip in range(self.NSTEPS):
            model.propagate(ip)
            self.experiment_data[ip, :] = [p[key] for p in model.particles]


    def launch_filter(self, props="state"):
        t_it = 0
        for _ in range(self.NSTEPS//self.k_meas):
            # predict
            self.ensemble.propagate(range(t_it, t_it + self.k_meas + 1))
            t_it += self.k_meas
            print(t_it)

            # ensemble covariance a priori
            sample = [m_it.get_props(props) for m_it in self.ensemble]
            sample = np.array(sample).T
            C = np.cov(sample, ddof=1)

            # correction
            denom = self.Cee + self.M @ C @ self.M.T
            for m_it in self.ensemble:
                # measurement
                # states = [particle[props] for particle in self.experiment_data[t_it]]
                d = random.multivariate_normal(self.experiment_data, self.Cee)

                # u
                state = m_it.get_props(props)
                dstate = d - self.M @ state
                dstate = C @ self.M.T @ la.solve(denom, dstate)
                m_it.update(props, state + dstate)

            # ensemble covariance a posteriori
            sample = [m_it.get_props(props) for m_it in self.ensemble]
            sample = np.array(sample).T
            C_post = np.cov(sample, ddof=1)

            # ensemble covariance a posteriori theoretical
            dC = C @ self.M.T @ la.solve(denom, self.M @ C)
            C_th = C - dC