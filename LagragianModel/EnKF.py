# pylint: disable=import-error, unbalanced-tuple-unpacking

import numpy as np
import numpy.random as random
import numpy.linalg as la
from copy import deepcopy
from lagrange import Lagrange

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

    def __init__(self, grid, NSTEPS=200):
        self.model = Lagrange(deepcopy(grid), NSTEPS)
        self.model.start_fire(ITYPESPARK=1)
        self.states = np.empty((NSTEPS, grid.NY - 1))
        self.current = 0
    
    def propagate(self, steps):
        for istep in steps:
            self.model.propagate(istep)
            self.states[istep, :] = self.get_fire_line()
            self.current = istep
    
    def get_fire_line(self):
        return self.model.interpolate_fire_line()
    
    def update(self, values):
        self.states[self.current] = values
    
    def get_current_state(self):
        return self.states[self.current]


class EnKF:

    def __init__(self, grid, N=10, k_meas=25, sigma=0.1, NSTEPS=200):
        self.grid   = grid         # model grid 
        self.N      = N            # number of ensemble members
        self.k_meas = k_meas       # number of timesteps between measurements
        self.sigma  = sigma
        self.NSTEPS = NSTEPS       # number of time step for model
        # measurement covariance
        self.Cee = sigma**2 * np.eye(grid.NY - 1)
        # measurement matrix
        self.M = np.eye(grid.NY - 1)
        # initialize ensemble
        self.ensemble = Ensemble()
        for _ in range(N):
            self.ensemble.append(Member(grid, NSTEPS))
    

    def start_twin_experiment(self, key="state"):
        print('============================')
        print('Starting the twin experiment')
        print('============================')
        print(self.grid, f"NSTEPS = {self.NSTEPS}")
        np.random.seed(2021)
        # create model
        model = Lagrange(self.grid, self.NSTEPS)
        model.start_fire(ITYPESPARK=1)
        # twin experiment data
        shape = (self.NSTEPS, self.grid.NY - 1)
        self.experiment_data = np.empty(shape)
        for ip in range(self.NSTEPS):
            # propagate model
            model.propagate(ip)
            # store data
            self.experiment_data[ip, :] = model.interpolate_fire_line()
        print(f"Creation experiment_data\nwith shape: {shape}")
        print('-----------[END]------------')

    def launch_filter(self):
        if not isinstance(self.experiment_data, np.ndarray):
            print("Can not start data assimilation. experiment_data are undefined")
            return
        print('\n'*3)
        print('==============================')
        print('Starting the data assmiliation')
        print('==============================')
        print("Ensemble kalman filter paramaters are: ")
        print("Model grid: ", self.grid)
        print("Number N of ensemble members: N = ", self.N)
        print("Number of timesteps between measurements: k_meas = ", self.k_meas)
        print("Measurement covariance std: sigma = ", self.sigma)

        t_it = 0
        for _ in range(self.NSTEPS//self.k_meas):
            # predict
            self.ensemble.propagate(range(t_it, t_it + self.k_meas + 1))
            t_it += self.k_meas
            print(f"Measurement at t = {t_it}")

            # ensemble covariance a priori
            sample = [m_it.get_current_state() for m_it in self.ensemble]
            sample = np.array(sample).T
            C = np.cov(sample, ddof=1)

            # correction
            denom = self.Cee + self.M @ C @ self.M.T
            for m_it in self.ensemble:
                # measurement
                d = random.multivariate_normal(self.experiment_data[t_it], self.Cee)

                # u
                state = m_it.get_current_state()
                dstate = d - self.M @ state
                dstate = C @ self.M.T @ la.solve(denom, dstate)
                # UPDATE the state
                m_it.update(state + dstate)

            # ensemble covariance a posteriori
            sample = [m_it.get_current_state() for m_it in self.ensemble]
            sample = np.array(sample).T
            C_post = np.cov(sample, ddof=1)

            # ensemble covariance a posteriori theoretical
            dC = C @ self.M.T @ la.solve(denom, self.M @ C)
            C_th = C - dC