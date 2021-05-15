# pylint: disable=import-error, unbalanced-tuple-unpacking

import numpy as np
import numpy.random as random
import numpy.linalg as la
from copy import deepcopy
from lagrange import Lagrange


EXPERIMENT_START = """
============================
Starting the twin experiment
============================
Using :
grid = {grid},
LagrangeModel(grid, {nsteps})
"""

EXPERIMENT_END = """
Created twin experiment data with shape {shape}.
{save}
[END of twin experiement]
"""

DA_MESSAGE = """
==============================
Starting the data assimilation
==============================

Ensemble kalman filter parameters are:
--------------------------------------
Model grid: grid = {grid}
Number N of ensemble members: N = {n}
Number of timesteps between measurements: k_meas = {k_meas}
Measurement covariance std: sigma = {std}
"""

class Ensemble:
    """ Object to store all the member of our ensemble kalman filter

    Attributes
    ----------
    list (list) :  container storing all the member


    Methods
    -------
    append(member) : add a new member to the ensemble

    propagate(steps) : propagate the lagragian model for
                       each member of the ensemble
    """
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
    """ Implements a member for the ensemble Kalman filter (EnKF). 
    You are responsible for setting the various state variables 
    to reasonable values.


    Init Parameters
    ---------------
    grid (grid.Grid) : object of grid.Grid representing the area 
                       of the model.
    
    NSTEPS (int): number of simulation timestep used for the
                  lagragian model (lagrange.Lagrange)

    
    Attributes
    ----------
    model (lagrange.Lagrange) : lagrangian model using the given grid
                                For now, we only consider a fire line model
                                ie ITYPESPARK = 1
    
    states (np.ndarray): store all the fire_line for each step
                         object of np.ndarray
    
    current (int) : integer to store the index of the current
                    fire line
    


    Methods
    ----------
    propagate(steps) : propagate the lagrangian model
                       for each step in the iterable steps
    
    get_fire_line() : return the fire line using the current state
                      of the grid
    
    get_current_state() : return the current fire_line (the last one)
    
    update(values) : update the current fire_line
    """
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
    

    def start_twin_experiment(self, save=False, file_name=""):
        print(EXPERIMENT_START.format(grid=self.grid, nsteps=self.NSTEPS))
        np.random.seed(2021)
        # create model
        model = Lagrange(self.grid, self.NSTEPS)
        model.start_fire(ITYPESPARK=1)
        # twin experiment data
        shape = (self.NSTEPS, self.grid.NY - 1)
        self.experiment_data = np.empty(shape)
        for ip in range(self.NSTEPS):
            model.propagate(ip)
            # store data
            self.experiment_data[ip, :] = model.interpolate_fire_line()
        save_msg = ""
        if save:
            np.savetxt(file_name, self.experiment_data)
            save_msg = f"Saved data in {file_name}"
        print(EXPERIMENT_END.format(shape=shape, save=save_msg))
    
    def load_twin_experiment_data(self, file_name):
        self.experiment_data = np.loadtxt(file_name)

    def launch_filter(self):
        if not isinstance(self.experiment_data, np.ndarray):
            print("Can not start data assimilation. experiment_data are undefined")
            return
        print(DA_MESSAGE.format(grid=self.grid, n = self.N,
                                  k_meas=self.k_meas, std=self.sigma))

        t_it = 0
        for _ in range(self.NSTEPS//self.k_meas):
            # predict
            self.ensemble.propagate(range(t_it, t_it + self.k_meas))
            t_it += self.k_meas
            print(f"Measurement at t = {t_it}")
            t_it = min(t_it, self.NSTEPS - 1)

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