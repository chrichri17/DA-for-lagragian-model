# pylint: disable=import-error, unbalanced-tuple-unpacking

import numpy as np
import numpy.random as random
import numpy.linalg as la
# from copy import deepcopy
from .ensemble import Ensemble
from .member import Member
from LagrangianModel.lagrange import Lagrange


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

class EnKF:
    """Implements the Ensemble Kalman Filtering for the lagrangian model.

    Init Parameters
    ---------------


    Attributes
    ----------

    
    Examples
    --------
    >>> from LagrangianModel.grid import Grid
    >>> from KalmanFiltering.EnKF import EnKF
    >>> grid = Grid(NX=51, NY=51, XMAX=250, YMAX=250)
    >>> grid.set_non_flammable_cells(0.25)
    >>> enkf = EnKF(grid, NSTEPS=100, N=2)
    >>> enkf.start_twin_experiment(save=True, file_name="data.out")
    >>> enkf.launch_filter()
    

    After saving the experiment_data, we can load it whenever we come
    and start the data assimilation
    >>> enkf.load_twin_experiment_data("data.out")
    >>> enkf.launch_filter()
    """
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