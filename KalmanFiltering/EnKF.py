import numpy as np
import os
import sys
from pathlib import Path

print('start')
# Avoid import error
file = Path(__file__).resolve()
parent, root = file.parent, file.parents[1]
print(parent, root)
sys.path.append(str(root))
# Additionally remove the current file's directory from sys.path
try:
    sys.path.remove(str(parent))
except ValueError: # Already removed
    pass

# from KalmanFiltering.Member import Member
# from KalmanFiltering.Ensemble import Ensemble
from algorithms.LagrangianModel.Lagrange import Lagrange

# class EnKF:

#     def __init__(self, grid, N=10, k_meas=25, sigma=0.1, NSTEPS=200):
#         self.grid   = grid         # model grid 
#         self.N      = N            # number of ensemble members
#         self.k_meas = k_meas       # number of timesteps between measurements
#         self.sigma  = sigma
#         self.NSTEPS = NSTEPS       # number of time step for model
#         # measurement covariance
#         self.Cee = sigma**2 * np.eye(grid.NY - 1)
#         # measurement matrix
#         self.M = np.eye(grid.NY - 1)
#         # initialize ensemble
#         self.ensemble = Ensemble()
#         for _ in range(N):
#             self.ensemble.append(Member(grid, NSTEPS))
    

#     def start_twin_experiment(self, key="state"):
#         print('============================')
#         print('Starting the twin experiment')
#         print('============================')
#         print(self.grid, f"NSTEPS = {self.NSTEPS}")
#         np.random.seed(2021)
#         # create model
#         model = Lagrange(self.grid, self.NSTEPS)
#         model.start_fire(ITYPESPARK=1)
#         # twin experiment data
#         shape = (self.NSTEPS, self.grid.NY - 1)
#         self.experiment_data = np.empty(shape)
#         for ip in range(self.NSTEPS):
#             # propagate model
#             model.propagate(ip)
#             # store data
#             self.experiment_data[ip, :] = model.interpolate_fire_line()
#         print(f"Creation experiment_data\nwith shape: {shape}")
#         print('-----------[END]------------')

#     def launch_filter(self):
#         if not isinstance(self.experiment_data, np.ndarray):
#             print("Can not start data assimilation. 'experiment_data' are undefined")
#             return
#         print('\n'*2)
#         print('==============================')
#         print('Starting the data assmiliation')
#         print('==============================')
#         print("Ensemble kalman filter paramaters are: ")
#         print("Model grid: ", self.grid)
#         print("Model NSTEPS: ", self.NSTEPS)
#         print("Number N of ensemble members: N = ", self.N)
#         print("Number of timesteps between measurements: k_meas = ", self.k_meas)
#         print("Measurement covariance std: sigma = ", self.sigma)

#         t_it = 0
#         for _ in range(self.NSTEPS//self.k_meas):
#             # predict
#             self.ensemble.propagate(range(t_it, t_it + self.k_meas + 1))
#             t_it += self.k_meas
#             print(f"Measurement at t = {t_it}")

#             # ensemble covariance a priori
#             sample = [m_it.get_current_state() for m_it in self.ensemble]
#             sample = np.array(sample).T
#             C = np.cov(sample, ddof=1)

#             # correction
#             denom = self.Cee + self.M @ C @ self.M.T
#             for m_it in self.ensemble:
#                 # measurement
#                 d = random.multivariate_normal(self.experiment_data[t_it], self.Cee)

#                 # u
#                 state = m_it.get_current_state()
#                 dstate = d - self.M @ state
#                 dstate = C @ self.M.T @ la.solve(denom, dstate)
#                 # UPDATE the state
#                 m_it.update(state + dstate)

#             # ensemble covariance a posteriori
#             sample = [m_it.get_current_state() for m_it in self.ensemble]
#             sample = np.array(sample).T
#             self.C_post = np.cov(sample, ddof=1)

#             # ensemble covariance a posteriori theoretical
#             dC = C @ self.M.T @ la.solve(denom, self.M @ C)
#             self.C_th = C - dC