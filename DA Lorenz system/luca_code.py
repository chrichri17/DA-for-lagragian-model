import numpy as np
import numpy.linalg as la
import numpy.random as random
import matplotlib.pyplot as plt

# set seed
random.seed(0)

# lorenz parameters
rho   = 28.
sigma = 10.
beta  = 8./3.

# Time integration method
def odeint(func, state, t):
    res = np.empty((t.size, len(state)))
    res[0, :] = state[:]

    for i in range(1, t.size):
        dydt = func(res[i-1, :], t[i-1])
        res[i, :] = res[i-1, :] + dydt * (t[i] - t[i-1])
    
    return res

# Lorenz model
def Lorenz63(state, t):
    x, y, z = state # Unpack the state vector
    
    # return the derivatives
    return np.array([
        sigma*(y - x),
        x*(rho - z) - y,
        x*y - beta*z
    ])

### Resolution
t = np.linspace(0., 40., 16001)
state_0 = np.array([1., 1., 1.])
states = odeint(Lorenz63, state_0, t)

print(states[-1, :])


### Ensemble implementation

class Ensemble:

    def __init__(self):
        self.list = []
    
    def __iter__(self):
        return iter(self.list)
    
    def __getitem__(self, key):
        return self.list[key]
    
    def append(self, member):
        self.list.append(member)
    
    def odeint(self, func, t):
        for member in self.list:
            member.odeint(func, t)

class Member:

    def __init__(self, state_0):
        self.states = [state_0]
    
    def odeint(self, func, t):
        self.states = odeint(func, self.states[0], t)



N = 2 # number of nsemble members
cov_0 = 0.1**2 * np.eye(3)

ensemble = Ensemble()
for i in range(N):
    state = random.multivariate_normal(state_0, cov_0)
    ensemble.append(Member(state))

# ensemble.odeint(Lorenz63, t)

# for member in ensemble:
#     print(member.states[-1])

# Error analysis
def plot_error(state_1, state_2):
    error = np.empty(t.size)

    for i in range(t.size):
        dx = state_2[i, 0] - state_1[i, 0]
        dy = state_2[i, 1] - state_1[i, 1]
        dz = state_2[i, 2] - state_1[i, 2]

        error[i] = np.sqrt(dx**2 + dy**2 + dz**2)
    
    plt.plot(t, error)
    plt.show()

# state_1 = ensemble[0].states
# state_2 = ensemble[1].states
# plot_error(state_1, state_2)


### Ensemble Kalman Filter
N      = 100                    # number of ensemble members
k_meas = 1000                   # number of timesteps between measurements
Cee    = 0.1**2 * np.eye(3)     # measurement covariance
M      = np.eye(3)              # measurement matrix

ensemble = Ensemble()
for i in range(N):
    state = random.multivariate_normal(state_0, cov_0)
    ensemble.append(Member(state))

t_it = 0
for meas_it in range(t.size//k_meas):

    # predict
    ensemble.odeint(Lorenz63, t[t_it:t_it+k_meas+1])
    t_it += k_meas
    print(t_it)

    # ensemble covariance a priori
    sample = [m_it.states[-1] for m_it in ensemble]
    sample = np.array(sample).T
    C = np.cov(sample, ddof=1)

    # correction
    denom = Cee + M @ C @ M.T
    for m_it in ensemble:
        # measurement
        d = random.multivariate_normal(states[t_it], Cee)

        # u
        state = m_it.states[-1]
        dstate = d - M @ state
        dstate = C @ M.T @ la.solve(denom, dstate)
        m_it.states[-1] = state + dstate

    # ensemble covariance a posteriori
    sample = [m_it.states[-1] for m_it in ensemble]
    sample = np.array(sample).T
    C_post = np.cov(sample, ddof=1)

    # ensemble covariance a posteriori theoretical
    dC = C @ M.T @ la.solve(denom, M @ C)
    C_th = C - dC
    print(C_th - C_post)
