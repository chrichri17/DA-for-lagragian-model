import numpy as np

def euler(rhs, state, dt, *args):
    k1 = rhs(state, *args)
    new_state = state + dt*k1
    return new_state


def RK4(rhs, state, dt, *args):
    k1 = rhs(state, *args)
    k2 = rhs(state+k1*dt/2, *args)
    k3 = rhs(state+k2*dt/2, *args)
    k4 = rhs(state+k3*dt, *args)
    new_state = state + (dt/6)*(k1+2*k2+2*k3+k4)
    return new_state

def Lorenz63(state, *args): #Lorenz 96 model
    sigma = args[0]
    beta = args[1]
    rho = args[2]
    x, y, z = state # Unpack the state vector
    f = np.zeros(3) # Derivatives
    f[0] = sigma*(y - x)
    f[1] = x*(rho - z) - y
    f[2] = x*y - beta*z
    return f


# parameters
sigma = 10.0
beta = 8.0/3.0
rho = 28.0
dt = 0.01
tm = 10
nt = int(tm/dt)
t = np.linspace(0, tm, nt+1)

############################ Twin experiment ##################################
def h(u): # Observation operator
    w = u
    return w

def Dh(u): # Jacobian of observation operator
    n =len(u)
    D = np.eye(n)
    return D

u0True = np.array([1, 1, 1]) # True initial conditions
np.random.seed(seed=1)
sig_m = 0.15  # standard deviation for measurement noise
R = sig_m**2 * np.eye(3) # covariance matrix for measurement noise
dt_m = 0.2 # time period between observations
tm_m = 2 # maximum time for observations
nt_m = int(tm_m/dt_m) # number of observation instants
ind_m = (np.linspace(int(dt_m/dt), int(tm_m/dt), nt_m)).astype(int)
t_m = t[ind_m] # time integration
uTrue = np.zeros([3, nt+1])
uTrue[:, 0] = u0True
km = 0
w = np.zeros([3, nt_m])
for k in range(nt):
    uTrue[:, k+1] = RK4(Lorenz63, uTrue[:, k], dt, sigma, beta, rho)
    if (km < nt_m) and (k+1 == ind_m[km]):
        w[:, km] = h(uTrue[:, k+1]) + np.random.normal(0, sig_m, [3,])
        km += 1

############################



def Jeuler(rhs, Jrhs, state, dt, *args):
    n =len(state)
    # k1 = rhs(state, *args)
    dk1 = Jrhs(state, *args)
    DM = np.eye(n) + dt * dk1
    return DM


def JRK4(rhs, Jrhs, state, dt, *args):
    n =len(state)
    k1 = rhs(state, *args)
    k2 = rhs(state + k1*dt/2, *args)
    k3 = rhs(state + k2*dt/2, *args)
    dk1 = Jrhs(state, *args)
    dk2 = Jrhs(state + k1*dt/2, *args) @ (np.eye(n) + dk1*dt/2)
    dk3 = Jrhs(state + k2*dt/2, *args) @ (np.eye(n) + dk2*dt/2)
    dk4 = Jrhs(state + k3*dt, *args) @ (np.eye(n) + dk3*dt)
    DM = np.eye(n) + (dt/6) * (dk1 + 2*dk2 + 2*dk3 + dk4)
    return DM


def JLorenz63(state, *args): # Jacobian of Lorenz 96 model
    sigma, beta, rho = args[:3]
    x, y, z = state # Unpack the state vector
    df = np.array([
        [ -sigma, sigma,     0],
        [rho - z,    -1,    -x],
        [      y,     x, -beta]
    ])
    return df

def KF(ub, w, H, R, B):
    """
        Implementation of the KF with linear dynamics and observation operator
    """
    # The analysis step for the Kalman filter in the linear case
    # i.e., linear model M and linear observation operator H

    n = ub.shape[0]
    # compute Kalman gain
    D = H @ B @ H.T + R
    K = B @ H @ np.linalg.inv(D)
    # compute analysis
    ua = ub + K @ (w - H@ub)
    P = (np.eye(n) - K@H) @ B
    return ua, P


def EKF(ub, w, ObsOp, JObsOp, R, B):
    # The analysis step for the extended Kalman filter with nonlinear dynamics
    # and nonlinear observation operator

    n = ub.shape[0]
    # compute Jacobian of observation operator at ub
    Dh = JObsOp(ub)
    # compute Kalman gain
    D = Dh @ B @ Dh.T + R
    K = B @ Dh.T @ np.linalg.inv(D)

    # compute analysis
    ua = ub + K @ (w-ObsOp(ub))
    P = (np.eye(n) - K@Dh) @ B
    return ua, P


#%% Application: Lorenz 63
########################### Data Assimilation #################################
u0b = np.array([2.0, 3.0, 4.0])
sig_b = 0.1
B = sig_b**2 * np.eye(3)
Q = 0.0*np.eye(3) # time integration
ub = np.zeros([3, nt+1])
ub[:, 0] = u0b
ua = np.zeros([3, nt+1])
ua[:, 0] = u0b
km = 0
for k in range(nt):
    # Forecast Step
    # background trajectory [without correction]
    ub[:, k+1] = RK4(Lorenz63, ub[:, k], dt, sigma, beta, rho)
    # EKF trajectory [with correction at observation times]
    ua[:, k+1] = RK4(Lorenz63, ua[:, k], dt, sigma, beta, rho)
    #compute model Jacobian at t_k
    DM = JRK4(Lorenz63, JLorenz63, ua[:, k], dt, sigma, beta, rho)
    #propagate the background covariance matrix
    B = DM @ B @ DM.T + Q
    if (km < nt_m) and (k+1 == ind_m[km]):
        # Analysis Step
        ua[:, k+1], B = EKF(ua[:, k+1], w[:, km], h, Dh, R, B)
        km += 1


import matplotlib.pyplot as plt

fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(10,8))
ax = ax.flat
for k in range(3):
    ax[k].plot(t, uTrue[k,:], label="True", linewidth = 3)
    ax[k].plot(t, ub[k,:], ':', label="Background", linewidth = 3)
    ax[k].plot(t[ind_m], w[k,:], 'o', fillstyle="none", label="Observation", markersize = 8, markeredgewidth = 2)
    ax[k].plot(t, ua[k,:], '--', label="Analysis", linewidth = 3)
    ax[k].set_xlabel("$t$",fontsize=22)
    ax[k].axvspan(0, tm_m, color="y", alpha=0.4, lw=0)

ax[0].legend(loc="center", bbox_to_anchor=(0.5, 1.25),ncol =4, fontsize=15)
ax[0].set_ylabel("$x(t)$")
ax[1].set_ylabel("$y(t)$")
ax[2].set_ylabel("$z(t)$")

# X, Y, Z = uTrue[0,:], uTrue[1,:], uTrue[2, :]
# X, Y = np.meshgrid(X, Y)

# plt.contour(X, Y, Z)
plt.show()