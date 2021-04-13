

############################ Twin experiment ##################################
# npart = model.NPART
# dt = 1
# tm = model.NSTEPS
# nt = tm // dt
# t = np.linspace(0, tm, nt+1)


# def h(u): # Observation operator
#     w = u
#     return w

# def Dh(u): # Jacobian of observation operator
#     n = len(u)
#     D = np.eye(n)
#     return D


# np.random.seed(seed=1)
# sig_m = 0.15                          # standard deviation for measurement noise
# R     = sig_m**2 * np.eye(npart)      # covariance matrix for measurement noise
# dt_m  = 25                            # time period between observations
# tm_m  = tm #// 4                       # maximum time for observations
# nt_m  = tm_m // dt_m                  # number of observation instants
# ind_m = (np.linspace(dt_m//dt, tm_m//dt, nt_m)).astype(int)
# t_m   = t[ind_m]                      # time integration
# uTrue = np.empty((npart, nt+1), dtype=object)

# start = time()
# uTrue[:, 0] = deepcopy(model.particles)
# km = 0
# w = np.empty((npart, nt_m), dtype=object)
# for k in range(nt):
#     model.propagate(k)
#     uTrue[:, k+1] = deepcopy(model.particles)
#     if (km < nt_m) and (k+1 == ind_m[km]):
#         w[:, km] = h(uTrue[:, k+1]) + np.random.normal(0, sig_m, [npart,])
#         km += 1

# end = time()
# print(f"Time = {end - start}")
# Time = 307.09354758262634