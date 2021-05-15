import numpy as np
from copy import deepcopy
from LagrangianModel.lagrange import Lagrange

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
