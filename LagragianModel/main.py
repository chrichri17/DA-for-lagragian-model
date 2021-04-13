# pylint: disable=import-error, unbalanced-tuple-unpacking
import numpy as np
from time import time
from copy import deepcopy
from Grid import Grid
from Lagrange import Lagrange
from EnKF import EnKF


# INITIALIZATION OF THE GRID
grid = Grid(NX=51, NY=51, XMAX=50*2, YMAX=50*2)

# SET NON FLAMMABLE CELLS
# using ratio = 0.5 (ratio is the fraction of non flammable cells);
# this also put few roads.
#
grid.set_non_flammable_cells(0.5)

# PLOT THE GRID
#
# grid.contour_plot()

# INITIALIZE THE LAGRAGIAN MODEL ON THE GRID
# we can specify parameter NSTEPS if we want but by defaut it
# is set to 400
#
model = Lagrange(grid, NSTEPS=200)

# START FIRE
# ITYPESPARK = 0 (point)  1 (line)
#
model.start_fire(ITYPESPARK=1)

# LAUNCH SIMULATION
#
# model.launch()

for istep in range(70):
    if istep%10 == 0:
        model.get_fire_line()
    model.propagate(istep)

# model.get_fire_line()
