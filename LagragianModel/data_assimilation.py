# pylint: disable=import-error, unbalanced-tuple-unpacking
import numpy as np
from grid import Grid
from EnKF import EnKF

# INITIALIZATION OF THE GRID
grid = Grid(NX=51, NY=51, XMAX=250, YMAX=250)

# SET NON FLAMMABLE CELLS
# using ratio = 0.5 (ratio is the fraction of non flammable cells);
# this also put few roads.
#
grid.set_non_flammable_cells(0.25)


### Data assimilation with twin experiment
enkf = EnKF(grid, NSTEPS=100, N=2)
# enkf.start_twin_experiment()
enkf.load_twin_experiment_data("data.out")
enkf.launch_filter()
