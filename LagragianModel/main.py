from Grid import Grid
from Lagrange import Lagrange

# INITIALIZATION OF THE GRID
grid = Grid(NX=101, NY=101)

# SET NON FLAMMABLE CELLS
# using ratio = 0.5 (ratio is the fraction of non flammable cells);
# this also put few roads.

grid.set_non_flammable_cells(0.5)

# PLOT THE GRID

# grid.contour_plot()

# INITIALIZE THE LAGRAGIAN MODEL ON THE GRID
# we can specify parameter NSTEPS if we want but by defaut it
# is set to 400

model = Lagrange(grid)

# START FIRE
# ITYPESPARK = 0 (point)  1 (line)

model.start_fire(ITYPESPARK=1)

# model.launch()

n = 1

for istep in range(n):
    model.propagate(istep)

model.get_fire_line()

