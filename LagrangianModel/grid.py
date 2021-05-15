import numpy as np
from random import random
import matplotlib.pyplot as plt
from .cell import Cell

class Grid:
    """Object to model the domain grid.

    Params:
    -------

    IN REAL APPLICATION, this comes from outside UNITS: m
    NX     (int, optional): # of grid lines (NX-1) is the number of cells. Defaults to 101
    NY     (int, optional): # of grid columns (NY-1) is the number of cells. Defaults to 101
    XMIN (float, optional): lower bound of grid in x-dir. Defaults to 0
    XMAX (float, optional): upper bound of grid in x-dir. Defaults to 2000
    YMIN (float, optional): lower bound of grid in y-dir. Defaults to 0
    YMAX (float, optional): upper bound of grid in y-dir. Defaults to 2000
    """

    def __init__(self, NX=101, NY=101, XMIN=0, XMAX=2000, YMIN=0, YMAX=2000):
        self.NX = NX
        self.NY = NY
        self.XMIN = XMIN
        self.YMIN = YMIN
        self.XMAX = XMAX
        self.YMAX = YMAX
        self.DX = (XMAX - XMIN) / (NX-1)
        self.DY = (YMAX - YMIN) / (NY-1)
        # Cell centers
        self.XCELL = np.empty((NX-1, NY-1))
        self.YCELL = np.empty((NX-1, NY-1))
        # All Cells
        self.CELLS = np.empty((NX-1, NY-1), dtype=object)
        for i in range(NX-1):
            for j in range(NY-1):
                self.XCELL[i, j] = XMIN + self.DX*(i + 0.5)
                self.YCELL[i, j] = YMIN + self.DY*(j + 0.5)
                self.CELLS[i, j] = Cell()
    

    def set_non_flammable_cells(self, ratio, NROAD=5, IROAD=29):
        """Put a fraction `ratio` of non-flammable cells. This is to model the
        "vegetation coverage" or "concrete buildings" & "roads".
        """
        NNONFL = int(ratio*(self.NX - 1)*(self.NY - 1))
        for _ in range(NNONFL):
            IRAND = round(random()*(self.NX - 2))
            JRAND = round(random()*(self.NY - 2))
            self.CELLS[IRAND, JRAND].cannot_burn()
        
        # Put a few roads
        for i in range(IROAD, self.NX - 1):
            for j in range(NROAD):
                JROAD = round((j+1)*(self.NY - 1)/(NROAD + 1)) - 1
                self.CELLS[i, JROAD].cannot_burn()

        for j in range(self.NY - 1):
            self.CELLS[IROAD,   j].cannot_burn()
            self.CELLS[IROAD+1, j].cannot_burn()
        
        istart = int(0.4*(self.NX-1)) - 1
        iend   = int(0.5*(self.NX-1))
        jstart = int(0.45*(self.NY-1)) - 1
        jend   = int(0.55*(self.NY-1))
        for i in range(istart, iend):
            for j in range(jstart, jend):
                self.CELLS[i, j].cannot_burn()
    
    def get_cells_prop(self, prop_name):
        """Get an array of `prop_name` values for each cell of the grid."""
        DATA = np.empty((self.NX - 1, self.NY - 1))
        for i in range(self.NX - 1):
            for j in range(self.NY - 1):
                DATA[i, j] = self.CELLS[i, j][prop_name]
        return DATA

    def get_heat_released(self):
        return self.get_cells_prop("QMAXTR")

    def get_burn_state(self):
        return self.get_cells_prop("BURNSTAT")

    def get_npart(self):
        """Return the number of particle."""
        return 2*Cell.NPARTMAX*(self.NX - 1)*(self.NY - 1)

    def contour_plot(self):
        fig = plt.figure()
        ax = fig.gca()
        QMAXTR = self.get_heat_released()

        ax.contourf(self.XCELL, self.YCELL, QMAXTR)
        plt.show()
    
    def __repr__(self):
        return f"Grid(NX={self.NX}, NY={self.NY})"

