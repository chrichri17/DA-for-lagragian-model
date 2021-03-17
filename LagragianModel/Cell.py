
"""
    Class cell that represent a cell in our grid.

    #### PARAMETERS FOR EACH CELL ("PLOT OF LAND") ####

    TREEAREA  = fraction of cell area occupied by trees
    BUILDAREA = fraction of cell area occupied by buildings
    TREEAREA + BUILDAREA do not need to sum up to 1.
    CLOCK    = the time axis for the potential fire in each 
    NPARTTR  = total number of particles emitted from tree fire (0-10)
    NPARTBLD = total number of particles emitted from building fire
    QMAXTR   = Maximum heat release (MW/m2) from trees in each cell
    QMAXBLD  = Maximum heat release (MW/m2) from buildings in each cell
    TIGNTR   = Ignition time delay for tree (since arrival of burning particle
               in cell) for each cell (in s)
    TENDTR   = Fire duration for tree for each cell (s)
    TIGNBLD  = As for TIGNTR,      but for the building
    TENDBLD  = As for TENDTR,      but for the building
    BURNSTAT = 0 in the beginning, will become 1 if visited by a burning particle.
    BURNPROG = 0 in the beginning, counts how many burning particles arrived
                 at this cell
    BURNSTAT2 = 1 if cell has all its particles, 0 when it has emptied them.

    Test numbers used for the above for now, more proper parameters and
    their dependence on wind & tree & building type to be used later.
    Other models for tree & building fire can be used, by adjusting NPART
    and their relation with QMAX and Q(t)
"""
# Cell constant parameters
IGNTIME    = 10.0
BURNDUR    = 120

class Cell:
    # Maximum number of particle
    NPARTMAX   = 1
    
    def __init__(self):
        # Cell's parameters
        self.TREEAREA  = 0.5
        self.BUILDAREA = 0.5
        self.CLOCK     = 0.0
        self.QMAXTR    = 10.0
        self.QMAXBLD   = 1.0
        self.TIGNTR    = IGNTIME
        self.TENDTR    = self.TIGNTR + BURNDUR
        self.TIGNBLD   = IGNTIME*5
        self.TENDBLD   = self.TIGNBLD + BURNDUR*5
        self.BURNSTAT  = 0.0
        self.BURNPROG  = 0.0
        self.BURNSTAT2 = 1.0
        self.NPARTTR   = Cell.NPARTMAX
        self.NPARTRAD  = Cell.NPARTMAX

    def cannot_burn(self):
        """ Define a cell that cannot burn """
        self.QMAXTR  = 0.
        self.QMAXBLD = 0.
    
    