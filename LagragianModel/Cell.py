
# Cell constant parameters
IGNTIME    = 10.0
BURNDUR    = 120

class Cell:
        
    """Object to model a cell in our grid.

        Params:
        -------

        TREEAREA  (float): fraction of cell area occupied by trees
        BUILDAREA (float): fraction of cell area occupied by buildings
                           TREEAREA + BUILDAREA do not need to sum up to 1.
        CLOCK     (float): the time axis for the potential fire in each 
        NPARTTR     (int): total number of particles emitted from tree fire (0-10)
        NPARTBLD    (int): total number of particles emitted from building fire
        QMAXTR    (float): Maximum heat release (MW/m2) from trees in each cell
        QMAXBLD   (float): Maximum heat release (MW/m2) from buildings in each cell
        TIGNTR    (float): Ignition time delay for tree (since arrival of burning particle
                           in cell) for each cell (in s)
        TENDTR    (float): Fire duration for tree for each cell (s)
        TIGNBLD   (float): As for TIGNTR, but for the building
        TENDBLD   (float): As for TENDTR, but for the building
        BURNSTAT    (int): 0 in the beginning, will become 1 if visited by a burning particle.
        BURNPROG    (int): in the beginning, counts how many burning particles arrived
                           at this cell
        BURNSTAT2   (int): 1 if cell has all its particles, 0 when it has emptied them.

        Test numbers used for the above for now, more proper parameters and
        their dependence on wind & tree & building type to be used later.
        Other models for tree & building fire can be used, by adjusting NPART
        and their relation with QMAX and Q(t)
    """
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
        """ Define a cell that cannot burn."""
        self.QMAXTR  = 0.0
        self.QMAXBLD = 0.0

    def __getitem__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        raise KeyError(f"Cell doesn't have key '{key}'")
