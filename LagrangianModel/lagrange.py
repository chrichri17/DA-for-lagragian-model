import numpy as np
import matplotlib.pyplot as plt
# from copy import deepcopy
from .particle import Particle
from .cell import Cell

# Message to show when the simulation is launched

MESSAGE = """
==========================================================================================
Lagrangian stochastic model for fire propagation in forests and the forest/urban interface
==========================================================================================

Model parameters : 
{grid}, NSTEPS = {nsteps}
"""

# INITIALISATIONS
#
# NOT ALL ARE USED IN CURRENT VERSION. MUST WORK FURTHER ON PROPAGATION MODEL
#
# TMEM = Particle memory (timescale of decay) in s. (Physically: time of
# flight beyond which gas parcel or firebrand not causing ignitions). Will
# determine the front thickness, to some extent.
#
# STTHR =  Threshold to decide if particle is "on". THIS COULD BE USED TO
# MODEL QUENCHING. OR, if made f(space) to model ignitability of the region.
#
# LANGFACTOR = Maximum factor (0-1) of the air Langevin walk step that the flame
# particles experience.
# RELT = randomize a bit the ignition time criterion so the propagation is
# not too jumpy
#
TMEM       = 10.0
STTHR      = 0.05
RELT       = 0.2
LANGFACTOR = 0.15
NSTEPS     = 400

# WIND SPEED AND LANGEVIN PARAMETERS
# All in SI units
# UX = mean wind speed in X-dir (later, it can be f(Z))
# UY = mean wind speed in Y-dir
# A = relative turbulence intensity (u'/<U>) - will be adjusted for proper
#     Langevin walk later
#
# CAREFUL: this may be a fraction of the real wind speed given empirical
# data of flame spread vs. wind speed. UPRIME & TTURB ARE TURBULENCE
# QUANTIIES. NEEDS WORK: I AM MIXING PURE AIR QUANTITIES WITH INFERRED
# MEAN FIRE PROPAGATION SPEED (FIRECONST * AIR SPEED), MUST FIX LATER.
#
# Just putting this down from Wikipedia. Not using it yet.
# BEAUFORT=5
# WIND SPEED = 0.836*(BEAUFORT)^(3/2)  m/s
#
# ULAM : a constant to take care of finite propagation with no wind
# TMEMRAD : the "memory" of the radiation particles calculated as
# Length/ULAM, Length = 10m or something else
#
ULAM      = 0.1
LRAD      = 10.0
TMEMRAD   = LRAD / ULAM
UX        = 25.0
VY        = 0.0
FIRECONST = 1.0
UXM       = UX*FIRECONST
VYM       = VY*FIRECONST
USPEED    = np.sqrt(UXM**2 + VYM**2)
A         = 0.3
UPRIME    = USPEED*A
TLENGTH   = 30.0
TFREQ     = UPRIME / TLENGTH
CLANG     = np.sqrt(2.0*UPRIME**3 / TLENGTH)


def normal(mu=0, sigma=1):
    return np.random.normal(mu, sigma)

class Lagrange:
    """Implement the lagragian model for fire propagation in forests and the forest/urban interface

    Params:
    -------
    grid       (grid.Grid): grid representing the surface area
    NSTEPS (int, optional): number of time to propagate the model to/
                            number of timesteps of the simulation
    

    Attributes:
    -----------
    DT       (float): timestep in (s)
    SQRDT    (float): square root of DT
    NPART      (int): number of particles in the grid 
    ITYPESPARK (int): number to specify the type of the initial spark
                      1 for fire line, 0 for a point

    TIME  (np.ndarray[float]): the time axis for the whole simulation
    BURNX (np.ndarray[float]): burning progress averaged in y-dir as a function of time & x-dir


    Methods:
    --------
    start_fire(ITYPESPARK): start the fire by setting either a front line
                            or a point fire.
    propagate(istep): propagate the model for the given timestep.


    """
    def __init__(self, grid, NSTEPS=NSTEPS):
        # self.grid_copy = deepcopy(grid)
        self.grid       = grid
        self.NSTEPS     = NSTEPS
        self.ITYPESPARK = None
        if USPEED > 0:
            self.DT = 2.0 * (self.grid.DX / USPEED)
        else:
            self.DT = 1.0 * (self.grid.DX / ULAM)
        self.SQRTDT = np.sqrt(self.DT)
        self.TIME = np.empty(self.NSTEPS)
        for i in range(self.NSTEPS):
            self.TIME[i] = (i + 1) * self.DT
        
        self.BURNX = np.empty((NSTEPS, self.grid.NX-1))
        # SAVE BURNSTAT AT EVERY TIMESTEP
        self.BURNEVOL = np.empty((NSTEPS, self.grid.NX-1, self.grid.NY-1))

        self.NPART = 0
        self.particles = None # particles array
        self._init_particles()
    
    def _init_particles(self):
        """Initialisation of particle's global state.
            
            Use large memory for the time being (all particles stay lit for a long
            time)
            Give all particles some random velocity component, so the drift in the
            later Langevin walk makes sense.
        """
        self.NPART = self.grid.get_npart()
        self.particles = np.empty(self.NPART, dtype=object)
        for i in range(self.NPART):
            tmem = TMEM
            ux = UXM + UPRIME*normal()*LANGFACTOR
            vy = VYM + UPRIME*normal()*LANGFACTOR
            self.particles[i] = Particle(tmem=tmem, ux=ux, vy=vy)
        #
        # PUT THE PARTICLES IN THE CELLS.
        # LOOP OVER CELLS AND DEFINE THEIR PARTICLES.
        # FOR NOW, ONLY POSITION DEPENDS ON SPACE HEIGHT & MEMORY DO NOT.
        # FIRST THE TREE PARTICLES, THEN THE BUILDING PARTICLES.
        #
        NX = self.grid.NX
        NY = self.grid.NY
        icounter = 0
        for i in range(NX - 1):
            for j in range(NY - 1):
                cell = self.grid.CELLS[i, j]
                x    = self.grid.XCELL[i, j]
                y    = self.grid.YCELL[i, j]
                for k in range(cell.NPARTTR):
                    self.particles[k + icounter].update(x=x, y=y, type=1)
                for k in range(cell.NPARTRAD):
                    self.particles[k + cell.NPARTTR + icounter].update(x=x, y=y, type=2)
                icounter += cell.NPARTTR + cell.NPARTRAD
    
    # def clear(self):
    #     self.__init__(self.grid_copy, self.NSTEPS)

    def start_fire(self, ITYPESPARK=1):
        """Start the fire
            
            ITYPESPARK = 0 (point)  1 (line)
        """
        self.ITYPESPARK = ITYPESPARK
        if ITYPESPARK == 0:
            IF = 4
            JRANGE = [49]
        elif ITYPESPARK == 1:
            IF = 1
            JRANGE = range(self.grid.NY - 1)
        for JF in JRANGE:
            self.grid.CELLS[IF, JF].BURNSTAT = 1
            for ip in range(self.NPART):
                particle = self.particles[ip]
                X = self.grid.XCELL[IF, JF]
                Y = self.grid.YCELL[IF, JF]
                x, y = particle.x, particle.y
                if (abs(X - x) < self.grid.DX) and (abs(Y - y) < self.grid.DY):
                    particle.update(state=1, factor=LANGFACTOR)
    
    def _move_particle(self, ip):
        """
            RANDOM WALK & DECAY STATUS TO ACCOUNT FOR MEMORY.
            IF PARTICLES LEAVE DOMAIN, IGNORE THEM.
            MOVE ONLY PARTICLES THAT ARE "ON".
            RANDOM WALK FOLLOWS CNF PAPER.
            
            NEW FEATURE: Multiply step by a FACTOR, which reflects the "age" of the
            particle relative to the clock of the fire in each cell. Hence, 0 causes
            no propagation, 1 full movement by the wind turbulence (note: in this
            implementation, FIRECONST = 1). 
            
            NEW FEATURE: THE "BUILDING" PARTICLES MODEL RADIATION, HENCE MOVE
            DIFFERENTLY
        """
        # RANDOM WALK
        # ADVANCE ONLY THE PARTICLES THAT ARE "ON" (i.e. ABOVE STTHR).
        #
        particle = self.particles[ip] # get particle
        props = ["state", "type", "x", "y", "ux", "vy", "factor", "tmem"]
        state, pType, x, y, ux, vy, factor, tmem = particle.get_from_keys(props)
        if state > STTHR and pType == 1:
            DU  = -(ux - UXM)*2.0*TFREQ*self.DT + CLANG*self.SQRTDT*normal()
            DV  = -(vy - VYM)*2.0*TFREQ*self.DT + CLANG*self.SQRTDT*normal()
            UXP = ux + DU
            VYP = vy + DV
            XP  = x + UXP*self.DT*factor
            YP  = y + VYP*self.DT*factor
            STP = state*np.exp(-self.DT/tmem)
            particle.update(ux=UXP, vy=VYP, x=XP, y=YP, state=STP)
        elif (state > STTHR) and pType == 2:
            DU = ULAM*normal()
            DV = ULAM*normal()
            XP = x + DU*self.DT
            YP = y + DV*self.DT
            STP = state*np.exp(-self.DT/ TMEMRAD)
            particle.update(x=XP, y=YP, state=STP)
        if x > self.grid.XMAX - self.grid.DX:
            particle.update(x=self.grid.XMAX - self.grid.DX, state=0.)
        elif x < self.grid.XMIN + self.grid.DX:
            particle.update(x=self.grid.XMIN + self.grid.DX, state=0.)
        if y > self.grid.YMAX - self.grid.DY:
            particle.update(y=self.grid.YMAX - self.grid.DY, state=0.)
        elif y < self.grid.YMIN + self.grid.DY:
            particle.update(y=self.grid.YMIN + self.grid.DY, state=0.)
    
    def _ignite_cells(self, istep, ip):
        """ 
            IGNITE CELLS VISITED BY BURNING PARTICLES, BUT ONLY IF VISITED
            FOR THE FIRST TIME.  
            NOTE: WORKS ONLY FOR CONSTANT NUMBER OF PARTICLES (NPARTMAX), OR MUST
            RE-CALCULATE THE INDEX indp DIFFERENTLY.
            NOTE: STUPID WAY TO GET INDY & INDX, BUT USE OF "ROUND" DID NOT WORK. WORKS
            FOR STRUCTURED GRIDS ONLY.
            
            To see random walk alone, put state > STTHR*10 below
            to avoid doing any ignitions.
            
            8/3/20: In this version, release all particles from newly-ignited cell
            after a delay TIGNTR
            
            The commented out "elif" can be used to completely stop particles, e.g. from a tall building.
            Land without flammables is modelled by not releasing new particles.
        """
        particle = self.particles[ip] # get particle
        state, x, y = particle.get_from_keys(["state", "x", "y"])
        if state > STTHR:
            for i in range(self.grid.NX-1):
                if abs(x - self.grid.XCELL[i, 0]) < self.grid.DX/2:
                    INDX = i
            for j in range(self.grid.NY-1):
                if abs(y - self.grid.YCELL[0, j]) < self.grid.DY/2:
                    INDY = j
            cell = self.grid.CELLS[INDX, INDY]
            cell.BURNPROG += 1
            if (cell.QMAXTR > 0 or cell.QMAXBLD > 0) and cell.BURNSTAT == 0:
                cell.BURNSTAT = 1
                cell.CLOCK    = self.TIME[istep]
            # elif cell.QMAXTR == 0 or cell.QMAXBLD == 0:
            #     particle.update(state=0.0, factor=0.0)
            # if pType == 2:
            #     particle.update(state=0.0)
    
    def _launch_particles(self, istep):
        """
            Launch new particles from the ignited cells. Randomise a bit the
            ignition delay time.
        """
        for i in range(self.grid.NX-1):
            for j in range(self.grid.NY-1):
                INDX = i
                INDY = j
                cell = self.grid.CELLS[INDX, INDY]
                TLOCAL = self.TIME[istep] - cell.CLOCK
                TCRIT  = cell.TIGNTR * (1 + RELT*normal())
                if cell.BURNSTAT == 1 and TLOCAL > TCRIT and cell.BURNSTAT2 == 1:
                    LOCALF = LANGFACTOR
                    indp = (INDX*(self.grid.NY - 1) + INDY)*2*Cell.NPARTMAX - 1
                    for k in range(cell.NPARTTR):
                        self.particles[k + indp].update(state=1.0, factor=LOCALF)
                    for k in range(cell.NPARTRAD):
                        self.particles[k + cell.NPARTTR + indp].update(state=1.0, factor=LOCALF)
                    cell.BURNSTAT2 = 0
    
    def save_burning_state(self, istep):
        if self.ITYPESPARK == 1:
            for i in range(self.grid.NX-1):
                for j in range(self.grid.NY-1):
                    self.BURNX[istep, i] += self.grid.CELLS[i, j].BURNPROG
                self.BURNX[istep, i] /= self.grid.NY - 1
        self.BURNEVOL[istep, :, :] = self.grid.get_burn_state()[:, :]

    def propagate(self, istep):
        for ip in range(self.NPART):
            self._move_particle(ip)
            self._ignite_cells(istep, ip)
        self._launch_particles(istep)
    
    def get_particles_props(self, *props, array=None):
        """Get all particles props given the props name.

        Args:
            array (np.ndarray or None, optional): the array of particles from
                                                  which we want to fetch the props. 
                                                  Defaults to None.

        Returns:
            properties (np.ndarray): the properties of the particles in array
        """
        n = len(props)
        if not isinstance(array, np.ndarray):
            array = self.particles
        if n == 0:
            return []
        elif n == 1:
            key = props[0]
            properties = np.empty(len(array))
        else:
            properties = np.empty((n, len(array)))
        
        for ip, particle in enumerate(array):
            if n == 1:
                properties[ip] = particle[key]
            else:
                for i, key in enumerate(props):
                    properties[i, ip] = particle[key]
        return properties



    def scatter_plot(self, speed=0.001):
        """Plot the current grid state."""
        FPARTX, FPARTY, FPARTST = self.get_particles_props('x', 'y', 'state')
        plt.clf()
        plt.scatter(FPARTX, FPARTY, FPARTST*5 + 0.01, c=FPARTST, cmap="jet")
        plt.clim(0, 1)
        plt.colorbar()
        plt.pause(speed)

    def launch(self):
        """Launch the simulation."""
        print(MESSAGE.format(grid=self.grid, nsteps=self.NSTEPS))

        iplotcount = 0
        for istep in range(self.NSTEPS):
            self.propagate(istep)
            self.scatter_plot()
            
            iplotcount += 1
            if iplotcount == 10:
                iplotcount = 0
            
            # self.save_burning_state(istep)
        print("[END of simulation]")
        plt.show()
    
    def get_fire_line(self, epsilon=0.01, line_only=False):
        """Extract the fire line."""
        # Get active particles
        ACTIVE_PARTICLES = self.particles[self.particles > STTHR]
        # get maximum burning state
        MAX = max(self.particles, key=lambda p: p.state)
        Max = MAX.state

        # fetch  fire  line
        FIRE_LINE = ACTIVE_PARTICLES[Max - epsilon <= ACTIVE_PARTICLES]
        return FIRE_LINE if line_only else ACTIVE_PARTICLES, FIRE_LINE, MAX
    
    def interpolate_fire_line(self):
        """Interpolate fire line to get a fire particle state at each cell center."""
        fire_line = self.get_fire_line(line_only=True)
        LINEX, LINEY = self.get_particles_props('x', 'y', array=fire_line)
        coef = np.polyfit(LINEY, LINEX, deg=min(10, len(LINEX)))
        fit  = np.poly1d(coef)
        return fit(self.grid.YCELL[0])


    def plot_fire_line(self, epsilon=0.01, speed=10, stop=200):
        """ Plot the fire line at specific time :
            Launch a simulation and plot the fire line at each step.
        """
        stop = min(stop, self.NPART)
        # PLOTS setup
        _, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 5))
        
        for istep in range(stop):
            self.propagate(istep)
            if istep % speed == 0:
                print(f"{istep = }")
                active, line, _ = self.get_fire_line(epsilon)
                XP, YP, ST = self.get_particles_props('x', 'y', 'state', array=active)
                LINEX, LINEY, LINEST = self.get_particles_props('x', 'y', 'state', array=line)

                coef = np.polyfit(LINEY, LINEX, deg=min(10, len(LINEX)))
                fit = np.poly1d(coef)
                data = fit(self.grid.YCELL[0])

                STinterp = np.interp(self.grid.YCELL[0], LINEY, LINEST)
                plt.cla()
                z1 = ax1.scatter(XP, YP, c=ST, cmap="jet")
                z2 = ax2.scatter(LINEX, LINEY, c=LINEST, cmap="jet")
                z3 = ax3.scatter(data, self.grid.YCELL[0], c=STinterp, cmap="jet")

                if istep == 0: plt.colorbar(z1)
                for ax, z in zip((ax1, ax2, ax3), (z1, z2, z3)):
                    ax.set_xlim(0, self.grid.XMAX)
                    z.set_clim(0, 1)
                plt.pause(1)
            if istep == stop - 1:
                print("[END of SIMULATION]")
                plt.pause(3)
                plt.close()
        plt.show()

        