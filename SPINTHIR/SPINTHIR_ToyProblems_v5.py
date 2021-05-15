import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import random

# get random value from normal distribution mu = 0 sigma = 1

def normal(mu=0, sigma=1):
    return np.random.normal(mu, sigma)

def uniform():
    return random.uniform(0, 1)

# DEFINE DOMAIN AND GRID
# IN REAL APPLICATION, this comes from outside
# UNITS: m
# NX = # of grid lines (NX-1) is the number of cells
# NY = same in y-dir

NX = NY = 101
XMIN, XMAX = 0, 2000
YMIN, YMAX = 0, 2000
DX = (XMAX - XMIN) / (NX-1)
DY = (YMAX - YMIN) / (NY-1)
GRIDX = np.empty(NX)
GRIDY = np.empty(NY)
for i in range(NX):
    GRIDX[i] = XMIN + DX*i
for i in range(NY):
    GRIDY[i] = YMIN + DY*i

# Cell centers
XCELL = np.zeros((NX-1, NY-1))
YCELL = np.zeros((NX-1, NY-1))
for i in range(NX-1):
    for j in range(NY-1):
        XCELL[i,j] = (GRIDX[i+1] + GRIDX[i]) / 2
        YCELL[i,j] = (GRIDY[j+1] + GRIDY[j]) / 2

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
ULAM    = 0.1
LRAD    = 10.0
TMEMRAD = LRAD / ULAM
#
#UX       = 25.0
#UX       = 7.77
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
CLANG     = np.sqrt(2*UPRIME**3 / TLENGTH)
#
# DEFINE PARAMETERS FOR EACH CELL ("PLOT OF LAND")
#
# TREEAREA = fraction of cell area occupied by trees
# BUILDAREA = fraction of cell area occupied by buildings
# TREEAREA + BUILDAREA do not need to sum up to 1.
# CLOCK = the time axis for the potential fire in each cell.
# NPARTTR = total number of particles emitted from tree fire (0-10)
# NPARTBLD = total number of particles emitted from building fire
# QMAXTR = Maximum heat release (MW/m2) from trees in each cell
# QMAXBLD = Maximum heat release (MW/m2) from buildings in each cell
# TIGNTR = Ignition time delay for tree (since arrival of burning particle
#          in cell) for each cell (in s)
# TENDTR = Fire duration for tree for each cell (s)
# TIGNBLD = As for TIGNTR, but for the building
# TENDBLD = As for TENDTR, but for the building
# BURNSTAT = 0 in the beginning, will become 1 if visited by a burning particle.
# BURNPROG = 0 in the beginning, counts how many burning particles arrived
#            at this cell
# BURNSTAT2 = 1 if cell has all its particles, 0 when it has emptied them.
#
# Test numbers used for the above for now, more proper parameters and
# their dependence on wind & tree & building type to be used later.
# Other models for tree & building fire can be used, by adjusting NPART
# and their relation with QMAX and Q(t)
#
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
NPARTMAX   = 1
TMEM       = 10.0
STTHR      = 0.05
IGNTIME    = 10.0
RELT       = 0.2
BURNDUR    = 120
LANGFACTOR = 0.15
NSTEPS     = 400


TREEAREA  = np.empty((NX-1,NY-1))
BUILDAREA = np.empty((NX-1,NY-1))
CLOCK     = np.empty((NX-1,NY-1))
QMAXTR    = np.empty((NX-1,NY-1))
QMAXBLD   = np.empty((NX-1,NY-1))
TIGNTR    = np.empty((NX-1,NY-1))
TENDTR    = np.empty((NX-1,NY-1))
TIGNBLD   = np.empty((NX-1,NY-1))
TENDBLD   = np.empty((NX-1,NY-1))
BURNSTAT  = np.empty((NX-1,NY-1))
NPARTTR   = np.empty((NX-1,NY-1), dtype=int)
NPARTRAD  = np.empty((NX-1,NY-1), dtype=int)
BURNPROG  = np.empty((NX-1,NY-1))
BURNSTAT2 = np.empty((NX-1,NY-1))
for i in range(NX-1):
    for j in range(NY-1):
        TREEAREA [i,j] = 0.5
        BUILDAREA[i,j] = 0.5
        CLOCK    [i,j] = 0.0
        QMAXTR   [i,j] = 10.0
        QMAXBLD  [i,j] = 1.0
        TIGNTR   [i,j] = IGNTIME
        TENDTR   [i,j] = TIGNTR[i,j] + BURNDUR
        TIGNBLD  [i,j] = IGNTIME*5
        TENDBLD  [i,j] = TIGNBLD[i,j] + BURNDUR*5
        BURNSTAT [i,j] = 0.0
        BURNPROG [i,j] = 0.0
        BURNSTAT2[i,j] = 1.0
        NPARTTR  [i,j] = NPARTMAX
        NPARTRAD [i,j] = NPARTMAX


# Put a fraction FNON of non-flammable cells this is to model the
# "vegetation coverage" or concrete buildings & roads.
#
IFLAG = 0
if IFLAG == 1:
    NNONFL = int(0.09*(NX-1)*(NY-1))
    for _ in range(NNONFL):
        IRAND = round(uniform()*(NX - 2))
        JRAND = round(uniform()*(NY - 2))
        QMAXTR [IRAND, JRAND] = 0.0
        QMAXBLD[IRAND, JRAND] = 0.0


IFLAG = 1
if IFLAG == 1:
    NNONFL = int(0.5*(NX-1)*(NY-1))
    for _ in range(NNONFL):
        IRAND = round(uniform()*(NX - 2))
        JRAND = round(uniform()*(NY - 2))
        QMAXTR [IRAND, JRAND] = 0.0
        QMAXBLD[IRAND, JRAND] = 0.0

    # Put a few roads

    NROAD = 5
    IROAD = 29
    for i in range(IROAD, NX-1):
        for j in range(NROAD):
            JROAD = round((j+1)*(NY - 1)/(NROAD + 1)) - 1
            QMAXTR [i, JROAD] = 0.0
            QMAXBLD[i, JROAD] = 0.0
            # QMAXTR [i, JROAD+1] = 0.0
            # QMAXBLD[i, JROAD+1] = 0.0
    for j in range(NY-1):
        QMAXTR [IROAD, j] = 0.0
        QMAXBLD[IROAD, j] = 0.0
        QMAXTR [IROAD+1, j] = 0.0
        QMAXBLD[IROAD+1, j] = 0.0

    for i in range(39, 50):
        for j in range(44, 55):
            QMAXTR [i,j] = 0.0
            QMAXBLD[i,j] = 0.0



fig = plt.figure()
ax = fig.gca()

#surf = ax.plot_surface(XCELL, YCELL, QMAXTR, cmap=cm.winter)
ax.contourf(XCELL, YCELL, QMAXTR)
#ax.view_init(azim=0, elev=90)

#plt.contour(XCELL, YCELL, QMAXTR, 10)

#Show the plot
plt.show()

#
# PARAMETERS FOR FORWARD INTEGRATION
# UNITS: S
# DT = timestep (in s)
# NTSTEPS = # of timesteps of the simulation
# TIME = the time axis for the whole simulation
#
# NOTE: Proper diffusion process needs care on DT it depends on # of
# new particles emitted AND on DX/U. Must experiment. See Appendix of 2012
# CNF paper. Here, I'm playing with the DT as a function of a mean cell transit time
#
# DT=2.0

if USPEED > 0:
    DT = 2.0 * (DX / USPEED)
else:
    DT = 1.0 * (DX / ULAM)

SQRTDT = np.sqrt(DT)
TIME = np.zeros(NSTEPS)
for i in range(NSTEPS):
    TIME[i] = (i + 1)*DT

#
# INITIALISE PARTICLE-RELATED ARRAYS
#
# FPARTST = For each particle, a value 0 or 1 ("cold" or "ignited"). 
#            If "ignited", it moves until it dies according to a decay.
# FPARTMEM = Timescale of the decay (s) for each particle -  for now, this
#            reflects the QMAX
# FPARTH = Height of release (m) for each particle - for a 3-D implementation
# FPARTX = x-position of each particle
# FPARTY = y-position of each particle
#
# FACTOR = From 0 to 1, according to a Gaussian, centred at max burning
# time, with spread so as to account for IGNTIME & BURNDUR.
# In this version, all cells emit equal number of particles, but they
# could have different MEMORY to account for different QMAX
#
# TOTAL NUMBER OF PARTICLES
#

NPART = np.sum(NPARTTR) +  np.sum(NPARTRAD)

# 
# INITIALISATION
# Use large memory for the time being (all particles stay lit for a long
# time)
# Give all particles some random velocity component, so the drift in the
# later Langevin walk makes sense.
#
FPARTST   = np.zeros(NPART)
FPARTMEM  = np.zeros(NPART)
FPARTH    = np.zeros(NPART)
FPARTX    = np.zeros(NPART)
FPARTY    = np.zeros(NPART)
FPARTUX   = np.zeros(NPART)
FPARTVY   = np.zeros(NPART)
FACTOR    = np.zeros(NPART)
FPARTTYPE = np.zeros(NPART)

for i in range(NPART):
    FPARTMEM[i] = TMEM
    FPARTUX [i] = UXM + UPRIME*normal()*LANGFACTOR
    FPARTVY [i] = VYM + UPRIME*normal()*LANGFACTOR
    FACTOR  [i] = 0.0

#
# PUT THE PARTICLES IN THE CELLS.
# LOOP OVER CELLS AND DEFINE THEIR PARTICLES.
# FOR NOW, ONLY POSITION DEPENDS ON SPACE HEIGHT & MEMORY DO NOT.
# FIRST THE TREE PARTICLES, THEN THE BUILDING PARTICLES.
#
icounter = 0
for i in range(NX-1):
    for j in range(NY-1):
        for k in range(NPARTTR[i, j]):
            FPARTX   [k + icounter] = XCELL[i, j]
            FPARTY   [k + icounter] = YCELL[i, j]
            FPARTTYPE[k + icounter] = 1
        for k in range(NPARTRAD[i, j]):
            FPARTX   [k + NPARTTR[i, j] + icounter] = XCELL[i, j]
            FPARTY   [k + NPARTTR[i, j] + icounter] = YCELL[i, j]
            FPARTTYPE[k + NPARTTR[i, j] + icounter] = 2
        icounter += NPARTTR[i, j] + NPARTRAD[i, j]

#
# START FIRE
#
# ITYPESPARK = 0 (point)  1 (line)
#
ITYPESPARK = 1
#
if ITYPESPARK == 0:
    # IF=(NX-1)/2
    IF = 4  
    # JF=(NY-1)/2
    JF = 49
    BURNSTAT[IF, JF] = 1
    for ip in range(NPART):
        if (abs(FPARTX[ip] - XCELL[IF, JF]) < DX) and (abs(FPARTY[ip] - YCELL[IF, JF]) < DY):
            FPARTST[ip] = 1.0
            FACTOR [ip] = LANGFACTOR
elif ITYPESPARK == 1:
    IF = 1
    for JF in range(NY-1):
        BURNSTAT[IF, JF] = 1
        for ip in range(NPART):
            if (abs(FPARTX[ip] - XCELL[IF, JF]) < DX) and (abs(FPARTY[ip] - YCELL[IF, JF]) < DY):
                FPARTST[ip] = 1.0
                FACTOR [ip] = LANGFACTOR

#
# RANDOM WALK
# ADVANCE ONLY THE PARTICLES THAT ARE "ON" (i.e. ABOVE STTHR).
#
# Burning progress averaged in y-dir as a function of time & x-dir
#

BURNX = np.zeros((NSTEPS, NX-1))
#
# SAVE BURNSTAT AT EVERY TIMESTEP
#
BURNEVOL = np.zeros((NSTEPS, NX-1, NY-1))
# ----------------------
# Setup for plotting
# ----------------------
# plt.colorbar()
# DO IT HERE



iplotcount = 0
for istep in range(NSTEPS):
    # Begin timesteps
    #
    # RANDOM WALK & DECAY STATUS TO ACCOUNT FOR MEMORY.
    # IF PARTICLES LEAVE DOMAIN, IGNORE THEM.
    # MOVE ONLY PARTICLES THAT ARE "ON".
    # RANDOM WALK FOLLOWS CNF PAPER.
    #
    # NEW FEATURE: Multiply step by a FACTOR, which reflects the "age" of the
    # particle relative to the clock of the fire in each cell. Hence, 0 causes
    # no propagation, 1 full movement by the wind turbulence (note: in this
    # implementation, FIRECONST = 1). 
    #
    # NEW FEATURE: THE "BUILDING" PARTICLES MODEL RADIATION, HENCE MOVE
    # DIFFERENTLY

    for ip in range(NPART):
        if (FPARTST[ip] > STTHR) and (FPARTTYPE[ip] == 1):
            DU  = -(FPARTUX[ip] - UXM)*2.0*TFREQ*DT + CLANG*SQRTDT*normal()
            DV  = -(FPARTVY[ip] - VYM)*2.0*TFREQ*DT + CLANG*SQRTDT*normal()

            FPARTUX[ip] += DU
            FPARTVY[ip] += DV
            FPARTX [ip] += FPARTUX[ip]*DT*FACTOR[ip]
            FPARTY [ip] += FPARTVY[ip]*DT*FACTOR[ip]
            FPARTST[ip]  = FPARTST[ip]*np.exp(-DT/FPARTMEM[ip])
        elif (FPARTST[ip] > STTHR) and (FPARTTYPE[ip] == 2):
            DU = ULAM*normal()
            DV = ULAM*normal()
            
            FPARTX [ip] = FPARTX[ip] + DU*DT
            FPARTY [ip] = FPARTY[ip] + DV*DT
            FPARTST[ip] = FPARTST[ip]*np.exp(-DT/ TMEMRAD)
        if FPARTX[ip] > XMAX - DX:
            FPARTX [ip] = XMAX - DX
            FPARTST[ip] = 0.0
        elif FPARTX[ip] < XMIN + DX:
            FPARTX [ip] = XMIN + DX
            FPARTST[ip] = 0.0
        if FPARTY[ip] > YMAX - DY:
            FPARTY [ip] = YMAX - DY
            FPARTST[ip] = 0.0
        elif FPARTY[ip] < YMIN + DY:
            FPARTY [ip] = YMIN + DY
            FPARTST[ip] = 0.0
    #  
    # NOW, IGNITE CELLS VISITED BY BURNING PARTICLES, BUT ONLY IF VISITED
    # FOR THE FIRST TIME.  
    # NOTE: WORKS ONLY FOR CONSTANT NUMBER OF PARTICLES (NPARTMAX), OR MUST
    # RE-CALCULATE THE INDEX indp DIFFERENTLY.
    # NOTE: STUPID WAY TO GET INDY & INDX, BUT USE OF "ROUND" DID NOT WORK. WORKS
    # FOR STRUCTURED GRIDS ONLY.
    #
    # To see random walk alone, put FPARTST(ip)>STTHR*10 below
    # to avoid doing any ignitions.
    #
    # 8/3/20: In this version, release all particles from newly-ignited cell
    # after a delay TIGNTR
    #
    # The commented out "elif" can be used to completely stop particles, e.g. from a tall building.
    # Land without flammables is modelled by not releasing new particles.
    #
    for ip in range(NPART):
        if FPARTST[ip] > STTHR:
            for i in range(NX-1):
                if abs(FPARTX[ip] - XCELL[i, 1]) < DX/2:
                    INDX = i
            for j in range(NY-1):
                if abs(FPARTY[ip] - YCELL[1, j]) < DY/2:
                    INDY = j
            BURNPROG[INDX, INDY] += 1
            if (QMAXTR[INDX, INDY] > 0 or QMAXBLD[INDX, INDY] > 0) and BURNSTAT[INDX, INDY] == 0:
                BURNSTAT[INDX, INDY] = 1
                CLOCK   [INDX, INDY] = TIME[istep]
            # elif QMAXTR[INDX, INDY] == 0 or QMAXBLD[INDX, INDY] == 0:
            #     FPARTST[ip] = 0.
            #     FACTOR [ip] = 0.
            # if FPARTTYPE[ip] == 2:
            #     FPARTST[ip] = 0.0
    #
    # NOW, launch new particles from the ignited cells. Randomise a bit the
    # ignition delay time.
    #

    for i in range(NX-1):
        for j in range(NY-1):
            INDX = i
            INDY = j
            TLOCAL = TIME[istep] - CLOCK[INDX, INDY]
            TCRIT = TIGNTR[INDX, INDY] * (1 + RELT*normal())
            if BURNSTAT[INDX, INDY] == 1 and TLOCAL > TCRIT and BURNSTAT2[INDX, INDY] == 1:
                LOCALF = LANGFACTOR
                indp = ((INDX)*(NY- 1) + INDY)*2*NPARTMAX - 1
                for k in range(NPARTTR[INDX, INDY]):
                    FPARTST[k + indp] = 1.0
                    FACTOR [k + indp] = LOCALF
                for k in range(NPARTRAD[INDX, INDY]):
                    FPARTST[k + NPARTTR[INDX, INDY] + indp] = 1.0
                    FACTOR [k + NPARTTR[INDX, INDY] + indp] = LOCALF
                BURNSTAT2[INDX, INDY] = 0
    
    
    ### plots here
    # if istep == 0:
    plt.clf()
    plt.scatter(FPARTX, FPARTY, FPARTST*5 + 0.01, c=FPARTST, cmap="jet")

    plt.pause(0.001)
    iplotcount += 1
    if iplotcount == 10:
        iplotcount = 0
    
    if ITYPESPARK == 1:
        for i in range(NX-1):
            for j in range(NY-1):
                BURNX[istep, i] += BURNPROG[i, j]
            BURNX[istep, i] /= NY - 1
    
    BURNEVOL[istep, :, :] = BURNSTAT[:, :]


print('[END]')
plt.show()