%
clear
%
pd1=makedist('Normal');
pd2=makedist('Normal');
pduni=makedist('Uniform');
%
% DEFINE DOMAIN AND GRID
% IN A REAL APPLICATION, this comes from outside
% UNITS: m
% NX = # of grid lines; (NX-1) is the number of cells
% NY = same in y-dir
%
NX=101;
NY=101;
XMAX=2000.0;
XMIN=0.0;
YMAX=2000.0;
YMIN=0.0;
DX=(XMAX-XMIN)/(NX-1);
DY=(YMAX-YMIN)/(NY-1);
GRIDX = zeros(1,NX);
GRIDY = zeros(1,NY);
for i=1:NX
    GRIDX(i)=XMIN+DX*(i-1);
end
for i=1:NY
    GRIDY(i)=YMIN+DY*(i-1);
end
%
% CELL CENTRES
%
XCELL = zeros(NX-1,NY-1);
YCELL = zeros(NX-1,NY-1);
for i=1:NX-1
    for j=1:NY-1
        XCELL(i,j)=(GRIDX(i+1)+GRIDX(i))/2;
        YCELL(i,j)=(GRIDY(j+1)+GRIDY(j))/2;
    end
end
%
% WIND SPEED AND LANGEVIN PARAMETERS
% All in SI units
% UX = mean wind speed in X-dir (later, it can be f(Z))
% UY = mean wind speed in Y-dir
% A = relative turbulence intensity (u'/<U>) - will be adjusted for proper
%     Langevin walk later
% 
% CAREFUL: this may be a fraction of the real wind speed; given empirical
% data of flame spread vs. wind speed. UPRIME & TTURB ARE TURBULENCE
% QUANTIIES. NEEDS WORK: I AM MIXING PURE AIR QUANTITIES WITH INFERRED 
% MEAN FIRE PROPAGATION SPEED (FIRECONST * AIR SPEED), MUST FIX LATER.
%
% Just putting this down from Wikipedia. Not using it yet.
% BEAUFORT=5;
% WIND SPEED = 0.836*(BEAUFORT)^(3/2)  m/s;
%
% ULAM : a constant to take care of finite propagation with no wind
% TMEMRAD : the "memory" of the radiation particles; calculated as
% Length/ULAM, Length = 10m or something else
%
ULAM = 0.1;
LRAD = 10.0;
TMEMRAD = LRAD/ULAM;
%
%UX=25.0;
%UX=7.77;
UX=25.0;
VY=0.0;
FIRECONST=1.0;
UXM=UX*FIRECONST;
VYM=VY*FIRECONST;
USPEED=sqrt((UXM*UXM)+(VYM*VYM));
A=0.3;
UPRIME=(USPEED*A);
TLENGTH=30.0;
TFREQ=UPRIME/TLENGTH;
CLANG=sqrt(2.0*UPRIME*UPRIME*UPRIME/TLENGTH);
%
% DEFINE PARAMETERS FOR EACH CELL ("PLOT OF LAND")
%
% TREEAREA = fraction of cell area occupied by trees
% BUILDAREA = fraction of cell area occupied by buildings
% TREEAREA + BUILDAREA do not need to sum up to 1.
% CLOCK = the time axis for the potential fire in each cell. 
% NPARTTR = total number of particles emitted from tree fire (0-10)
% NPARTBLD = total number of particles emitted from building fire
% QMAXTR = Maximum heat release (MW/m2) from trees in each cell
% QMAXBLD = Maximum heat release (MW/m2) from buildings in each cell
% TIGNTR = Ignition time delay for tree (since arrival of burning particle
%          in cell) for each cell (in s)
% TENDTR = Fire duration for tree for each cell (s)
% TIGNBLD = As for TIGNTR, but for the building
% TENDBLD = As for TENDTR, but for the building
% BURNSTAT = 0 in the beginning, will become 1 if visited by a burning particle.
% BURNPROG = 0 in the beginning, counts how many burning particles arrived
%            at this cell
% BURNSTAT2 = 1 if cell has all its particles, 0 when it has emptied them.
% 
% Test numbers used for the above for now, more proper parameters and
% their dependence on wind & tree & building type to be used later.
% Other models for tree & building fire can be used, by adjusting NPART
% and their relation with QMAX and Q(t)
%
% INITIALISATIONS
%
% NOT ALL ARE USED IN CURRENT VERSION. MUST WORK FURTHER ON PROPAGATION MODEL 
% 
% TMEM = Particle memory (timescale of decay) in s. (Physically: time of
% flight beyond which gas parcel or firebrand not causing ignitions). Will
% determine the front thickness, to some extent.
%
% STTHR =  Threshold to decide if particle is "on". THIS COULD BE USED TO
% MODEL QUENCHING. OR, if made f(space) to model ignitability of the region.
% 
% LANGFACTOR = Maximum factor (0-1) of the air Langevin walk step that the flame
% particles experience.
% RELT = randomize a bit the ignition time criterion so the propagation is
% not too jumpy
%
NPARTMAX = 1;
TMEM = 10.0;
STTHR = 0.05;
IGNTIME = 10.0;
RELT = 0.2;
BURNDUR = 120;
LANGFACTOR = 0.15;
%
NSTEPS=400;
%
TREEAREA = zeros(NX-1,NY-1);
BUILDAREA = zeros(NX-1,NY-1);
CLOCK = zeros(NX-1,NY-1);
QMAXTR = zeros(NX-1,NY-1);
QMAXBLD = zeros(NX-1,NY-1);
TIGNTR = zeros(NX-1,NY-1);
TENDTR = zeros(NX-1,NY-1);
TIGNBLD = zeros(NX-1,NY-1);
TENDBLD = zeros(NX-1,NY-1);
BURNSTAT = zeros(NX-1,NY-1);
NPARTTR = zeros(NX-1,NY-1);
NPARTRAD =zeros(NX-1,NY-1);
BURNPROG = zeros(NX-1,NY-1);
BURNSTAT2 = zeros(NX-1,NY-1);
%
for i=1:NX-1
    for j=1:NY-1
        TREEAREA(i,j)=0.5;
        BUILDAREA(i,j)=0.5;
        CLOCK(i,j)=0.0;
        QMAXTR(i,j)=10.0;
        QMAXBLD(i,j)=1.0;
        TIGNTR(i,j)=IGNTIME;
        TENDTR(i,j)= TIGNTR(i,j)+BURNDUR;
        TIGNBLD(i,j)=IGNTIME*5;
        TENDBLD(i,j)=TIGNBLD(i,j)+BURNDUR*5;
        BURNSTAT(i,j)=0.0;
        BURNPROG(i,j)=0.0;
        BURNSTAT2(i,j)=1.0;
        NPARTTR(i,j)=NPARTMAX;
        NPARTRAD(i,j)=NPARTMAX;
    end
end
%
% FOR TEST: DEFINE SOME REGIONS THAT CANNOT BURN.
%
%for i=45:55
%  for j=45:55
%    QMAXTR(i,j)=0.0;
%    QMAXBLD(i,j)=0.0;
%  end
%end
%
%NLINE=25;
%for i=NLINE:NX-2
%  j=(NY-1)-(i-NLINE);
%  QMAXTR(i,j)=0.0;
%  QMAXBLD(i,j)=0.0;
%  QMAXTR(i+1,j)=0.0;
%  QMAXBLD(i+1,j)=0.0;
%end
%
% Put a fraction FNON of non-flammable cells; this is to model the
% "vegetation coverage" or concrete buildings & roads.
% 
IFLAG=0;
if (IFLAG==1)
  NNONFL=0.09*(NX-1)*(NY-1);
  IRAND = zeros(1,NNONFL);
  JRAND = zeros(1,NNONFL);
  for ibuild=1:NNONFL
    IRAND=round(random(pduni)*(NX-2)+1);
    JRAND=round(random(pduni)*(NY-2)+1);
    QMAXTR(IRAND,JRAND)=0.0;
    QMAXBLD(IRAND,JRAND)=0.0;
  end
end
%
IFLAG2=1;
if (IFLAG2==1)
  NNONFL=0.5*(NX-1)*(NY-1);
  IRAND = zeros(1,NNONFL);
  JRAND = zeros(1,NNONFL);
  for ibuild=1:NNONFL
    IRAND=round(random(pduni)*(NX-2)+1);
    JRAND=round(random(pduni)*(NY-2)+1);
    QMAXTR(IRAND,JRAND)=0.0;
    QMAXBLD(IRAND,JRAND)=0.0;
  end
% 
% Put a few roads
% 
  NROAD=5;
  IROAD=30;
  for i=IROAD:NX-1
    for j=1:NROAD
      JROAD=round(j*(NY-1)/(NROAD+1));
      QMAXTR(i,JROAD)=0.0;
      QMAXBLD(i,JROAD)=0.0;
%      QMAXTR(i,JROAD+1)=0.0;
%      QMAXBLD(i,JROAD+1)=0.0;    
    end
  end
  for j=1:NY-1
    QMAXTR(IROAD,j)=0.0;
    QMAXBLD(IROAD,j)=0.0;
    QMAXTR(IROAD+1,j)=0.0;
    QMAXBLD(IROAD+1,j)=0.0;    
  end
  for i=40:50
    for j=45:55
      QMAXTR(i,j)=0.0;
      QMAXBLD(i,j)=0.0;
    end
  end
end
%
%surf(XCELL,YCELL,QMAXTR); view([0 90]);
%pause;
%
% PARAMETERS FOR FORWARD INTEGRATION
% UNITS: S
% DT = timestep (in s)
% NTSTEPS = # of timesteps of the simulation
% TIME = the time axis for the whole simulation
%
% NOTE: Proper diffusion process needs care on DT; it depends on # of 
% new particles emitted AND on DX/U. Must experiment. See Appendix of 2012
% CNF paper. Here, I'm playing with the DT as a function of a mean cell transit time
%
%DT=2.0;
if (USPEED > 0)
  DT=2.0*(DX/USPEED);
else
  DT=1.0*(DX/ULAM);
end
SQRTDT=sqrt(DT);
TIME = zeros(1,NSTEPS);
for i=1:NSTEPS
  TIME(i)=i*DT;
end
%
% INITIALISE PARTICLE-RELATED ARRAYS
%
% FPARTST = For each particle, a value 0 or 1 ("cold" or "ignited"). 
%            If "ignited", it moves until it dies according to a decay.
% FPARTMEM = Timescale of the decay (s) for each particle -  for now, this
%            reflects the QMAX
% FPARTH = Height of release (m) for each particle - for a 3-D implementation
% FPARTX = x-position of each particle
% FPARTY = y-position of each particle
%
% FACTOR = From 0 to 1, according to a Gaussian, centred at max burning
% time, with spread so as to account for IGNTIME & BURNDUR.
% In this version, all cells emit equal number of particles, but they
% could have different MEMORY to account for different QMAX
%
% TOTAL NUMBER OF PARTICLES
%
NPART=sum(NPARTTR,'all')+sum(NPARTRAD,'all');
% 
% INITIALISATION
% Use large memory for the time being (all particles stay lit for a long
% time)
% Give all particles some random velocity component, so the drift in the
% later Langevin walk makes sense.
%
FPARTST = zeros(1,NPART);
FPARTMEM = zeros(1,NPART);
FPARTH = zeros(1,NPART);
FPARTX = zeros(1,NPART);
FPARTY = zeros(1,NPART);
FPARTUX = zeros(1,NPART);
FPARTVY = zeros(1,NPART);
FACTOR = zeros(1,NPART);
FPARTTYPE = zeros(1,NPART);
for i=1:NPART
  FPARTMEM(i)=TMEM;
  FPARTUX(i)=UXM+UPRIME*random(pd1)*LANGFACTOR;
  FPARTVY(i)=VYM+UPRIME*random(pd2)*LANGFACTOR;
  FACTOR(i)=0.0;
end
%
% PUT THE PARTICLES IN THE CELLS.
% LOOP OVER CELLS AND DEFINE THEIR PARTICLES.
% FOR NOW, ONLY POSITION DEPENDS ON SPACE; HEIGHT & MEMORY DO NOT.
% FIRST THE TREE PARTICLES, THEN THE BUILDING PARTICLES.
%
icounter=0;
for i=1:NX-1
  for j=1:NY-1
    for k=1:NPARTTR(i,j)
      FPARTX(k+icounter)=XCELL(i,j);
      FPARTY(k+icounter)=YCELL(i,j);
      FPARTTYPE(k+icounter)=1;
    end
    for k=1:NPARTRAD(i,j)
      FPARTX(k+NPARTTR(i,j)+icounter)=XCELL(i,j);
      FPARTY(k+NPARTTR(i,j)+icounter)=YCELL(i,j);
      FPARTTYPE(k+NPARTTR(i,j)+icounter)=2;
    end
    icounter=icounter+NPARTTR(i,j)+NPARTRAD(i,j);
  end
end
%
% START FIRE
%
% ITYPESPARK = 0 (point) ; 1 (line)
%
ITYPESPARK = 1;
%
if (ITYPESPARK == 0) 
%  IF=(NX-1)/2;
  IF=5;
%  JF=(NY-1)/2;
  JF=50;
  BURNSTAT(IF,JF)=1;
  for ip=1:NPART
    if (abs(FPARTX(ip)-XCELL(IF,JF))<DX) && (abs(FPARTY(ip)-YCELL(IF,JF))<DY)
        FPARTST(ip)=1.0;
        FACTOR(ip)=LANGFACTOR;
    end
  end
elseif (ITYPESPARK == 1)
  IF=2;
  for JF=1:(NY-1)
    BURNSTAT(IF,JF)=1;
    for ip=1:NPART
       if (abs(FPARTX(ip)-XCELL(IF,JF))<DX) && (abs(FPARTY(ip)-YCELL(IF,JF))<DY)
         FPARTST(ip)=1.0;
         FACTOR(ip)=LANGFACTOR;
       end
    end
  end
end
%
% RANDOM WALK
% ADVANCE ONLY THE PARTICLES THAT ARE "ON" (i.e. ABOVE STTHR).
%
% Burning progress averaged in y-dir as a function of time & x-dir
%
BURNX = zeros(NSTEPS,NX-1);
%
% SAVE BURNSTAT AT EVERY TIMESTEP
%
BURNEVOL = zeros(NSTEPS,NX-1,NY-1);
% ----------------------
% Setup figure for plotting later
% Savvas 18/03/20 
figure;
colormap(jet); % Try colormap(jet(N)) for a coarser colorbar (N can be anything)
% other colormaps are jet, parula, autumn, spring, hsv and others
colorbar; 
% ----------------------
iplotcount=0;
for istep=1:NSTEPS
%
% Begin timesteps
%
% RANDOM WALK & DECAY STATUS TO ACCOUNT FOR MEMORY.
% IF PARTICLES LEAVE DOMAIN, IGNORE THEM.
% MOVE ONLY PARTICLES THAT ARE "ON".
% RANDOM WALK FOLLOWS CNF PAPER.
%
% NEW FEATURE: Multiply step by a FACTOR, which reflects the "age" of the
% particle relative to the clock of the fire in each cell. Hence, 0 causes
% no propagation, 1 full movement by the wind turbulence (note: in this
% implementation, FIRECONST = 1). 
%
% NEW FEATURE: THE "BUILDING" PARTICLES MODEL RADIATION, HENCE MOVE
% DIFFERENTLY
% 
  for ip=1:NPART
    if (FPARTST(ip)>STTHR) && (FPARTTYPE(ip)==1)
        DU=-(FPARTUX(ip)-UXM)*2.0*TFREQ*DT+CLANG*SQRTDT*random(pd1);
        DV=-(FPARTVY(ip)-VYM)*2.0*TFREQ*DT+CLANG*SQRTDT*random(pd2);
        UXP=FPARTUX(ip)+DU;
        VYP=FPARTVY(ip)+DV;
        FPARTUX(ip)=UXP;
        FPARTVY(ip)=VYP;
        FPARTX(ip)=FPARTX(ip)+UXP*DT*FACTOR(ip);
        FPARTY(ip)=FPARTY(ip)+VYP*DT*FACTOR(ip);
        FPARTST(ip)=FPARTST(ip)*exp(-DT/FPARTMEM(ip));
%        if (FPARTST(ip)<=STTHR) 
%            FPARTST(ip)=0.0;
%        end
    elseif (FPARTST(ip)>STTHR) && (FPARTTYPE(ip)==2)
        DU=ULAM*random(pd1);
        DV=ULAM*random(pd2);
        FPARTX(ip)=FPARTX(ip)+DU*DT;
        FPARTY(ip)=FPARTY(ip)+DV*DT;
%        FPARTX(ip)=FPARTX(ip)+DX*random(pd1);
%        FPARTY(ip)=FPARTY(ip)+DY*random(pd2);
%        phi1=random(pduni)*360;
%        phi2=random(pduni)*360;
%        FPARTX(ip)=FPARTX(ip)+DX*sind(phi1);
%        FPARTY(ip)=FPARTY(ip)+DY*cosd(phi1);
%        FPARTX(ip)=FPARTX(ip)+LRAD*sind(phi1);
%        FPARTY(ip)=FPARTY(ip)+LRAD*cosd(phi1);
        FPARTST(ip)=FPARTST(ip)*exp(-DT/TMEMRAD);
%        FPARTST(ip)=FPARTST(ip)*exp(-DT/5);
%        if (FPARTST(ip)<=STTHR) 
%            FPARTST(ip)=0.0;
%        end
    end
    if (FPARTX(ip)>(XMAX-DX)) 
        FPARTX(ip)=XMAX-DX;
        FPARTST(ip)=0.0;
    elseif (FPARTX(ip)<(XMIN+DX))
        FPARTX(ip)=XMIN+DX;
        FPARTST(ip)=0.0;
    end
    if (FPARTY(ip)>(YMAX-DY)) 
        FPARTY(ip)=YMAX-DY;
        FPARTST(ip)=0.0;
    elseif (FPARTY(ip)<(YMIN+DY))
        FPARTY(ip)=YMIN+DY;
        FPARTST(ip)=0.0;
    end
  end
%  
% NOW, IGNITE CELLS VISITED BY BURNING PARTICLES, BUT ONLY IF VISITED
% FOR THE FIRST TIME.  
% NOTE: WORKS ONLY FOR CONSTANT NUMBER OF PARTICLES (NPARTMAX), OR MUST
% RE-CALCULATE THE INDEX indp DIFFERENTLY.
% NOTE: STUPID WAY TO GET INDY & INDX, BUT USE OF "ROUND" DID NOT WORK. WORKS
% FOR STRUCTURED GRIDS ONLY.
%
% To see random walk alone, put FPARTST(ip)>STTHR*10 below
% to avoid doing any ignitions.
%
% 8/3/20: In this version, release all particles from newly-ignited cell
% after a delay TIGNTR
%
% The commented out "elseif" can be used to completely stop particles, e.g. from a tall building.
% Land without flammables is modelled by not releasing new particles.
%
  for ip=1:NPART
    if (FPARTST(ip)>STTHR) 
%      INDX=round(((FPARTX(ip)-XMIN)/(XMAX-XMIN))*(NX-1));
%      INDY=round(((FPARTY(ip)-YMIN)/(YMAX-YMIN))*(NY-1));
      for i=1:NX-1
        if (abs(FPARTX(ip)-XCELL(i,1))<DX/2)
            INDX=i;
        end
      end
      for j=1:NY-1
        if (abs(FPARTY(ip)-YCELL(1,j))<DY/2)
            INDY=j;
        end
      end
      BURNPROG(INDX,INDY)=BURNPROG(INDX,INDY)+1;
      if ((QMAXTR(INDX,INDY)>0) || (QMAXBLD(INDX,INDY)>0)) && (BURNSTAT(INDX,INDY)==0)
        BURNSTAT(INDX,INDY)=1;
        CLOCK(INDX,INDY)=TIME(istep);
      end
%      elseif (QMAXTR(INDX,INDY)==0) && (QMAXBLD(INDX,INDY)==0)
%          FPARTST(ip)=0.0;
%          FACTOR(ip)=0.0;
%      end
%      if (FPARTTYPE(ip)==2)
%        FPARTST(ip)=0.0;
%      end
    end
  end
%
% NOW, launch new particles from the ignited cells. Randomise a bit the
% ignition delay time.
%
  for i=1:NX-1
    for j=1:NY-1
      INDX=i;
      INDY=j;
      TLOCAL=TIME(istep)-CLOCK(INDX,INDY);
      TCRIT=TIGNTR(INDX,INDY)*(1+RELT*random(pd1));
      if ((BURNSTAT(INDX,INDY)==1) && (TLOCAL > TCRIT) && (BURNSTAT2(INDX,INDY)==1))
        LOCALF = LANGFACTOR;
        indp=((INDX-1)*(NY-1)+INDY-1)*2*NPARTMAX;
        for k=1:NPARTTR(INDX,INDY)
          FPARTST(k+indp)=1.0;
          FACTOR(k+indp)=LOCALF;
        end
        for k=1:NPARTRAD(INDX,INDY)
          FPARTST(k+NPARTTR(INDX,INDY)+indp)=1.0;
          FACTOR(k+NPARTTR(INDX,INDY)+indp)=LOCALF;
        end
        BURNSTAT2(INDX,INDY)=0.0;   
      end
    end
  end
%
% PLOT CURRENT SOLUTION
% 
% ----------------------
% Savvas 18/03/20  
%  scatter (FPARTX,FPARTY,FPARTST*10+1);
  colorVar = FPARTST; colLim = [0 1];
  scatter (FPARTX,FPARTY,FPARTST*10+1,colorVar,'filled');
  colorbar; caxis(colLim);
% ----------------------  
  hold off; 
%  contour(XCELL,YCELL,BURNSTAT);
  pause(0.01);
  iplotcount=iplotcount+1;
  if (iplotcount==10)
%      pause;
      iplotcount=0;
  end
%
  if (ITYPESPARK == 1)
    for i=1:NX-1
      for j=1:NY-1
      BURNX(istep,i)=BURNX(istep,i)+BURNPROG(i,j);
      end
      BURNX(istep,i)=BURNX(istep,i)/(NY-1);
    end
%  plot (BURNX);
%  hold on;
%  pause;
  end
  BURNEVOL(istep,1:NX-1,1:NY-1)=BURNSTAT(1:NX-1,1:NY-1);
%
% End timesteps
% 
end
%
% PLOTTING
%
contour(XCELL,YCELL,BURNSTAT);
figure;
contourf(XCELL,YCELL,BURNPROG);
surf(XCELL,YCELL,BURNPROG); view([0 90]);
pause;
figure;
plot(XCELL(1:99,50),BURNX(50,1:99)); hold on; plot(XCELL(1:99,50),BURNX(100,1:99)); hold on;
plot(XCELL(1:99,50),BURNX(200,1:99)); hold on; plot(XCELL(1:99,50),BURNX(500,1:99));
pause;
figure;
plot(XCELL(1:NX-1,50),BURNEVOL(NSTEPS,1:NX-1,50));
pause;

