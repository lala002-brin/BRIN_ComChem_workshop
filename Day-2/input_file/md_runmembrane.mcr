# YASARA MACRO
# TOPIC:       3. Molecular Dynamics
# TITLE:       Running a molecular dynamics simulation of a membrane protein with normal or fast speed
# REQUIRES:    Dynamics
# AUTHOR:      Elmar Krieger
# LICENSE:     GPL
# DESCRIPTION: This macro sets up and runs a simulation of a membrane protein. It scans the protein for secondary structure elements with hydrophobic surface residues, orients it accordingly and embeds it in a membrane of adjustable lipid composition. Finally a 250 ps restrained equilibration simulation is run, which ensures that the membrane can adapt to the newly embedded protein. Then the real simulation starts.

# Include library functions
include md_library

# Parameter section - adjust as needed, but NOTE that some changes only take
# effect if you start an entirely new simulation, not if you continue an existing one. 
# ====================================================================================

Processors CPUThreads=24,GPU=0
Antialias 0
Console Off

# The structure to simulate must be present with a .pdb or .sce extension.
# If a .sce (=YASARA scene) file is present, the membrane and cell must have been added.
# You can either set the target structure by clicking on Options > Macro > Set target,
# by providing it as command line argument (see docs at Essentials > The command line),
# or by uncommenting the line below and specifying it directly.
#MacroTarget 'c:\MyProject\1crn'
MacroTarget '/home/curie/MyData/agus_ananto/BRIN/MD/4hjo_001.pdb'

# Extension of the cell on each side of the protein in the membrane plane (=XZ plane)
# '15' means that the membrane will be 30 A larger than the protein
memextension=15

# Extension of the cell on each side of the protein along the third (water) axis (=Y-axis)
# '10' means that the cell will be 20 A higher than the protein
waterextension=10

# Flag to use a square membrane. This makes sure that also elongated proteins
# embedded in the membrane can rotate freely during very long simulations. If
# only a short simulation is planned, it can be speeded up by setting the flag
# to 0, creating a rectangular membrane that fits the solute more tightly.
square=1

# Membrane composition: The three letter names of phosphatidyl-ethanolamine (PEA),
# phosphatidyl-choline (PCH, also known as POPC), phosphatidyl-serine (PSE),
# phosphatidyl-glycerol (PGL) and cholesterol (CLR) are each followed by the percentage for
# each membrane side, and must sum up to 100. All lipids are 1-palmitoyl, 2-oleoyl by default.
# The first percentage is for the bottom side of the membrane, the second is for the top side.
# When YASARA shows you the suggested membrane embedding, you need to check that the protein
# orientation matches the membrane composition. If not, flip first and second percentages below and
# rerun the macro. Note that PCH has a large headgroup which cannot form hydrogen bonds, and
# thus reduces membrane stability. PEA is the most stable membrane lipid.
memcomplist()='PEA',100,100,'PCH',0,0,'PSE',0,0,'PGL',0,0,'CLR',0,0

# Or uncomment below to use your own membrane template with 10x10 lipids on each side,
# see membrane simulation recipes for details. (If usermemname='YourChoice', the membrane
# must be saved as yasara/yob/membrane_YourChoice.yob)
usermemname=''
#usermemsize=77.21,73.24

# pH at which the simulation should be run, by default physiological pH 7.4.
ph=7.4

# The ion concentration as a mass fraction, here we use 0.9% NaCl (physiological solution)
ions='Na,Cl,0.9'

# Forcefield to use (this is a YASARA command, so no '=' used)
ForceField AMBER14

# Simulation temperature, which also serves as the random number seed (see Temp command).
# If you increase the temperature significantly by X%, you also need to reduce the timestep by X%
# by changing the 'tslist' that matches your speed below.
temperature='310K'

# Pressure at which the simulation should be run [bar].
pressure=1 

# Cutoff
cutoff=8

# Equilibration period in picoseconds:
# During this initial equilibration phase, the membrane is artificially stabilized
# so that it can repack and cover the solute, while solvent molecules are kept outside.
equiperiod=250

# Delay for animations, 1=maximum speed
delay=100

# The format used to save the trajectories: YASARA 'sim', GROMACS 'xtc' or AMBER 'mdcrd'.
# If you don't pick 'sim', a single *.sim restart file will be saved too, since the other
# two formats don't contain velocities, only positions.
format='sim'

# Duration of the complete simulation, must be longer than equiperiod above.
# Alternatively use e.g. duration=5000 to simulate for 5000 picoseconds
# 'if !count duration' simply checks if variable 'duration' as been defined previously (e.g. by an including macro)
if !count duration
  duration='100000'

# The simulation speed, either 'slow' (2*1 fs timestep), 'normal' (2*1.25 fs timestep) or
# 'fast' (maximize performance with 2*2.5 fs timestep and constraints)
# Do not use 'fast' if you simulate incorrect molecules (that would not be stable in reality) 
# 'if !count speed' simply checks if variable 'speed' as been defined previously (e.g. by an including macro) 
if !count speed
  speed='normal'

# The save interval for snapshots. Normally you don't need more than 500-1000 snapshots
# of your simulation, since that's the resolution limit of a typical figure in a journal.
if speed=='fast'
  # Fast speed, save simulation snapshots every 250000 fs, i.e. 250 ps.
  saveinterval=250000
else  
  # Slow or normal speed, save simulation snapshots every 100000 fs, i.e. 100 ps.
  saveinterval=100000

# Flag if the protein's membrane embedding needs to be confirmed
confirmneeded=1

# Normally no change required below this point
# ============================================

RequireVersion 15.1.1

# Treat all simulation warnings as errors that stop the macro
WarnIsError On

# Membrane simulations are always periodic
Boundary periodic

# Do we have a target?
if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target structure"

# When run as a macro in text mode, add configuration details to log file
if runWithMacro and ConsoleMode
  Processors

Clear
Console off
SimSteps 1,1
# Do we already have a scene with water?
waterscene = FileSize (MacroTarget)_water.sce
if waterscene
  LoadSce (MacroTarget)_water
else
  # No scene with water present yet
  # Do we have a scene with the protein embedded in the membrane?
  scene = FileSize (MacroTarget).sce
  if scene
    # Yes, protein with membrane is already present
    LoadSce (MacroTarget)
    # Search for the membrane
    c = CountObj Membrane
    if c!=1
      # No membrane object, this must be a user-provided initial scene to prevent modifications like Clean/OptHyd, load later
      scene=0
  if !scene
    # No membrane scene present yet
    # Has the user already oriented the protein inside the membrane interactively?
    oriscene = FileSize (MacroTarget)_ori.sce
    if oriscene
      LoadSce (MacroTarget)_ori
      Unselect
    else
      # Load the protein, assuming it's a SCE, PDB or YOB file
      for filetype in 'sce','yob','pdb'
        size = FileSize (MacroTarget).(filetype)
        if size
          break
      if !size
        RaiseError 'Initial structure not found, expected (MacroTarget).pdb or .yob. Make sure to create a project directory and place the structure there in PDB or YOB format'
      # Load structure
      Load(filetype) (MacroTarget)
      DelObj SimCell
      nameclash = CountObj Membrane
      if nameclash
        NameObj all,Solute
      # In case user accidentally provided a YOb file with selected atoms
      Unselect
      # Orient the protein in the membrane
      OriInMem 1,(filetype!='sce'),ph,waterextension,memextension,square,0,confirmneeded,delay
      SaveSce (MacroTarget)_ori
    # Next step is to build the membrane
    AutoPosObj 1,X=-200,Z=100,Steps=(delay),Wait=No
    AutoPosObj MembPreview,X=200,Z=100,Steps=(delay),Wait=No
    AutoPosObj SimCell,X=0,Y=200,Z=100,Steps=(delay)
    DelObj not SimCell
    membranename = BuildMem MacroTarget,memcomplist,usermemname,delay,0
    # Load the orientation scene again and include the membrane
    LoadSce (MacroTarget)_ori
    # Insert the protein into the membrane
    InsertIntoMem 1,membranename,0,0,''
    SaveSce (MacroTarget)
  # Solvate system with membrane, get membrane Y-coordinate
  SolvateMem ph,ions,cutoff
  # Save scene with water
  SaveSce (MacroTarget)_water
  HideMessage
Wait 1
# Don't keep selected atoms, LoadXTC/LoadMDCRD would load only selected ones
Unselect

# Choose timestep and activate constraints
if speed=='fast'
  # Fast simulation speed
  # Constrain bonds to hydrogens
  FixBond all,Element H
  # Constrain certain bond angles involving hydrogens
  FixHydAngle all
  # Choose a multiple timestep of 2*2.5 = 5 fs
  # For structures with severe errors, 2*2 = 4 fs is safer (tslist=2,2)
  tslist=2,2.5
else
  # Slow or normal simulation speed
  # Remove any constraints
  FreeBond all,all
  FreeAngle all,all,all
  if speed=='slow'
    # Choose a multiple timestep of 2*1.00 = 2.0 fs
    tslist=2,1.0
  else
    # Choose a multiple timestep of 2*1.25 = 2.5 fs
    tslist=2,1.25
    # With this timestep, atoms may get too fast in very rare circumstances (only
    # in a specific protein, only once every few nanoseconds). The command below
    # slows down atoms moving faster than 13000 m/s. Such a 'random collision' every
    # few nanoseconds has no more impact than the random number seed. You can comment
    # it out for most proteins, or use the smaller timestep with speed 'slow' above:
    Brake 13000
# During equilibration update the pairlist every 10 steps
SimSteps Screen=10,Pairlist=10
# Calculate total timestep, we want a float, so tslist2 is on the left side
ts=tslist2*tslist1
# Snapshots are saved every 'savesteps'
savesteps=saveinterval/ts
# Set final simulation parameters
TimeStep (tslist)
Temp (temperature)
Cutoff (cutoff)
Longrange Coulomb
# Make sure all atoms are free to move
FreeAll
# Alread a snapshot/trajectory present?
i=00000
if format=='sim'
  trajectfilename='(MacroTarget)(i).sim'
else  
  trajectfilename='(MacroTarget).(format)'
  restartfilename='(MacroTarget).sim'  
  # Backwards compatibility: Starting with YASARA version 12.8.1, XTC trajectories no longer contain a number in the filename
  old = FileSize (MacroTarget)(i).xtc
  if old
    RenameFile (MacroTarget)(i).xtc,(trajectfilename)
running = FileSize (trajectfilename)
if not running
  # Perform energy minimization
  Experiment Minimization
  Experiment On
  Wait ExpEnd
  # And now start the real simulation
  Sim On
else
  # Simulation has been running before
  ShowMessage "Simulation has been running before, loading last snapshot..."
  # Switch console off to load the snapshots quickly
  Console Off
  if format=='sim'
    # Find and load the last SIM snapshot
    do
      i=i+1
      found = FileSize (MacroTarget)(i).sim
    while found
    i=i-1
    LoadSim (MacroTarget)(i)
    # Adjust savesteps to save snapshots in the same interval as previously
    if i>0
      t = Time
      savesteps=0+t/(ts*i)
  else
    # Do we have a restart file with atom velocities?
    found = FileSize (restartfilename)
    if found
      # Yes. First determine the savesteps if possible by loading the 2nd XTC/MDCrd snapshot
      last,t = Load(format) (trajectfilename),1
      if !last
        last,t = Load(format) (trajectfilename),2
        savesteps=0+t/ts
      # Then load the restart file
      LoadSim (restartfilename)
    else
      # No restart file found, load the last snapshot in the XTC/MDCrd trajectory
      do
        i=i+1
        last,t = Load(format) (trajectfilename),(i)
        ShowMessage 'Searching (format) trajectory for last snapshot, showing snapshot (i) at (0+t) fs'
        Sim Pause
        Wait 1
      while !last
      savesteps=0+t/(ts*(i-1))
      Sim Continue
HideMessage

# Set temperature and pressure control
TempCtrl Rescale
PressureCtrl Manometer2D,Pressure=(pressure)

# Now the simulation is running, here you can make changes to the force field

# Uncomment to fix certain atoms in space
# FixAtom Backbone Mol B

# Uncomment to add distance constraints
# AddSpring O Res Lys 80,H Res Glu 84,Len=1.9
# AddSpring 2204,4667,Len=1.9,SFC=100
# Uncomment to modify charges, e.g. let Trp 12 in Mol A lose an electron:
# ChargeRes Trp 12 Mol A,+1

# And finally, make sure that future snapshots are saved
Save(format) (trajectfilename),(savesteps)
if format!='sim'
  # We save an XTC/MDCrd trajectory plus a single Sim restart file
  SaveSim (restartfilename),(savesteps),Number=no

CorrectDrift off
    
# At the beginning of the simulation, the membrane is not ideally packed yet and needs
# about 250ps equilibration time. During this period it is still 'vulnerable' to
# water molecules that are squeezed in. We thereforce keep these waters out.
t = Time
if t<equiperiod*1000
  EquilibrateMem equiperiod,ts

# Membrane has been equilibrated, now keep the membrane protein from diffusing around and crossing periodic boundaries
CorrectDrift on

# After equilibration update the pairlist every 10 (CPU) or 25 (GPU) steps
_,_,gpu = Processors
if gpu
  SimSteps Screen=25,Pairlist=25
else    
  SimSteps Screen=10,Pairlist=10

if duration=='forever'
  Console On
  if ConsoleMode
    # In the console, we need to wait forever to avoid a prompt for user input
    Wait forever
else
  Console Off
  measurements=0
  # Wait for given number of picoseconds
  do
    # Tabulate properties you want to monitor during the simulation,
    # e.g. the speeds and velocity vectors of atoms 4, 5 and 7:
    #Tabulate SpeedAtom 4 5 7
    # Or apply a pulling force, in this example downwards
    #AccelRes GLI 1027,Y=-10000
    if ConsoleMode
      # In console model, display the ns/day
      Console On
    # Wait for one screen update
    Wait 1
    measurements=measurements+1
    t = Time
  while t<1000.*duration+1
  # Did we create a table with measurements?
  vallist() = Tab Default
  if count vallist
    # Yes, save the table
    SaveTab default,(MacroTarget)_duringsim,Format=Text,Columns=(count vallist/measurements),Header='Insert your own header here'
  Sim Off
# Exit YASARA if this macro was provided as command line argument in console mode and not included from another macro
if runWithMacro and ConsoleMode and !IndentationLevel
  Exit
