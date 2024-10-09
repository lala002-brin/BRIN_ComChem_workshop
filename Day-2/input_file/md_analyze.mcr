# YASARA MACRO
# TOPIC:       3. Molecular Dynamics
# TITLE:       Analyzing a molecular dynamics trajectory
# REQUIRES:    Dynamics
# AUTHOR:      Elmar Krieger and Kornel Ozvoldik
# LICENSE:     GPL
# DESCRIPTION: This macro analyzes a simulation and creates a detailed report with a large number of plots, e.g. energies, RMSDs, hydrogen bonds. It also tries to identify the main ligand and provides ligand-specific data. All results are additionally written to a simple text table, which can be imported into your favorite spreadsheet program. Your own analysis can often be added with just one line of code, search for 'Example:'.

RequireVersion 20.1.1

# MD report initialization parameters and flags
# =============================================

Processors CPUThreads=8,GPU=0
Antialias 0
Console Off



# The structure to analyze must be present with a .sce extension.
# You can either set the target structure by clicking on Options > Macro > Set target,
# by providing it as command line argument (see docs at Essentials > The command line),
# or by uncommenting the line below (=remove the '#') and specifying it directly.
#MacroTarget 'c:\MyProject\1crn'
MacroTarget /mgpfs/home/wsbrinapctp03/workshop_yasara/4hjo_001.pdb,Remove=Extension

# Set common beginning for all result filenames. By default,
# this is the same as the macro target, but you can change
# it to run multiple analyses at the same time.
resultbase=MacroTarget            # Default
#resultbase='(MacroTarget)_run1'  # Example

# Forcefield to use for analysis, should be the same as the one used to run the simulation
ForceField AMBER14,SetPar=Yes   # Default
#ForceField YASARA2,SetPar=Yes  # Example: Add a quality Z-score in YASARA Structure

# Number of the solute object whose RMSDs from the starting conformation will be calculated
# If the protein is an oligomer, check the documentation of the 'Sup' command at 'analyzing a simulation' to avoid pitfalls.
soluteobj=1

# Flag to convert the entire trajectory to PDB format (solute object only)
pdbsaved=0

# The B-factors calculated from the root-mean-square fluctuations can be too large to fit them
# into the PDB file's B-factor column. Replace e.g. 1.0 with 0.1 to scale them down to 10%
bfactorscale=1.0

# Trajectory block to be analyzed. The 'if not count block' skips this part if this macro is included
# by the md_analyzeblock macro, that analyzes the trajectory in blocks (see 'Analyzing a trajectory' in the docs).
if not count block
  # First snapshot to be analyzed, increase number to ignore an equilibration period.
  firstsnapshot=0
  # Number of snapshots to be analyzed
  snapshots='all'
  # Set snapshotstep > 1 if you don't want to analyze every snapshot.
  # (e.g. snapshotstep=10 means that only every 10th snapshot will be analyzed)
  snapshotstep=1
  
# All snapshots will be superposed on this reference snapshot to calculate RMSDs etc.
# The starting structure is snapshot 0. Having run the macro once, you can also change
# refsnapshot=X to refsnapshot='average' to superpose on the time average structure.
refsnapshot=0

## In case you want to use commands like GroupCenter inside a periodic cell, you will
# get a warning about unexpected results due to periodic boundary effects. In this
# case set 'central' to an atom in the middle of your region of interest. This atom will
# be kept at the center of the cell, which will hopefully also keep your region of interest
# away from the cell boundaries. Then you can safely ignore the warning.
central=0

# In case you want to cluster the trajectory, set the minimum heavy-atom RMSD between different clusters.
# A representative structure for each cluster will be saved as YourProtein_cluster*.yob,
# all representative structures will be at least 'rmsdmin' Angstrom apart.
rmsdmin=0.0

# Width of figures in pixels
figurewidth=1024

# Flag to save high resolution versions of plots in 4:3 format
hiresplotted=1

# Set maximum number of table rows that will be shown in the report
tabrowsmax=10     # Default: Show only first 10 rows of each table
#tabrowsmax='all' # Example: Show the complete tables

# Selection of atoms to include for 'Calpha' RMSD calculation
# (includes C1* of nucleic acids to consider DNA/RNA)
casel='CA Protein or C1* NucAcid'

# Ligand atom selection for RMSD analysis
ligandsel=''                   # Default: Identify the ligand automatically
#ligandsel='Res ATP 501 Mol A' # Example: Select residue ATP 501 in molecule A as ligand
#ligandsel='Atom Zn'           # Example: Select Zinc ion as ligand
# Select solute residues as ligand to analyze e.g. protein domain movement or protein-protein contacts of a complex
#ligandsel='Mol A'             # Example: Select protein Mol A of an oligomeric solute as ligand
#ligandsel='Res 1-50'          # Example: Select solute protein residue 1 to 50 as ligand

# Atom selection to list individual hydrogen bonds
hbosel=''                      # Default: Choose ligand, i.e. list ligand hydrogen bonds
#hbosel='Res Arg 17 Mol A'     # Example: List H-bonds of residue Arg 17 in molecule A

# Maximum number of hydrogen bonds that will be listed
hbondsmax=30                   # Default: List only the first 30 hydrogen bonds
#hbondsmax='all'               # Example: List all hydrogen bonds
#hbondsmax=0                   # Example: Don't list hydrogen bonds

# Set minimum and/or maximum residue number which will be included in the residue RMSF and per-residue plots.
resnummin=-999
resnummax=9999

# Selection to calculate and visualize the dynamic cross-correlation matrix (DCCM)
dccmsel='Atom CA Protein or C1* NucAcid' # Default: Calculate the DCCM for protein Calpha atoms or nucleic acid C1* atoms
#dccmsel='Res Protein'                   # Example: Calculate the DCCM for protein residue centers
#dccmsel=''                              # Example: Don't calculate the DCCM

# Set the color of the minimum and maximum value of the DCCM visualization.
dccmcol()='blue','yellow'

# Save a file *_dccm.yob, which visualizes the DCCM by joining atoms with
# red (correlation>=dccmcut) and blue (correlation<=-dccmcut) lines.
dccmcut=0.9

# Method for calculation of the electrostatic potential(ESP) in the cell. The ESP for each snapshot is saved as a
# Gaussian cube map and then visualized automatically by the md_play macro. The grid resolution is set with the SolvPar command.
espmethod=''                             # Default: Don't calculate the ESP
#espmethod='PBS'                         # Example: Runs the APBS program to solve the Poisson-Boltzmann equation (slow)
#espmethod='PME'                         # Example: Uses Particle Mesh Ewald electrostatic potential (fast)

# Selection to calculate the radial distribution function (RDF)
# Syntax: rdfsellist()='Atom selection 1','Atom selection 2',Bins,BinWidth
# Note that you may have to save more snapshots than usually to avoid problems with sparse
# data and noisy RDF results.
# Default: Don't calculate an RDF
rdfsellist()=''
# Example: Calculate the RDF of water in 40 bins, each 0.25 A wide (thus up to 10 A).
#rdfsellist()='O Res HOH','O Res HOH',40,0.25
# Example: Calculate the RDF between two specific atoms in 20 bins, each 0.5 A wide
#rdfsellist()='CG Res Asp 120','ND1 Res His 200',20,0.5

# Definition of analyses to perform
# =================================
#
# Please see the user manual at Recipes > Run a molecular dynamics simulation > Analyzing a trajectory
# (scroll to the end) for detailed instructions how to add your own analyses.

# Define the analyses to perform inside the simulation cell, considering periodic boundaries
# ==========================================================================================
def AnalyzeInsideCell
  global ligandsel,hbosel,hbondsmax,resnummin,resnummax
  if hbosel==''
    hbosel=ligandsel
  
  # Plot simulation cell lengths
  celllist1,celllist2,celllist3=Cell
  Plot celllist,'Simulation cell lengths','Length in Angstrom','CellLengthX CellLengthY CellLengthZ'
  
  # Plot the total energy...
  elist()=EnergyAll All
  Plot sum elist,'Total potential energy of the system','Energy in (EnergyUnit)','TotalEnergy'
  # ...and the individual components
  Plot elist,'Potential energy components','Energy in (EnergyUnit)','Bond Angle Dihedral Planarity Coulomb VdW Packing1D Packing3D'
  
  # Plot the VdW, molecular and solvent accessible surfaces of the solute
  Plot "SurfObj Solute",'Surface areas of the solute','Surface areas in Angstrom^2','SurfVdW SurfMol SurfAcc'

  # Plot the number of hydrogen bonds inside the solute
  # We need to divide by two since bonds are listed in both directions 
  hbolist()=ListHBoAtom Obj Solute,Obj Solute
  Plot count hbolist/2,'Number of hydrogen bonds in the solute','Hydrogen bonds','SoluteHBonds'

  # Plot the number of hydrogen bonds between solute and solvent
  hbolist()=ListHBoAtom Obj Solute,Obj Solvent
  Plot count hbolist,'Number of hydrogen bonds between solute and solvent','Hydrogen bonds','SltSlvHBonds'
  
  # Plot the number of hydrogen bonds of the selected atoms (Default: ligand)
  if hbosel!=''
    # Internal bonds are counted twice, so the total number stays the same if alternative bonds are formed
    hbosellist()=ListHBoAtom (hbosel),(hbosel)
    hboreclist()=ListHBoAtom Obj Solute and not (hbosel),(hbosel)
    hboslvlist()=ListHBoAtom Obj Solvent Membrane,(hbosel)
    hbolist()=count hbosellist,count hboreclist,count hboslvlist
    hbolist4=sum hbolist
    Plot hbolist,'Number of hydrogen bonds made by (hbosel)','Hydrogen bonds','SelHBonds RecHBonds SlvHBonds TotHBonds'
  
  # Tabulate all hydrogen bonds made by the selected atoms (Default: ligand) excluding waters.
  # We first need to estimate the maximum number of H-bonds so that the table can be padded accordingly.
  # For every H-bond, the table then contains atom names, energy and hydrogen-acceptor distance.
  if hbosel!='' and hbondsmax!=0
    acceptors=CountAtom Element N O S P and (hbosel)
    donors=CountAtom Element H with bond to Element N O S P and (hbosel)
    bonds=acceptors*2+donors
    if hbondsmax=='all' or hbondsmax>bonds
      hbondsmax=bonds
    if hbondsmax
      hbolist()=ListHBoAtom (hbosel),!Water,Results=6
      hbonds=count hbolist/6
      hboname=''
      for j=1 to hbondsmax
        hboname=hboname+'HB(j)Atm1 HB(j)Atm2 HB(j)E HB(j)D '
        if j<=hbonds
          # Get data of atom 1 and atom 2.
          for k=0 to 1
            hboresultlist(4*j-3+k)=ListAtom (hbolist(j*6-5+k)),Format='ATOMNAME.RESNAME1RESNUM.MOLNAME'
          if hbolist(j*6-2)=='acc'
            # First listed atom is always the donor
            swap=hboresultlist(4*j-3)
            hboresultlist(4*j-3)=hboresultlist(4*j-2)
            hboresultlist(4*j-2)=swap
          # Get energy and distance.
          hboresultlist(4*j-1)=hbolist(j*6-3)
          hboresultlist(4*j)=hbolist(j*6)
        else
          # If no bond then fill cells with '-' to facilitate parsing.
          for k=1 to 4
            hboresultlist(4*j-4+k)='-'
      # Write data to table in report.
      WriteTable hboresultlist,'Hydrogen bonds made by (hbosel)','(hboname)'
  
  # Plot the secondary structure content
  proteins=CountMol Protein
  if proteins
    Plot "SecStr",'Protein secondary structure content','Secondary structure fractions in percent','Helix Sheet Turn Coil Helix310 HelixPi'
  
  # Plot a per-residue map of the protein secondary structure
  sel='Protein Obj Solute and Res (resnummin)-(resnummax)'
  # Define secondary structure types to analyze, must match the SecStr command
  secstrnamelist()=  'Helix','Sheet','Turn','Coil','Helix310','HelixPi'
  secstrletterlist()='H',    'E',    'T',   'C',   'G',       'I'
  # Make a letter->value dictionary and get default secondary structure colors
  for i=1 to count secstrletterlist
    secstrnum(secstrletterlist(i))=i
    secstrcolorlist(i)=ColorPar SecStr,(secstrnamelist(i))
  # Create the residue list
  secstrreslist()=ListRes (sel) and SecStr (join secstrnamelist)
  if count secstrreslist
    # Get the secondary structure of each residue
    secstrlist()=SecStrRes (join secstrreslist())
    for i=1 to count secstrlist
      secstrlist(i)=secstrnum(secstrlist(i))
    # Plot per-residue data as a function of simulation time
    PlotRes 'secstr',secstrreslist,secstrlist,1,6,'Per-residue protein secondary structure','(join secstrnamelist)','(join secstrcolorlist)'
  
  # Plot a per-residue map of the number of contacts
  # Create a solute residue list
  conreslist()=ListRes Protein or NucAcid and Obj Solute Res (resnummin)-(resnummax)
  if count conreslist
    # Count number of contacts for each residue
    conlist()=CountConRes (join conreslist),Obj Solute,Cutoff=0.5,Subtract=HBoRadii,Exclude=5,Occluded=No
    # Plot per-residue number of contacts as a function of simulation time
    PlotRes 'con',conreslist,conlist,1,15,'Per-residue number of contacts','','Blue Yellow'
  
  # Plot a per-residue map of interactions of the receptor with the ligand
  if ligandsel!=''
    sel='Protein or NucAcid and Obj Solute Res (resnummin)-(resnummax) and not (ligandsel)'
    ligintreslist()=ListRes (sel)
    if count ligintreslist
      # Create a list of empty interaction values (=0)
      for i=1 to count ligintreslist
        vallist(i)=0
      # Loop over the interactions to analyze, ignoring occluding hydrophobic interactions
      intlist='HBonds','Hydrophobic','Ionic'
      occlist='',      'No',         'Yes'
      for i=1 to count intlist
        # Get the list of residues interacting with the ligand
        if i==1
          intreslist() = ListHBoRes (sel),(ligandsel),Results=1
        else
          intreslist() = ListIntRes (sel),(ligandsel),Type=(intlist(i)),Exclude=5,Occluded=(occlist(i)),Results=1
        idx=1
        intresidues=count intreslist
        for j=1 to count ligintreslist
          # Stop if we tagged all interacting residues
          if idx>intresidues
            break
          if ligintreslist(j)==intreslist(idx)
            # Tag residue j with interaction i
            vallist(j)=vallist(j)|(1<<(i-1))
            idx=idx+1
      # Increment vallist, the heatmap starts at 1
      vallist()=vallist+1
      # Now vallist contains values from 1 (no interactions)..8 (all three interactions).
      # Each interaction is an RGB color component: HBonds=Red,Hydrophobic=Green,Ionic=Blue, so the legend entries are
      # HB  Hyd   Hyd+HB Ion  Ion+HB  Ion+Hyd Ion+Hyd+HB , and the corresponding colors are
      # Red Green Yellow Blue Magenta Cyan    Gray
      PlotRes 'ligresint',ligintreslist,vallist,1,8,'Per-residue ligand interactions of the receptor','None HB Hyd Hyd+HB Ion Ion+HB Ion+Hyd Ion+Hyd+HB','White Red Green Yellow Blue Magenta Cyan Gray'

  # Plot a per-atom map of interactions of the ligand with the receptor
  if ligandsel!=''
    sel='Protein or NucAcid and Obj Solute Res (resnummin)-(resnummax) and not (ligandsel)'
    ligatomlist()=ListAtom (ligandsel)
    if count ligatomlist
      # Create a list of empty interaction values (=0)
      for i=1 to count ligatomlist
        vallist(i)=0
      # Loop over the interactions to analyze, ignoring occluding hydrophobic interactions
      intlist='HBonds','Hydrophobic','Ionic'
      occlist='',      'No',         'Yes'
      for i=1 to count intlist
        # Get the list of ligand atoms interacting with the receptor
        if i==1
          intatomlist() = ListHBoAtom (ligandsel),(sel),Results=1
        else
          intatomlist() = ListIntAtom (ligandsel),(sel),Type=(intlist(i)),Exclude=5,Occluded=(occlist(i)),Results=1
        idx=1
        intatoms=count intatomlist
        for j=1 to count ligatomlist
          # Stop if we tagged all interacting atoms
          if idx>intatoms
            break
          if ligatomlist(j)==intatomlist(idx)
            # Tag atom j with interaction i
            vallist(j)=vallist(j)|(1<<(i-1))
            idx=idx+1
      # How vallist values correspond to interactions see 'Plot a per-residue map of interactions with the ligand' above
      vallist()=vallist+1
      # A screen shot of the ligand with atom numbers will be added to the report (see 'PrintLigAtomIntInfo' below)
      PlotAtom 'ligatomint',ligatomlist,vallist,1,8,'Per-atom receptor interactions of the ligand','None HB Hyd Hyd+HB Ion Ion+HB Ion+Hyd Ion+Hyd+HB','White Red Green Yellow Blue Magenta Cyan Gray'
  
  # The following examples provide a few hints for other things to analyze while
  # the simulation is active (i.e. considering periodic boundaries).
  # Just delete the '#' character to uncomment the command(s) of each example.

  # Example: Plot the total energy of the ligand and the individual components
  #elist()=EnergyRes (ligandsel),All
  #Plot sum elist,'Total potential energy of the ligand','Energy in (EnergyUnit)','LigTotEnergy'
  #Plot elist,'Potential energy components of the ligand','Energy in (EnergyUnit)','LigBond LigAngle LigDihedral LigPlanarity LigCoulomb LigVdW LigPacking1D LigPacking3D'
  
  # Example: Plot the distance between the carboxyl group of Glu 123 (Cdelta) in molecule/chain A
  #          and the guanidinium group of Arg 345 (Czeta) in molecule/chain B:
  # NOTE: This distance measurement considers periodic boundaries and may yield unexpected results
  #       for distances larger than half the smallest cell side length. To measure such large distances,
  #       look for the Distance example further below in the 'AnalyzeOutsideCell' section.
  #Plot "Distance CD Res Glu 123 Mol A, CZ Res Arg 345 Mol B",'Distance between Glu 123 [Cdelta] and Arg 345 [Czeta]','Distance','DGLUARG'
  
  # Example: Plot the angle between three Calpha atoms:
  #Plot "Angle CA Res Glu 10, CA Res Asp 30, CA Res Lys 40",'Angle between three Calpha atoms','Angle','Phi'
  
  # Example: Plot the number of salt-bridges involving Lys,Arg and Asp,Glu
  #lysbridges=CountRes Lys Atom NZ with distance<4 from Asp Glu Atom OD? OE?
  #argbridges=CountRes Arg Atom NE NH? with distance<3.5 from Asp Glu Atom OD? OE?
  #Plot (lysbridges+argbridges),'Number of salt-bridges involving Lys,Arg and Asp,Glu','Salt-bridges','Bridges'

  # Example: Plot the number of residues involved in hydrophobic interactions with Phe 13
  #intlist()=ListIntRes Phe 13,all,Type=Hydrophobic,Exclude=5
  #Plot count intlist,'Residues in hydrophobic interactions with Phe 13','Interactions','IAResPhe13'
  
  # Example: Tabulate ionic interactions inside the solute in an extra tab file
  #intlist()=ListIntRes Obj Solute,Obj Solute,Type=Ionic,Exclude=5,Results=2
  #result=''
  #for i=1 to count intlist step 2
  #  pair()=ListRes (intlist(i)) (intlist(i+1)),Format=RES
  #  result=result+'(pair(1))-(pair(2)),'
  #WriteExtraTable 'IONIC',result,'Ionic'

  # Example: Tabulate the distance and interaction strength between ligand
  # and receptor atoms participating in pi-pi interactions in an extra tab file.
  #intlist()=ListIntAtom (ligandsel),Obj Solute,Type=PiPi,Exclude=5,Results=5
  #result=''
  #for i=1 to count intlist step 5
  #  for j=1 to 2
  #    atom(j)=ListAtom (intlist(i+j-1)),Format='ATOMNUM.ATOMNAME.RESNAME3RESNUM.MOLNAME'
  #  result=result+'|(atom1),(atom2),(0.00+intlist(i+2)),(0.00+intlist(i+4))'
  #WriteExtraTable 'PIPIATOM',result,'A1,A2,D[A],S'
  
  # Example: Plot the potential energy of residue Glu 123:
  #Plot "EnergyRes Glu 123",'Potential energy of residue Glu 123','Energy','EpotGLU123'

  # Example: Plot the solvent accessible surface areas of methionine sulfur atoms to predict their
  # susceptibility for oxidation. Add entire solute to surface environment for partial surface calculation.
  #RemoveEnvAll
  #AddEnvObj Solute
  #reslist()=ListRes Met Obj Solute,Format=RESNAME1RESNUM
  #Plot "SurfAtom SD Res Met Obj Solute,Type=Accessible,Unit=Atom",'SASA of Met S-delta atoms','SASA in A^2','(join reslist)'
  
  # Example: Plot the solvent accessible surface area and the electrostatic surface potential of the residues Thr 1, Ser 11 and Arg 17
  #ressel='Thr 1 Ser 11 Arg 17'
  #RemoveEnvAll
  #AddEnvObj Solute
  #reslist()=ListRes (ressel),Format='RESNAME1RESNUM.MOLNAME'
  #datalist()=SurfESPRes (ressel),accessible,PME,Res
  #for i=1 to count reslist
  #  saslist(i)=datalist(2*i-1)
  #  esplist(i)=datalist(2*i)
  #Plot saslist,'Solvent accessible surface area of solute residues','SASA in A^2','(join reslist)'
  #Plot esplist,'Electrostatic surface potential of solute residues','ESP in (EnergyUnit)','(join reslist)'
  
  # Example: Plot a per-residue map and the mean per-residue solvent accessible surface area
  #sasreslist()=ListRes Obj Solute Res (resnummin)-(resnummax)
  #RemoveEnvAll
  #AddEnvObj Solute
  #saslist()=SurfRes (join sasreslist),Accessible,Res
  #PlotRes 'sasres',sasreslist,saslist,'','','Per-residue SASA in A^2','Min Max','Blue Yellow'
  
  # Example: Plot the distance from the nearest water molecule to the atoms selected in dissel1 and dissel2
  #dissel1='C Res Gly 2'
  #dissel2='NE2 Res Hid 254'
  #oatom=ListAtom Element O Res HOH with distance<3.2 from (dissel1) and with distance<3 from (dissel2)
  #for i=1 to 2
  #  if oatom
  #    distance(i)=Distance (oatom),(dissel(i))
  #  else
  #    distance(i)=5
  #Plot (distance),'Water distance to 1 = (dissel1) and 2 = (dissel2)','Distance in Angstrom','WaterDis1 WaterDis2'
  
  # Example: Count water molecules inside a membrane pore
  #poreressel='Ile 263'
  #_,bottom=GroupCenter Element P Segment MEM1
  #_,top=GroupCenter Element P Segment MEM2
  #bottom=bottom+2
  #top=top-2
  #olist()=ListAtom Element O Res HOH LocalY<(top) LocalY>(bottom) with distance<3.5 from Res (poreressel)
  #do
  #  waters=count olist
  #  olist()=ListAtom (join olist) or Element O Res HOH LocalY<(top) LocalY>(bottom) with distance<3.5 from (join olist)
  #while waters<count olist
  #Style Ball
  #StickRes not Atom (join olist)
  #ColorRes Atom (join olist),Blue
  #Plot count olist,'Number of water molecules in membrane pore (poreressel)','Number of water molecules','PoreWaters'
  
  # Example: Tabulate the distance of every protein residue to the common center of all proteins
  #solcen()=GroupCenter Protein Obj Solute
  #residlist()=ListRes Protein Obj Solute
  #resnamelist()=ListRes Protein Obj Solute,Format='RESNAME1RESNUM.MOLNAME'
  #for i=1 to count residlist
  #  rescen=GroupCenter (residlist(i))
  #  dis(i)=norm (solcen-rescen)
  #WriteTable (dis),'Distance of protein residues to the solute center','(join resnamelist)'
  
  # Example: Tabulate the distances of the nearest 10 waters less than 5 A from Res Glu 123
  #sel='Res Glu 123'
  #dis=5
  #distances=10
  #dislist()=Distance (sel),O Water with distance<(dis) from (sel)
  #dislist()=sort dislist
  #disname=''
  #for i=1 to distances
  #  if i>count dislist
  #    disresultlist(i)='-'
  #  else
  #    disresultlist(i)=dislist(i)
  #  disname=disname+'ODis(i) '
  #WriteTable disresultlist,'Distances of the nearest (distances) waters less than (dis) A from (sel)','(disname)'
  
  # Example: Save all ligand-receptor contacts with a distance smaller than 5 A in an extra tab file
  #contactdatalist()=ListConAtom (ligandsel) and Obj Solute,not (ligandsel) and Obj Solute,Cutoff=5,Results=3
  #contact='None'
  #if count contactdatalist
  #  contact=''
  #  for i=1 to count contactdatalist step 3
  #    for j=1 to 2
  #      atom(j)=ListAtom (contactdatalist(i+j-1)),Format='ATOMNUM.ATOMNAME.RESNAME3RESNUM.MOLNAME'
  #    contact=contact+'|(atom1),(atom2),(0.00+contactdatalist(i+2))'
  #WriteExtraTable 'HCON',contact,'At1,At2,D[A]'
 
  # Example: Scatter plot of dihedral angles Chi1 (X-axis) vs Chi2 (Y-axis) of Asp 22 in molecule A
  #resid='Asp 22 Mol A'
  #chi1=Dihedral N Res (resid),CA,CB,CG,Bound=Yes
  #chi2=Dihedral CA Res (resid),CB,CG,CD,Bound=Yes
  #ScatterPlot chi,1,'Dihedral angles Chi1 vs Chi2 of (resid)','Dihedral angle Chi1 in degrees','Dihedral angle Chi2 in degrees','Chi1 Chi2'
 
  # Example: Scatter plot of dihedral angles Chi1, Chi2 and Chi4 (Y-axis) vs Chi3 (X-axis) of Lys 22 in molecule B
  #resid='Lys 22 Mol B'
  #chi1=Dihedral N Res (resid),CA,CB,CG,Bound=Yes
  #chi2=Dihedral CA Res (resid),CB,CG,CD,Bound=Yes
  #chi3=Dihedral CB Res (resid),CG,CD,CE,Bound=Yes
  #chi4=Dihedral CG Res (resid),CD,CE,NZ,Bound=Yes
  #ScatterPlot chi,3,'Dihedral angles Chi1,Chi2,Chi4 vs Chi3 of (resid)','Dihedral angle Chi3 in degrees','Dihedral angle in degrees','Chi1 Chi2 Chi3 Chi4'
  
  # Example: Plot the distance between the atom with the highest and the one with the lowest Y-Coordinate
  #_,ly,_,_,y=GroupBox Obj Solute,Nuclear
  #atom1=ListAtom Obj Solute and LocalY<(y-ly/2+.001)
  #atom2=ListAtom Obj Solute and LocalY>(y+ly/2-.001)
  #dis=Distance (atom1),(atom2)
  #Plot dis,'Distance between atoms with highest/lowest Y','Distance in [A]','DisHiLoY'
  
  # Example: Plot the total MMPBSA energy of the system
  # IMPORTANT: This function must be the last in AnalyzeInsideCell,
  # as removing the Solvent turns off the running simulation.
  #surfcost=(0.65e0/6.02214199e20)*JToUnit
  #RemoveObj Solvent
  #epot=Energy
  #esolcoulomb,esolvdw=SolvEnergy PBS
  #molsurf=Surf molecular
  #esol=esolcoulomb+esolvdw+molsurf*surfcost
  #etot=epot+esol
  #AddObj Solvent
  #Plot etot,'Total MMPBSA energy of the system','MMPBSA Energy','MMPBSAenergy'
  
# Define the analyses to perform outside the simulation cell, without periodic boundaries
# =======================================================================================
def AnalyzeOutsideCell
  global fof
  
  # If the force field is YASARA2, plot the model quality
  if fof=='YASARA2'
    for checktype in 'dihedrals','packing1d','packing3d'
      zscore(checktype)=CheckObj Solute,(checktype)
    zscore=zscoredihedrals*0.145+zscorepacking1d*0.390+zscorepacking3d*0.465 
    Plot zscore,'Structure quality of the solute','Quality Z-score','Quality'
  
  # Plot the radius of gyration of the solute
  Plot "RadiusObj Solute,Center=Mass,Type=Gyration",'Radius of gyration of the solute','Radius in Angstrom','RadGyration'
  
  # The following examples provide a few hints for other things to analyze while the
  # simulation is *not* active (periodic boundaries removed), but the starting structure
  # has not yet been added to the soup to perform RMSD calculations,
  # just delete the '#' character to uncomment the command(s) of each example.
  
  # Example: Plot the distance between the carboxyl group of Glu 123 (Cdelta) in molecule/chain A
  #          and the guanidinium group of Arg 345 (Czeta) in molecule/chain B:
  # NOTE: This distance measurement does not consider periodic boundaries.
  #       Distances measured between separate molecules that cross periodic boundaries
  #       and are not connected with covalent bonds may be incorrect,
  #       but normally the boundary crossing is avoided by the CorrectDrift command.
  #       Alternatively, look at the Distance example further above in the 'AnalyzeInsideCell' section.
  #Plot "Distance CD Res Glu 123 Mol A, CZ Res Arg 345 Mol B",'Distance between Glu 123 [Cdelta] and Arg 345 [Czeta]','Distance','DGLUARG'
  
  # Example: Plot the angle between the transmembrane helix formed by residues 35-65
  # in Mol A and the membrane normal (which is parallel to the Y-axis, since the membrane is
  # parallel to the XZ-plane) using formula angle=acos(dotpro(dir,Y)/(len(dir)*len(Y))) 
  #_,_,_,dir()=GroupLine CA Protein Res 35-65 Mol A
  #Plot (acos dir2),'Angle between the transmembrane helix and the membrane normal','Angle','Phi'

  # Example: Plot the angle between the secondary structure elements formed by
  # residues 106-140 and 149-169 in Mol B, in the range -180 to 180 (see GroupAngle docs)
  #Plot "GroupAngle CA Protein Res 106-140 Mol B,CA Protein Res 149-169 Mol B,Range=360",'Angle between the secondary structure elements','Angle','Psi'
  
  # Example: Plot the distance between two centers of mass, e.g. the loop from
  #          residue Ala 205 to Glu 210, and the ligand NAD: 
  #cenA()=GroupCenter Res Ala 205 - Glu 210 Obj Solute, Type=Mass
  #cenB()=GroupCenter Res NAD Obj Solute, Type=Mass
  #Plot norm (cenA-cenB),'Distance between two centers of mass','Distance','DAB'
  
  # Example: Plot the Coulomb and VdW binding energies between the solute and the solvent
  # EBind=ESolute + ESolvent - ETotal. For ligand binding look at md_analyzebindenergy macro.
  #etotvdw,etotcoulomb=Energy VdW,Coulomb
  #RemoveObj Solvent
  #esoluvdw,esolucoulomb=Energy VdW,Coulomb
  #AddObj Solvent
  #RemoveObj Solute
  #esolvvdw,esolvcoulomb=Energy VdW,Coulomb
  #AddObj Solute
  #ebindlist()=(esoluvdw+esolvvdw-etotvdw),(esolucoulomb+esolvcoulomb-etotcoulomb)
  #Plot ebindlist,'Coulomb and VdW binding energies between the solute and the solvent','Binding Energy','CoulombEb VdWEb'

# Define the analyses to perform with respect to the starting structure, also outside the simulation cell
# =======================================================================================================
# The reference structure (object SoluteRef) has been added to the soup now
def AnalyzeChange
  global casel,caselref,ligandsel,ligandselref
  
  # Plot Calpha (or C1* of nucleic acids), backbone and all-atom RMSDs
  # Trick: If the solute contains neither CA nor backbone atoms,
  # we simply assign the SupAtom error code (=0)
  resultlist1=SupAtom (casel),(caselref)
  resultlist2=SupAtom Backbone Obj Solute,Backbone Obj SoluteRef
  resultlist3=SupAtom Element !H Obj Solute,Element !H Obj SoluteRef
  Plot resultlist,'Solute RMSD from the starting structure','RMSD in Angstrom','RMSDCa RMSDBb RMSDAll'
  
  # Plot ligand movement RMSD after superposing on the receptor
  if ligandsel!=''
    if casel!='None'
      SupAtom (casel) and not (ligandsel),(caselref) and not (ligandselref)
    else
      SupAtom Element !H Obj Solute and not (ligandsel),Element !H Obj SoluteRef and not (ligandselref)
    result=RMSDAtom (ligandsel),(ligandselref)
    Plot result,'Ligand movement RMSD after superposing on the receptor','RMSD in Angstrom','RMSDLigMove'
  
  # Plot ligand conformation RMSD after superposing on the ligand
  if ligandsel!=''
    result=SupAtom (ligandsel),(ligandselref)
    Plot result,'Ligand conformation RMSD after superposing on the ligand','RMSD in Angstrom','RMSDLigConf'
  
  # The following examples provide a few hints for other things to analyze while the
  # simulation is *not* active (i.e. things that can't be analyzed with periodic boundaries)
  # with respect to the starting structure, just delete the '#' character to uncomment
  # the command(s) of each example.

  # Example: Plot the backbone RMSD for residues 1-120:
  #Plot "SupAtom Backbone Res 1-120 Obj Solute,Backbone Res 1-120 Obj SoluteRef",'Backbone RMSD','RMSD','BBRMSD'

  # Example: Plot the sidechain heavy-atom RMSD:
  #Plot "SupAtom Sidechain Element !H Obj Solute,Sidechain Element !H Obj SoluteRef",'Sidechain heavy-atom RMSD','RMSD','SCHARMSD'

  # Example: Plot how far residue HP6 86 moved since the start (after superposing on Calphas) 
  #SupAtom (casel),(caselref)
  #Plot "GroupDistance Res HP6 86 Obj Solute,Res HP6 86 Obj SoluteRef",'Drift HP6 86','Distance','DHP6'

# Optional text printed in the report in front of each analysis
# =============================================================
#
# Please see the user manual at Recipes > Run a molecular dynamics simulation > Analyzing a trajectory
# (scroll to the end) for detailed instructions how to add your own texts.

# Get the force field
fof=ForceField

# AnalyzeInsideCell:

# Text paragraph for the 'Simulation cell lengths' plot
CellLengthXText='Conformational changes of the simulated solute molecules lead to fluctuations in density. '
                'If the simulation box has a constant size, changes in density lead to changes in pressure. '
                'This is not realistic, because molecules normally "live" in a constant pressure environment. '
                'During the simulation the cell is therefore rescaled to maintain a constant cell '
                'pressure. Depending on the chosen pressure control mode, '
                'the three cell axes are either rescaled together [Manometer1D], partly together '
                '[X- and Z-axes, Manometer2D, used for membrane simulations], independently [Manometer3D], '
                'or not at all [Off]. You can deduce the pressure control mode from the plot below.'

# Text paragraph for the 'Total potential energy of the system' plot
TotalEnergyText='The total potential energy of the system is plotted, '
                'according to the (fof) force field. If you ran the simulation '
                'with a different force field, you need to adapt the `ForceField` command '
                'at the top of this macro accordingly.\n\n'
                'When the simulation is started from an energy-minimized "frozen" conformation, '
                'there is usually a sharp increase in energy during the first picoseconds, '
                'since the added kinetic energy is partly stored as potential energy. '
                'Also on a larger time-scale, the potential energy will often not decrease. '
                'A common reason are counter ions. These are initially placed at the '
                'positions with the lowest potential energy, usually close to charged solute groups, '
                'from where they detach to gain entropy, but also potential energy. '

# Text function for the 'Potential energy components' plot
def PrintBondInfo
  global fof
  text='The following individual components of the total potential energy are plotted: '
       'bond energies [Bond], bond angle energies [Angle], dihedral angle energies [Dihedral], '
       'planarity or improper dihedral energies [Planarity], Van der Waals energies [VdW]'
  if fof=='YASARA2'
    text=text+', electrostatic energies [Coulomb], 1D packing energies [Packing1D] and 3D packing energies [Packing3D]. '
  else
    text=text+' and electrostatic energies [Coulomb]. '
  text=text+'Force field energies help to judge the structural quality of a protein: '
            'distortions of local covalent geometry can be found by looking at the bond, angle and planarity energies. '
            'Unrealistically close contacts [bumps] lead to a high Van der Waals energy, '
            'just like a large number of hydrogen bonds [since they pull the atoms closer than '
            'their normal Van der Waals contact distance]. The Coulomb energy is the least '
            'informative, because it strongly depends on the amino acid composition '
            '[e.g. proteins with a net charge have a higher Coulomb energy].'
  WriteReport Paragraph,'(text)'

# Text paragraph for the 'Protein secondary structure' plot
HelixText='The total percentages of alpha helices, beta sheets, turns, coils, 3-10 helices and pi helices are '
          'calculated and plotted. For clarification, a turn is simply a stretch of four residues that are not '
          'part of other secondary structure elements and form a hydrogen bond between the O of the first and '
          'the NH of the last residue. A coil is anything that does not fit into the other categories. '
          'Note that pi-helices [helices with hydrogen bonds between residues N and N+5] are rather unstable and '
          'thus do not normally occur in proteins, except for short bulges in alpha helices '
          '[which are often the result of single residue insertions and prolines]. '

# Text paragraph for the 'Per-residue protein secondary structure' plot
SecStrText='The following plots show the protein secondary structure per residue as a function of '
           'simulation time. They are helpful to monitor protein folding and all other kinds of structural '
           'changes. The default secondary structure colors are used, you can change them at '
           'View > Color > Parameters > Secondary structure colors. One plot per protein molecule is shown.'

# Text paragraph for the 'Per-residue number of contacts' plot
ConText='The number of contacts per residue as a function of simulation time is shown in the following plots. '
        'There is one plot for each protein or nucleic acid molecule. Even though contacts between atoms separated '
        'by up to four chemical bonds are excluded, neighboring residues in the molecule usually have enough close '
        'atoms to be counted as a contact. Consequently residues with zero contacts are very rare and often glycines. '
        'The number of contacts tells you how densely a certain residue range is packed and allows to identify structurally '
        'very important residues, e.g. a phenylalanine in the hydrophobic core can contact 15 or more other residues. '

# Text paragraph for the 'Per-residue ligand interactions of the receptor' plot
LigResIntText='The following plots show the types of interactions by receptor residues with the ligand as a function of simulation time. '
            'There is one plot for each protein or nucleic acid molecule. Three types of interactions are shown: Hydrogen bonds [red] ',
            'hydrophobic [green] and ionic interactions [blue]. Also mixtures of these three colors can show up if a '
            'certain residue is involved in more than one type of interaction with the ligand [see plot legend]. '

# Text function for the 'Per-atom receptor interactions of the ligand' plot
def PrintLigAtomIntInfo
  global resultbase,block
  WriteReport Paragraph,
   'The following plots show the types of interactions by atoms of the ligand with the receptor as a function of simulation time. '
   'There is one plot for each ligand molecule. Three types of interactions are shown: Hydrogen bonds [red] '
   'hydrophobic [green] and ionic interactions [blue]. Also mixtures of these three colors can show up if a '
   'certain ligand atom is involved in more than one type of interaction [see plot legend]. '
   'Atom numbers in the plots match the atom numbers shown in the screenshot of the ligand in the figure below. '
  delimage='Yes'
  if count block
    # We want to show the ligand with atom numbers next to the plot for every block for convenience
    delimage='No'
  WriteReport Image,Filename=(resultbase)_ligandatomnum.png,Style=Figure,
    Caption='A ray-traced picture of the ligand with atoms labeled by their respective atom number',Delete=(delimage)

# Text paragraph for the 'Surface areas of the solute' plot
SurfVdWText='The Van der Waals [SurfVdW], molecular [SurfMol] and solvent accessible [SurfAcc] surface areas of the '
            'solute in A^2 are plotted. The difference between these surface types can be summarized as follows:\n\n'
            '`Van der Waals surface`: if you think of atoms as spheres with a given Van der Waals radius, '
            'then the Van der Waals surface consists of all the points on these spheres that are not inside another sphere. '
            'In practice, the Van der Waals surface is of limited use, because it can be found throughout a protein and '
            'does not tell much about the interaction with the solvent.\n\n'
            '`Molecular surface`: this is the Van der Waals surface from the viewpoint of a solvent molecule, '
            'which is a much more useful concept. The water is assumed to be a sphere of a given radius '
            '[also called the water probe], that rolls over the solute. '
            'Those parts of the Van der Waals surface that the water probe can touch are simply copied to the molecular surface '
            '[and called the contact surface]. Clefts in the Van der Waals surface that are too narrow for the water probe to enter '
            'are replaced by the Van der Waals surface of the water probe itself [and called the reentrant surface]. '
            'So the molecular surface is a smooth composition of two Van der Waals surfaces: '
            'the one of the solute and the one of the solvent molecule while it traces the contours of the solute. '
            'Other common names for the molecular surface are solvent excluded surface or Connolly surface.\n\n'
            '`Solvent accessible surface`: this surface consists of all the points that the center of the water probe '
            '[i.e. the nucleus of the oxygen atom in the water molecule] can reach while rolling over the solute. '
            'The shortest possible distance between the water oxygen nucleus and a solute atom is simply '
            'the sum of the Van der Waals radii of the solute atom and the water probe.'

# Text function for the 'Number of hydrogen bonds in solute' plot
def PrintSoluteHBondsInfo
  WriteReport Paragraph,
    'The number of hydrogen bonds inside the solute is plotted below. '
    'One hydrogen bond per hydrogen atom is assigned at most, '
    'picking the better one if two acceptors are available.'
    'The following formula yields the bond energy in [kJ/mol] '
    'as a function of the Hydrogen-Acceptor distance in [A] and two scaling factors:'
  WriteReport Image,Filename=(YASARADir)/doc/ListHBoAtomResMolObj1.png,Style=Center,Name="formula_energyhbo0"
  WriteReport Paragraph,
    'The first scaling factor depends on the angle formed by Donor-Hydrogen-Acceptor:'
  WriteReport Image,Filename=(YASARADir)/doc/ListHBoAtomResMolObj2.png,Style=Center,Name="formula_energyhbo1"
  WriteReport Paragraph,
    'The second scaling factor is derived from the angle formed by Hydrogen-Acceptor-X, '
    'where X is the atom covalently bound to the acceptor. If X is a heavy atom:'
  WriteReport Image,Filename=(YASARADir)/doc/ListHBoAtomResMolObj3.png,Style=Center,Name="formula_energyhbo2"
  WriteReport Paragraph,
    'If X is a hydrogen, slightly smaller angles are allowed:'
  WriteReport Image,Filename=(YASARADir)/doc/ListHBoAtomResMolObj4.png,Style=Center,Name="formula_energyhbo3"
  WriteReport Paragraph,
    'A hydrogen bond is counted if the hydrogen bond energy obtained with this formula '
    'is better than 6.25 kJ/mol [or 1.5 kcal/mol], which is 25% of the optimum value 25 kJ/mol. '

# Text paragraph for the 'Number of hydrogen bonds between solute and solvent' plot
SltSlvHBondsText='The plot shows the number of hydrogen bonds between solute and solvent. '
                 'Together with the plot above, it is a good indicator for successful protein folding, '
                 'indicated by a decreasing number of bonds with the solvent and a growing '
                 'number of bonds within the solute.'

# Text paragraph for the 'Number of hydrogen bonds made by (hbosel)' plot
SelHBondsText='The plot shows three kinds of hydrogen bonds made by atoms selected with `hbosel`: '
              'Internal hydrogen bonds [blue], hydrogen bonds with the rest of the solute [green], '
              'those with the water [red] including membrane molecules if any, and the sum of all three [gray]. '
              'Internal bonds will be counted twice, once in each direction, so that the total number of bonds '
              'stays the same if alternative hydrogen bonds are formed. '
              'The green graph is an indicator for the completeness of a dissociation or docking event, '
              'assuming that there are hydrogen bonds involved at all. '
              'It is noteworthy that hydrogen bonds do not contribute significantly to the free energy of binding, '
              'as long as receptor and ligand can form alternative hydrogen bonds internally or with the solvent. '
              'Only if unsatisfied hydrogen bond donors or acceptors get '
              'buried, i.e. the sum of all hydrogen bonds shown in gray gets smaller, '
              'there is an adverse impact on the binding energy.'

# Text function for the 'Hydrogen bonds made by (hbosel)' table
def PrintHB1Atm1Info
  global hbosel,hbondsmax,EnergyUnit
  # The 'caller3' keyword allows to access variables of the 3rd last calling function
  caller3 acceptors,donors,bonds
  WriteReport Paragraph,
    'The following table shows all hydrogen bonds made by (hbosel), excluding water molecules. '
    'With (acceptors) acceptors and (donors) donors, '
    'a total number of (bonds) hydrogen bonds are possible. '
    '(hbondsmax) hydrogen bonds are listed - labeled HB1 to HB(hbondsmax). '
    'The donor of the bonding pair is labeled as Atm1 and the acceptor as Atm2, respectively. '
    'Internal bonds of (hbosel) will be listed twice, once in each direction. '
    'The atom ID separates atom name, residue ID and molecule name with dots. A lower-case "h" '
    'indicates hetgroups. E and D are short for the hydrogen bonding energy in [(EnergyUnit)] '
    'and the distance between the bonding partners in [A]. '
    'To list other hydrogen bonds, edit the `hbosel` variable at the beginning of this macro. '
    'To list more or fewer hydrogen bonds, edit the `hbondsmax` variable at the beginning of this macro.'

# AnalyzeOutsideCell:

# Text function for the 'Structure quality of the solute' plot
def PrintQualityInfo
  WriteReport Paragraph,
    'When validating a structure, one can use a gold standard of highest resolution reference '
    'structures to obtain estimates for the expected average energy and its standard deviation, '
    'and then calculate how many standard deviations the actual energy is away from the average, '
    'thereby obtaining a Z-score:'
  WriteReport Image,Filename=(YASARADir)/doc/CheckAtomResObjAll1.png,Style=Center,Name="formula_zscore"
  WriteReport Paragraph,
    'In the formula above, `x` is the energy of the current structure, '
    '`my` and `sigma` are the average value and the standard deviation '
    'of the energies in the gold standard population. '
    'Assuming that the gold standard population has a standard distribution, '
    '95.4% of all proteins have Z-scores between -2 to +2. '
    'The others can be called outliers. '
    'Positive outliers are usually small perfect proteins like a single alpha helix, '
    'and negative outliers are proteins with serious errors.'
  MakeTab zscores,2,2
  descriptionlist()='disgusting','terrible','bad','poor','satisfactory','good'
  for i=-5 to 0
    Tabulate '< (i)','(descriptionlist(i+6))'
  Tabulate '> 0','optimal'
  WriteReport Table,zscores,"Mapping between Z-score and human language.",.0f,InfoColumnName='Z-score',DataColumnName='Description'

# Text function for the 'Radius of gyration of the solute' plot
def PrintRadGyrationInfo()
  WriteReport Paragraph,
    'After determining the center of mass of the solute, the radius '
    'of gyration is calculated and plotted according to this formula:'
  WriteReport Image,Filename=(YASARADir)/doc/RadiusAtomResMolObjAll2.png,Style=Center,Name="formula_gyrrad"
  WriteReport Paragraph,
    'In this formula, `C` is the center of mass, and `Ri` is the position of atom `i` of `N`.'

# AnalyzeChange:

# Text function for the 'Solute RMSD from the starting structure' plot
def PrintRMSDCaInfo
  global casel,calphas
  WriteReport Paragraph,
    'The plot shows Calpha [RMSDCa], backbone [RMSDBb] and all-heavy atom [RMSDAll] RMSDs calculated '
    'according to this formula, where `Ri` is the vector linking the positions of atom `i` [of `N` atoms] '
    'in the reference snapshot and the current snapshot after optimal superposition: '
  WriteReport Image,Filename=(YASARADir)/doc/RMSDAtomResMolObj1.png,Style=Center,Name="formula_rmsd"
  if casel=='None'
    text='Less than three atoms matched the Calpha selection `(casel)`, therefore the Calpha RMSD plot '
         'graph is set to flat zero because at least three atoms are needed for structure superposition. '
  else
    text='The selection for the Calpha RMSD calculation is `(casel)`, this matched (calphas) atoms. '
    if casel=='CA Protein or C1* NucAcid and Obj Solute'
      text=text+'The Calpha selection thus includes the main backbone carbon C1* of nucleic acids, so '
                'the plot also shows a Calpha RMSD if you simulate just nucleic acids. In simulations '
                'of protein-DNA complexes, the Calpha RMSD therefore considers the DNA too. '
  text=text+'To change the Calpha selection, edit the `casel` variable at the beginning of this macro.'
  WriteReport Paragraph,'(text)'

# Text paragraph for the 'Ligand movement RMSD after superposing on the receptor' plot
RMSDLigMoveText='The following plot shows the RMSD of the ligand heavy atoms '
                'over time, measured after superposing the receptor on its reference structure. '
                'This procedure delivers information about the movement of the ligand in its binding pocket.'

# Text paragraph for the 'Ligand conformation RMSD after superposing on the ligand' plot
RMSDLigConfText='This plot displays the RMSD of the ligand atoms '
                'over time, measured after superposing on the reference structure of the ligand. '
                'The gained data summarize the conformational changes of the ligand. '

# Normally no change required below this point
# ============================================

# Do we have a target?
if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target structure"

# Check order of residue min/max user selection
if resnummax<resnummin
   RaiseError 'The maximum residue number (resnummax) is smaller than the minimum (resnummin), please swap resnummin and resnummax'

Clear
Console Off
SurfPar Molecular=Numeric

if not count block
  id=''

# Load scene with water or other solvent
waterscene=FileSize (MacroTarget)_water.sce
solventscene=FileSize (MacroTarget)_solvent.sce
if waterscene
  LoadSce (MacroTarget)_water
elif solventscene 
  LoadSce (MacroTarget)_solvent
else
  RaiseError 'Could not find initial scene file (MacroTarget)_water.sce. You must run a simulation with the macro md_run first'

UnselectAll

# Choose object names to make user-defined analyses easier
objs1=CountObj Solute and not (soluteobj)
objs2=CountObj SoluteRef
if sum objs
  RaiseError "The object names 'Solute' and 'SoluteRef' are reserved for use by this macro, please rename your objects and try again"
solutename = NameObj (soluteobj)
NameObj Water,Solvent
objs=CountObj Solvent
slvfirstatm=0
if !objs and Objects==2
  waters=CountRes Water
  if waters>9
    # Number of waters is larger than arbitrary cutoff
    _,last=SpanAtom all with bond to all and not Water
    _,last=SpanAtom all with arrow to all or 1-(last)
    if last<atoms
      # Solvent is part of the solute, split it off
      slvfirstatm=last+1
      SplitAtom (slvfirstatm)
      SplitObj (soluteobj),Yes,Atom (slvfirstatm)
      NameObj Atom (slvfirstatm),Solvent
      slvdetected=1
NameObj (soluteobj),Solute

# Avoid the old empty molecule names
NameMol ' ','_' 

# Verify Calpha selection
if casel==''
  calphas=0
else
  calphas=CountAtom (casel)
if calphas<3
  # We cannot superpose 1 or 2 Calpha atoms
  casel='None'

ShowMessage "Preparing analysis, please wait..."
Wait 1

# Backwards compatibility: Starting with YASARA version 12.8.1, XTC trajectories no longer contain a number in the filename
old=FileSize (MacroTarget)00000.xtc
if old
  RenameFile (MacroTarget)00000.xtc,(MacroTarget).xtc

# Determine trajectory format
for format in 'xtc','mdcrd','sim'
  found=FileSize (MacroTarget).(format)
  if found
    break
  
# Identify the ligand: choose the largest hetgroup with >6 atoms if any
# Changes here also in md_runsteered.
LigandObjWarning=''
if ligandsel==''
  mols=CountMol Obj Solute
  if mols>1
    ligandreportstring='automatically by YASARA'
    reslist()=ListRes Hetgroup !Water Obj Solute with >0 bonds to all
    reslenmax=6
    ligandname=''
    carbohydlist=()
    for res in reslist
      reslen=CountAtom (res)
      resname=NameRes (res)
      carbons=CountAtom Element C (res) with bond to Element O
      if carbons>4 and resname not in carbohydlist
        carbohydlist(1+count carbohydlist)=resname
      if reslen>reslenmax
        reslenmax=reslen
        ligandname=resname
    if ligandname!=''
      # Build ligand selection, adding quotes to deal with unusual residue names
      if ligandname in carbohydlist
        # Ligand is a carbohydrate, treat all carbohydrate residues as ligands
        ligandsel='Res '+join carbohydlist
      else
        ligandsel='Res (ligandname)'
else
  ligandobjlist()=ListObj (ligandsel) and not Membrane Solvent
  if !count ligandobjlist
    RaiseError 'Your ligand selection "(ligandsel)" did not match any atoms. If you want to select a residue add "Res" at the beginning of your selection.'
  if count ligandobjlist>1
    RaiseError 'Your ligand "(ligandsel)" is present in objects (ligandobjlist), but only one ligand object is allowed'
  solfirstatm,sollastatm=SpanAtom Obj Solute
  ligfirstatm,liglastatm=SpanAtom Obj (ligandobjlist1)
  if ligfirstatm==sollastatm+1
    # Ligand is not part of the solute object, but atoms numbers are adjacent so we can simply join them
    ligatoms=CountAtom (ligandsel)
    ligandobjname=NameObj (ligandobjlist1)
    # Set segment name for the ligand
    SegAtom (ligandsel),LigO
    # We will reference the ligand by segment name
    ligandsel='Segment LigO'
    JoinObj (ligandobjlist),Solute
    if ligatoms!=liglastatm-ligfirstatm+1
      # Ligand is only part of the object, show a warning
      LigandObjWarning='`WARNING:` The selected ligand is only part of a larger object. The remaining atoms in '
                        'object (ligandobjlist1) with name `(ligandobjname)` were treated as solute atoms.'
  elif ligfirstatm!=solfirstatm 
    RaiseError 'Ligand (ligandsel) is not part of the solute, no ligand RMSD values can be calculated. '
               'You could click Edit > Join > Objects to join the ligand object (ligandobjlist1) with the solute object (soluteobj) '
               'in the original scene, but this would change the order of atoms and require to run the MD once again.'
  ligandreportstring='by the user'

if refsnapshot=='average'
  # We superpose on the time average structure
  filename='(resultbase)_average(id).pdb'
  exists=FileSize (filename)
  if !exists
    RaiseError "No time-average structure has been calculated yet, cannot superpose onto it. Run the macro once with refsnapshot=0 (or any other number), then run again with refsnapshot='average'"
  refobj=LoadPDB (filename)
  SupObj (refobj),Solute
else
  if refsnapshot
    # We superpose on a certain snapshot
    if format=='sim'
      LoadSim (MacroTarget)(00000+refsnapshot)
    else
      Load(format) (MacroTarget),(refsnapshot+1)
  # Duplicate the intial object for RMSD calculation
  refobj=DuplicateObj Solute
# Rename reference object for convenience
NameObj (refobj),SoluteRef
# Select atoms in the reference structure.
# We use ligandref to allow for joining of the ligand and soluteobj,
# and to be able to use atom selections.
ligandselref=''
caselref='None'
if ligandsel!=''
  ligandselref=ligandsel+' Obj SoluteRef'
  ligandsel=ligandsel+' Obj Solute'
if casel!='None'
  caselref=casel+' and Obj SoluteRef'
  casel=casel+' and Obj Solute'
# Remove reference structure
RemoveObj SoluteRef
# Remember cell size of the reference snapshot for e.g. ESP
cellsize1,cellsize2,cellsize3=Cell

if rmsdmin
  # Cluster the trajectory, get a list of representative structures
  clusterobjlist()=refobj

# Check DCCM selection
if dccmsel!=''
  dccmunits=Count(dccmsel) and Obj Solute
  if !dccmunits and dccmsel!='Atom CA Protein or C1* NucAcid'
    RaiseError 'The DCCM selection (dccmsel) did not match any atoms'

# Create the headers for saving the final tables and make the ray-traced pictures for the report.
# Do this already now to uncover errors in user-provided graph names
# and to determine the total number of table columns
# Calling all analysis functions with task='AddHeader' lets them just join their graph names to the headers
task='AddHeader'
tables=0
plotmaptabs=0

# First analyze inside the cell during a simulation
Sim On
AnalyzeInsideCell
# Then outside the cell
Sim Off
AnalyzeOutsideCell
# Then with respect to the reference structure
AddObj SoluteRef
AnalyzeChange
RemoveObj SoluteRef
if !count introwritten
  # Get ray-traced picture of the entire system
  ShowSystem
  SaveScreenshot 'sim','simulated system'
  # Get ray-traced picture of the solute
  ShowSolute
  Style Ribbon,Stick
  StickRadius 20
  SaveScreenshot 'solute','solute object'
  if ligandsel!=''
    # Get ray-traced picture of the ligand
    ligatomlist()=ListAtom (ligandsel)
    ligandobj=DuplicateRes (ligandsel)
    dupatomlist()=ListAtom Obj (ligandobj)
    SwitchAll off
    SwitchObj (ligandobj),on
    NiceOriObj (ligandobj)
    LabelAtom Obj (ligandobj) Element !H and Element !N,Format=ATOMNAME,Height=.3,Color=Black
    LabelAtom Obj (ligandobj) Element N,Format=ATOMNAME,Height=.3,Color=White
    LabelAtom Obj (ligandobj) Element H,Format=ATOMNAME,Height=.2,Color=Black
    ZoomObj (ligandobj),Steps=0
    Style BallStick
    SaveScreenshot 'ligand','ligand'
    UnlabelAtom Obj (ligandobj)
    for i=1 to count ligatomlist
      LabelAtom (dupatomlist(i)) and Element N,(ligatomlist(i)),Height=.2,Color=White
      LabelAtom (dupatomlist(i)) and Element !N,(ligatomlist(i)),Height=.2,Color=Black
    SaveScreenshot 'ligandatomnum','ligand with atom numbers'
    DelObj (ligandobj)
  # And show the entire system again
  ShowSystem

# Create the main tables with the now known number of columns
for i=1 to tables
  MakeTab (tabnamelist(i)),Dimensions=2,Columns=(tabcolslist(i))

# The data of the main table for the Y-axis can be found in columns 3+, column 1/2 contains the simulation time in [ps]/[ns]
ycolumn=3

# Set tabrowsmax
if tabrowsmax=='all'
  tabrowsmax=0

# Run the actual analysis
task='Tabulate'
i=00000+firstsnapshot
emin=1e99
last=0
while !last and snapshots!=0
  # Load next snapshot from SIM or XTC trajectory
  lastsnapshot=i
  if format=='sim'
    # Set last to 1 if last snapshot loaded
    sim=FileSize (MacroTarget)(i+1).sim
    if not sim
      last=1
    LoadSim (MacroTarget)(i)
  else
    last=Load(format) (MacroTarget),(i+1)
  Sim Pause
  if central
    # Keep a chosen atom at the center of the cell (Cell returns center as values 7-9) 
    _,_,_,_,_,_,cen()=Cell
    pos()=PosAtom (central)
    MoveAtom all,(cen-pos)
  # Add time in picoseconds and nanoseconds to table
  simtime=Time
  ShowMessage 'Analyzing snapshot (0+i) at (0+(simtime/1000)) ps'
  Wait 1
  for j=1 to tables
    SelectTab (tabnamelist(j))
    Tabulate (simtime/1000),(simtime/1000000)
  # Perform analysis inside cell
  AnalyzeInsideCell
  if count rdfsellist==4
    # Collect data for radial distribution function
    BinDistance (rdfsellist)
  # Prepare to save the minimum energy structure (ignoring solvent-solvent interactions)
  e=EnergyObj Solute
  # Stop simulation
  Sim Off
  if e<emin
    # Save minimum energy structure
    emin=e
    SaveSce (resultbase)_energymin(id)
    SavePDB Solute,(resultbase)_energymin(id)
  if last
    # Save last structure
    SaveSce (resultbase)_last
    SavePDB Solute,(resultbase)_last
  if pdbsaved
    # Save an entire trajectory of solute PDB files
    SavePDB Solute,(resultbase)_(i)
  # Perform analysis outside cell
  AnalyzeOutsideCell
  # Add reference structure to the soup to perform RMSD calculations etc.
  AddObj SoluteRef
  if espmethod!=''
    # Superpose solute with the reference structure
    TransferObj Solute Membrane,SoluteRef,Fix
    SupAtom Element !H Obj Solute,Element !H Obj SoluteRef
    # Remove reference structure and solvent, we use an implicit solvent model for PBS and vacuum for PME
    RemoveObj SoluteRef Solvent
    # Transform the membrane the same way
    move()=Transformation Translation
    rot()=Transformation Angles
    RotateAtom Obj Membrane,(rot)
    MoveAtom Obj Membrane,(move)
    # Set cell size to the one of the reference snapshot
    Cell (cellsize)
    # Sim On to wrap around atoms at the periodic boundary
    Sim On
    # Save electrostatic potential
    SaveESP (resultbase)(i).cube.gz,(espmethod),Gaussian
    # Load electrostatic potential to calculate the average
    tbl=LoadTab (resultbase)(i).cube.gz
    vallist()=Tab (tbl)
    if i==firstsnapshot
      esptbl=tbl
      espvallist=vallist
    else
      espvallist=espvallist+vallist
      DelTab (tbl)
    AddObj SoluteRef Solvent
  # Perform analysis outside cell, with reference structure present in the soup
  AnalyzeChange
  # Superpose again and add the current atom positions to internal table to obtain RMSF and average positions
  SupAtom Element !H Obj Solute,Element !H Obj SoluteRef
  AddPosAtom Obj Solute
  if rmsdmin
    # Add all representative cluster members to the soup, refobj is already present
    AddObj (join clusterobjlist) and not SoluteRef
    # Superpose the current structure on all representative cluster members
    for j=1 to count clusterobjlist
      aarmsd = SupAtom Element !H Obj Solute,Element !H Obj (clusterobjlist(j))
      if aarmsd<=rmsdmin
        break
    if aarmsd>rmsdmin
      # Found a new cluster member, add to list
      obj = DuplicateObj Solute
      members=count clusterobjlist+1
      clusterobjlist(members)=obj
      NameObj (obj),Cluster(members)
      # Store the snapshot number
      PropObj (obj),(i)
    # Remove cluster members from soup again
    RemoveObj (join clusterobjlist)
  else
    # No cluster analysis, remove initial structure again
    RemoveObj SoluteRef
  # Save last structure if last snapshot not already reached
  if !last
    for j=2 to snapshotstep
      # Set last (last file) to 1 if last snapshot loaded
      if format=='sim'
        sim=FileSize (MacroTarget)(i+j).sim
        if not sim
          LoadSim (MacroTarget)(i+j-1)
          last=1
      else
        # Set last (end of file) to 1 if last snapshot loaded
        last=Load(format) (MacroTarget),(i+j)
      if last
        # Save last structure
        Sim Off
        SaveSce (resultbase)_last
        SavePDB Solute,(resultbase)_last
        break
  # Next snapshot
  i=i+snapshotstep
  if snapshots!='all'
    snapshots=snapshots-1
if i==firstsnapshot
  RaiseError "This macro is meant to analyze a molecular dynamics trajectory created with md_run, but none was found in this directory"
snapshots=(-firstsnapshot+i)/snapshotstep

if rmsdmin
  # Save all cluster members
  if count block and block>1
    # Omit the refsnapshot except for the first block
    DelObj (clusterobjlist1)
    DelVar clusterobjlist1
  clustermembers=count clusterobjlist
  AddObj (join clusterobjlist)
  for j=1 to clustermembers
    # Get the snapshot number
    obj=clusterobjlist(j)
    snum=PropObj (obj)
    clusterfilenamelist(j)='(resultbase)_cluster(id)_(zeroed clustermembers+j)_snapshot(0+snum)'
    SaveYOb (obj),(clusterfilenamelist(j))
  DelObj (join clusterobjlist) and not SoluteRef
  RemoveObj SoluteRef

# Decide if the timescale shown in plots and tables should be picoseconds or nanoseconds
simtime=Time
if simtime<1000000
  xcolumn=1
  plottimestring='picoseconds'
  tabtimestring='Time [ps]'
  simtime=0.00+simtime/1000
else
  xcolumn=2
  plottimestring='nanoseconds'
  tabtimestring='Time [ns]'
  simtime=0.00+simtime/1000000

# Calculate time-average structure
AveragePosAtom Obj Solute
# Set B-factors, the dummy assignment '_ =' ensures that B-factors are not printed
_ = RMSFAtom Obj Solute,Unit=BFactor
if bfactorscale!=1.0
  # Scale B-factors so that they fit into the PDB format
  firstatm,lastatm=SpanAtom Obj Solute
  for i=firstatm to lastatm
    bf=BFactorAtom (i)
    BFactorAtom (i),(bf*bfactorscale)    
if refsnapshot!='average'
  # The time average structure has incorrect covalent geometry and should be energy minimized
  SavePDB Solute,(resultbase)_average(id)
# Additionally create an RMSF tables, in case B-factors are too large for the PDB format
MakeTab RMSF
firstatm,lastatm=SpanAtom Obj Solute
rmsflist()=RMSFAtom Obj Solute
for i=firstatm to lastatm
  res=ListAtom (i),Format="'ATOMNAME','RESNAME','RESNUM','MOLNAME'"
  Tabulate '(i)',(res),(rmsflist(i-firstatm+1))
SaveTab RMSF,(resultbase)_rmsf(id),Format=Text,Columns=6,NumFormat=8.2f,"Table of atomic Root Mean Square Fluctuations in [A]"
MakeTab RESRMSF
reslist()=ListRes Obj Solute,Format="'RESNAME','RESNUM','MOLNAME'"
resrmsflist()=RMSFRes Obj Solute
for i=1 to count reslist
  Tabulate (reslist(i)),(resrmsflist(i))
SaveTab RESRMSF,(resultbase)_rmsfres(id),Format=Text,Columns=4,NumFormat=8.2f,"Table of residue Root Mean Square Fluctuations in [A]" 

bound=Boundary
ycolumn=3
if !count introwritten
  introwritten=1
  # Make an info tab
  MakeTab Info,2,2
  # Tabulate info about proteins and nucleic acids
  infolist()='Protein molecules','Mol Protein','Protein residues','Res Protein','Protein atoms','Atom Protein',
             'Nucleic acid molecules','Mol NucAcid','Nucleic acid residues','Res NucAcid','Nucleic acid atoms','Atom NucAcid'
  for i=1 to count infolist step 2
    Tabulate '(infolist(i))'
    Tabulate Count(infolist(i+1))
  # Get list of unique residue names that are not part of the proteins,nucleic acids or the solvent
  resnamelist()=NameRes !Protein and !NucAcid and !Water
  if count resnamelist
    unilist=()
    for resname in resnamelist
      if resname not in unilist
        unilist(count unilist+1)=resname
    # Remove residues with only one atom
    for resname in unilist
      # Get number of atoms
      reslist()=ListRes (resname)
      residlist=()
      residueslist=()
      resatomslist=()
      for i=1 to count reslist
        resatoms=CountAtom Res (reslist(i))
        if resatoms not in resatomslist
          resatomslist(count resatomslist+1)=resatoms
          residueslist(count residueslist+1)=1
          residlist(count residlist+1)=reslist(i)
        else
          for j=1 to count resatomslist
            if resatoms==resatomslist(j)
              residueslist(j)=residueslist(j)+1
      for i=1 to count residueslist
        if resatomslist(i)==1
          # Single atom residue, get element
          element=ListAtom (residlist(i)),Format='ATOMElement'
          Tabulate 'Residue (resname) with 1 atom of element (element)',(residueslist(i))
        else 
          # Tabulate actual residues with number of atoms
          Tabulate 'Residue (resname) with (resatomslist(i)) atoms',(residueslist(i))
  # Tabulate water residues
  Tabulate 'Water residues'
  Tabulate CountRes Water
  # Tabulate total number of atoms
  Tabulate 'Total number of atoms'
  Tabulate (Atoms)
  
  # Start the report
  WriteReport Title,Filename=(resultbase)_report,Text='YASARA Molecular Dynamics Trajectory Analysis for (basename MacroTarget)'
  # Write about the system
  WriteReport Heading,2,'About the simulation'
  if not count block
    WriteReport Paragraph,
      'The trajectory `(MacroTarget)` has been analyzed with YASARA version (Version) over a period of '
      '(simtime) (plottimestring) with (snapshots) snapshots and the (fof) force field. Note that the MD '
      'simulation may have been run with a different force field, but (fof) was used to calculate the '
      'energies in this report. To change this, edit the ForceField setting at the start of this macro. '
  else
    WriteReport Paragraph,
      'The trajectory has been analyzed with YASARA version (Version) in blocks of (blocksnapshots) snapshots and the (fof) force field.'
  WriteReport Paragraph,
    'All plots and pictures in this report [like the simulated system below] are (figurewidth) pixels wide, you '
    'can change the `figurewidth` variable in this macro as needed.'
  # Include a screenshot of the simulated system
  caption='A ray-traced picture of the simulated system. The simulation cell boundary is set to (bound). '
          'Atoms that stick out of the simulation cell will '
  if bound=='periodic'
    caption=caption+'be wrapped to the opposite side of the cell during the simulation.'
  else
    caption=caption+'not be included in the simulation.'
  WriteReport Image,Filename=(resultbase)_sim.png,Style=Figure,Caption=(caption),Delete=Yes
  # Include table with system composition info
  WriteReport Heading,3,'Composition of the system'
  WriteReport Paragraph,'The components of the system are shown in the table below. '
  WriteReport Table,Info,'Composition of the simulated system',.0f,InfoColumnName="Type",DataColumnName="Number"
  # Show the solute
  WriteReport Paragraph,
    'Object (soluteobj) with name `(solutename)` has been identified as the solute and is shown below. '
    'If this is not the intended solute, please change the `soluteobj` variable in this macro. (LigandObjWarning)'
  WriteReport Image,Filename=(resultbase)_solute.png,Style=Figure,Caption='The solute oriented along the major axes.',Delete=Yes
  if slvfirstatm
    WriteReport Paragraph,
      'YASARA detected that the solvent was still part of the solute and thus automatically split them into separate objects at atom (slvfirstatm). '
      'If this is not the intended solvent, please split them manually in the "file://(basename resultbase)_water.sce" file.'
  if ligandsel!=''
    ligandresidues=CountRes (ligandsel)
    ligandatoms=CountAtom (ligandsel)
    WriteReport Heading,3,'The ligand'
    WriteReport Paragraph,
      'A special analysis has been performed for the ligand, chosen (ligandreportstring) with the selection `(ligandsel)`. '
      'The number of residues matching the ligand selection is (ligandresidues), with (ligandatoms) atoms. '
      'To change the ligand selection, edit the `ligandsel` variable at the beginning of this macro.'
    WriteReport Image,Filename=(resultbase)_ligand.png,Style=Figure,
      Caption='A ray-traced picture of the ligand (ligandsel). Bonds are colored by '
              'their order: Gray = 1, blue = 1.25, magenta = 1.33, red = 1.5, orange = 1.66, '
              'bright orange = 1.75, yellow = 2, lime green = 2.5, green = 3 and cyan = 4.',Delete=Yes
# Write about the analyses inside the simulation cell, with periodic boundaries
WriteReport Heading,2,'Analyses inside the simulation cell'
WriteReport Paragraph,
  'This section shows all analyses that have been performed inside the simulation cell, '
  'when all atoms share the common coordinate system of the simulation cell. '
if bound=='periodic'  
  WriteReport Paragraph,
    'Periodic boundaries are active and considered for distance measurements. Calculations '
    'that involve groups of atoms [center of mass, regression lines, enclosing spheres..] '
    'are ambiguous and should be placed in the next section, unless it is known that the '
    'atom group does not drift through a periodic boundary.'
# Write the data to the report
task='WriteReport'
Sim on
AnalyzeInsideCell
Sim off
# Write about the analyses outside the simulation cell, without periodic boundaries
WriteReport Heading,2,'Analyses outside the simulation cell'
WriteReport Paragraph,
  'The following section presents data gathered outside the simulation cell, where each object '
  'has its own local coordinate system and no periodic boundaries are present. Calculations '
  'that involve the interaction between objects [common surface areas, contacts between objects..] '
  'must be placed in the previous section.'
AnalyzeOutsideCell
# Write about the analyses done with respect to the reference structure
AddObj SoluteRef
if refsnapshot=='average'
  text='Analyses performed with respect to the time average structure'
else
  if !refsnapshot
    text='Analyses performed with respect to the starting structure'
  else
    text='Analyses performed with respect to snapshot (refsnapshot)'
WriteReport Heading,2,'(text)'
WriteReport Paragraph,
  '(text) are shown in this section. '
  'These are also done outside the simulation cell, where each object has its own local '
  'coordinate systems and no periodic boundaries are present. To choose another reference '
  'snapshot than (refsnapshot), edit the `refsnapshot` variable at the beginning of this macro.'
AnalyzeChange
RemoveObj SoluteRef
# Loop over the columns (except 1 and 2, the time) to calculate the mean, minimum and maximum values
# and append them to the main table.
for func in 'Mean','Min','Max'
  SelectTab Main
  Tabulate "_"
  Tabulate '(func)'
  for i=3 to tabcolslist1
    vallist()=Tab Main,Column=(i)
    Tabulate ((func) vallist)
for i=1 to tables
  tabname='(resultbase)_analysis'
  if i>1
    tabname=tabname+'_(tabnamelist(i))'
  SaveTab (tabnamelist(i)),(tabname)(id),Format=Text,Columns=(tabcolslist(i)),NumFormat=12.3f,(tabheaderlist(i))

# Special analysis functions
# ==========================
# Additionally calculate and plot the RMSF (Root Mean Square Fluctuation) for every residue per molecule
# Lists of solute residues for RMSF calculations
rmsfreslist()=RMSFRes Obj Solute
residlist()=ListRes Obj Solute
resnumlist()=ListRes Obj Solute,'RESNUMWIC'
sel='Protein or NucAcid and Obj Solute'
# Lists of solute macro molecules for RMSF calculations
molnamelist()=ListMol (sel),Format='Mol MOLNAME'
molnametitlelist()=ListMol (sel),Format='MOLNAME'
mols=count molnamelist
# List of solute macro molecule (MM) residues for RMSF calculations
mmresidlist()=ListRes (sel)
mmresidues=count mmresidlist
if mmresidues
  mmrmsflist()=ShortList rmsfreslist,residlist,mmresidlist
  mmmolnamelist()=ListRes (sel),'Mol MOLNAME'
  mmresnumlist()=ListRes (sel),'RESNUMWIC'
  mmresnumlist()=0+mmresnumlist
  # Get min and max residue number
  minlist()=min mmresnumlist,0+resnummin
  maxlist()=max mmresnumlist,0+resnummax
  rmsfresnummin=max minlist
  rmsfresnummax=min maxlist
  # Make a table of solute macro molecule (MM) residues RMSFs
  MakeTab MMRMSF,2,(mols+1),(rmsfresnummax-rmsfresnummin+1)
  Tab MMRMSF,Set=0
  for i=rmsfresnummin to rmsfresnummax
    Tab MMRMSF,1,(i-rmsfresnummin+1),Set=(i)
  # Enter RMSF values
  for i=1 to mmresidues
    if mmresnumlist(i)>=rmsfresnummin and mmresnumlist(i)<=rmsfresnummax
      for j=1 to mols
        if mmmolnamelist(i)==molnamelist(j)
          Tab MMRMSF,(1+j),(mmresnumlist(i)-rmsfresnummin+1),Set=(mmrmsflist(i))

# Lists of solute HET residues for RMSF calculations
sel='!Protein !NucAcid !Water Obj Solute'
hetresinfolist()=ListRes (sel),"'MOLNAME','RESName``RESNUM','ATOMNUM'"
hetresidues=count hetresinfolist
if hetresidues
  hetresidlist()=ListRes (sel)
  hetrmsflist()=ShortList rmsfreslist,residlist,hetresidlist
  # Make a table of solute HET residues RMSFs
  MakeTab HetRMSF,2,4
  for i=1 to hetresidues
    Tabulate (hetresinfolist(i)),(hetrmsflist(i))

if mmresidues or hetresidues
  WriteReport Heading,2,'Solute residue RMSF'
  WriteReport Paragraph,
    'The Root Mean Square Fluctuation [RMSF] per solute residue is calculated from the average RMSF of its constituting atoms. '
    'The RMSF of atom i with j runing from 1 to 3 for the x, y, and z coordinate of the atom position vector P and '
    'k runing over the set of N evaluated snapshots is given by following formula:'
  WriteReport Image,Filename=(YASARADir)/doc/RMSFAtomResMol1.png,Style=Center,Name="formula_rmsf"
  if mmresidues
    WriteReport Paragraph,'Each graph in the following plot represents one molecule, '
                          'so that you can easily see differences between molecules. '
                          'Note: Residue numbers are not unique, so graphs can overlap. '
    WriteReport Plot,'The Root Mean Square Fluctuation [vertical axis] per solute protein/nucleic acid residue [horizontal axis] '
                      'calculated from the average RMSF of the atoms constituting the residue. '
                      'A RMSF of exactly zero means that that residue number is not present in the molecule. '
                      'Atom RMSF table: "file://(basename resultbase)_rmsf(id).tab", '
                      'residue RMSF table: "file://(basename resultbase)_rmsfres(id).tab"',
                MMRMSF,Width=(figurewidth),Height=480,Title='Solute protein/nucleic acid residue RMSF',
                XColumn=1,YColumn=2,YColumns=(mols),XLabel='Residue number',
                YLabel="RMSF in Angstrom",LegendPos='Outside',Graphname=(molnamelist)
    if hiresplotted
      SavePlot Filename="LastReportPlot_hires",MMRMSF,Width=1600,Height=1200,Title='Solute protein/nucleic acid residue RMSF',
               XColumn=1,YColumn=2,YColumns=(mols),XLabel='Residue number',YLabel="RMSF in Angstrom",Graphname=(molnamelist)
    if mols>1
      # Additionally create per-molecule RMSF plots
      WriteReport Paragraph,'In case the plot above is too crowded, '
                            'the per-residue RMSF values are shown separately for all (mols) molecules in the following plots:'
      for i=1 to mols
        WriteReport Plot,'The Root Mean Square Fluctuation [vertical axis] per solute protein/nucleic acid residue [horizontal axis] '
                          'calculated from the average RMSF of the atoms constituting the residue. '
                          'A RMSF of exactly zero means that that residue number is not present in the molecule. '
                          'Atom RMSF table: "file://(basename resultbase)_rmsf(id).tab", '
                          'residue RMSF table: "file://(basename resultbase)_rmsfres(id).tab"',
                    MMRMSF,Width=(figurewidth),Height=480,
                    Title='Solute protein/nucleic acid residue RMSF of molecule (molnametitlelist(i))',
                    XColumn=1,YColumn=(1+i),YColumns=1,XLabel='Residue number',
                    YLabel="RMSF in Angstrom",LegendPos='Outside',Graphname=(molnamelist(i))
        if hiresplotted
          SavePlot Filename="LastReportPlot_hires",MMRMSF,Width=1600,Height=1200,
                   Title='Solute protein/nucleic acid residue RMSF of molecule (molnametitlelist(i))',
                   XColumn=1,YColumn=(1+i),YColumns=1,XLabel='Residue number',YLabel="RMSF in Angstrom",Graphname=(molnamelist(i))
  if hetresidues
    WriteReport Table,HetRMSF,'RMSF in Angstrom for non-protein/nucleic acid residues in the solute.',InfoColumnName='Mol',DataColumnName='Residue','First atom','RMSF[A]'

if count rdfsellist==4
  # Additionally calculate and include the radial distribution function (RDF) in the report.
  rdflist()=RDF
  MakeTab RDF,2,2
  for i=1 to count rdflist
    Tabulate (rdfsellist4*i),(rdflist(i))
  SaveTab RDF,(resultbase)_rdf(id),Format=Text,Columns=2,NumFormat=6.3f,
          'Radial distribution function with parameters (rdfsellist) as a function of the radial distance in A'
  # Write to report
  WriteReport Heading,2,'Radial Distribution Function'
  WriteReport Paragraph,
    'The Radial Distribution Function is calculated by first by determining distances '
    'between all (rdfsellist1) - (rdfsellist2) pairs and sorting them in (rdfsellist3) bins '
    'with a bin width of (rdfsellist4) A. The RDF is then computed with the following formula:'
  WriteReport Image,Filename=(YASARADir)/doc/RDF1.png,Style=Center,Name="formula_rdf"
  WriteReport Paragraph,
    'The RDF in bin `i` is thus calculated from the number of `CountsInBin i` divided by `Atoms1` '
    '[the number of atoms matching the first selection (rdfsellist1)] times the volume of the '
    'shell corresponding to bin i. To change the selection for the RDF edit the variable "rdfsel" '
    'at the beginning of this macro.'
  WriteReport Plot,'Radial Distribution Function [vertical axis] for (rdfsellist1) - (rdfsellist2) pairs '
                   'as a function of the radial distance r in units of A [horizontal axis] '
                   'calculated from (rdfsellist3) bins with (rdfsellist4) A bin width. '
                   'A table with the raw data is available here: "file://(basename resultbase)_rdf(id).tab"',
              RDF,Width=(figurewidth),Height=480,Title='RDF for (rdfsellist1) - (rdfsellist2) pairs',
              XColumn=1,YColumn=2,YColumns=1,XLabel='Radial distance r in A',YLabel="RDF(r)",LegendPos='Outside',Graphname='RDF'
  if hiresplotted
    SavePlot Filename="LastReportPlot_hires",RDF,Width=1600,Height=1200,Title='RDF for (rdfsellist1) - (rdfsellist2) pairs',
             XColumn=1,YColumn=2,YColumns=1,XLabel='Radial distance r in A',YLabel="RDF(r)"

if espmethod!=''
  # Calculate and save average ESP
  espvallist()=espvallist()/snapshots
  SelectTab (esptbl),Column=1
  Tabulate (espvallist)
  SaveTab (esptbl),(resultbase)_average(id).cube,Format=Gaussian,NumFormat=12.4e
  # Visualize average ESP
  AddObj SoluteRef
  SwitchObj not SimCell SoluteRef,off
  Cell (cellsize)
  esp=LoadESP (resultbase)_average(id).cube
  PointPar Radius=32.00
  SaveScreenshot 'esp(id)','electrostatic potential'
  DelObj (esp)
  Switch On
  RemoveObj SoluteRef
  # Add ESP to the report
  WriteReport Heading,2,'Electrostatic potential'
  WriteReport Paragraph,
    'The electrostatic potential along a grid inside the simulation cell has been calculated '
    'for each snapshot and saved in Gaussian cube format in "file://(basename resultbase)(00000+firstsnapshot).cube.gz" '
    'to "file://(basename resultbase)(lastsnapshot).cube.gz" and the average in "file://(basename resultbase)_average(id).cube".'
  if (espmethod=='PBS')
    WriteReport Paragraph,
      'The chosen method "PBS" uses APBS, the Adaptive Poisson-Boltzmann Solver [Baker et al., 2001], and '
      'thus the Poisson-Boltzmann equation to include solvent and counter ion effects implicitly. '
      'This method treats the cell as non-periodic. '
  else
    WriteReport Paragraph,
      'The chosen method "PME" is based on the Particle Mesh Ewald approach [Essman et al., 1995]: '
      'the atom charges are distributed on a grid, then a fast Fourier transform is made to obtain '
      'the reciprocal space portion of the Coulomb potential. Compared to the real potential, '
      'this is a smoothed representation without short-range noise and singularities [Krieger, Nielsen et al., 2006]. '
      'The resulting ESP considers all atoms in the cell and nothing else, so it is a vacuum potential without '
      'implicit solvent. If the cell is not neutral, the PME method adds a uniform "neutralizing plasma" '
      'throughout the cell to avoid artifacts. This method requires that the cell is periodic. '
  WriteReport Image,(resultbase)_esp(id).png,Style=Figure,
    'Visualization of the average electrostatic potential at each point on a grid throughout the simulation cell with transparent dots. '
    'Red dots indicate a negative, blue dots a positive potential. '

# IMPORTANT: This function must be the last. Place additional functions above this point.
if dccmsel!='' and casel!='None'
  # Additionally calculate and show the dynamic cross-correlation matrix (DCCM).
  # This matrix correlates the displacements from the time average structure,
  # see the documentation of the 'DCCM' command for details.
  # First get the number of selected units, i.e. the rows/columns in the matrix 
  if !dccmunits
    WriteReport Heading,2,'Dynamic Cross-Correlation Matrix'
    WriteReport Paragraph,
      'The DCCM selection (dccmsel) did not match any atoms and therefore the DCCM could not be calculated. '
      'You can change the selection by editing the `dccmsel` variable at the beginning of the macro.'
  else
    # Take the time average structure as the start object to superpose onto
    DelObj SoluteRef
    refobj=DuplicateObj Solute
    NameObj (refobj),SoluteRef
    RemoveObj SoluteRef
    # Loop over the snapshots a second time to calculate the displacements from the time average
    i=00000+firstsnapshot
    last=0
    if count block
      snapshots=blocksnapshots
    while !last and snapshots!=0
      # Load next snapshot from SIM or XTC trajectory
      if format=='sim'
        # Set last (last file) to 1 if last snapshot loaded
        sim=FileSize (MacroTarget)(i+1).sim
        if not sim
          last=1
        LoadSim (MacroTarget)(i)
      else
        last=Load(format) (MacroTarget),(i+1)
      Sim Pause
      ShowMessage 'Calculating dynamic cross-correlation matrix, analyzing snapshot (0+i)...'
      Wait 1
      Sim Off
      # Superpose snapshot on the time average structure
      AddObj SoluteRef
      SupAtom (casel),(caselref)
      # Add the current displacements to an internal table to obtain the DCCM
      AddDisp(dccmsel) and Obj Solute,(dccmsel) and Obj SoluteRef
      RemoveObj SoluteRef
      # Check for last snapshot
      if format!='sim' and !last
        for j=2 to snapshotstep
          last=Load(format) (MacroTarget),(i+j)
          if last
            break
      # Next snapshot
      if snapshots!='all'
        snapshots=snapshots-1
      i=i+snapshotstep
    HideMessage
    # Store the DCCM in a table named DCCM
    MakeTab DCCM,Dimensions=2,Columns=(dccmunits)
    Tabulate DCCM
    # Visualize the DCCM
    coordsys=CoordSys
    pointwidth=1.
    height=5.
    dccmobj1=ShowTab DCCM,Width=(pointwidth),Range=(height),Min=-1,MinCol=(dccmcol(1)),Max=1.0,MaxCol=(dccmcol(2))
    # By default, ShowTab shows the minimum at Z=0, move so that correlation 0 is at Z=0
    MoveMesh (dccmobj1),Z=(-height*0.5)
    # Visualize the zero level with a flat DCCM wireframe
    dccmobj2=ShowTab DCCM,Width=(pointwidth),Range=0,Min=-1,Max=1.0
    dccmobj3=ShowWireObj (dccmobj2),Static,Mesh=Solid
    DelObj (dccmobj2)
    PointPar Radius=0.5,Plastic=No
    NameObj (dccmobj3),ZeroLevel
    ScaleObj (dccmobj1) (dccmobj3),X=(coordsys)
    RotateObj (dccmobj1) (dccmobj3),X=180
    # Create a text object with the residue names and the table header
    textwidth=pointwidth*dccmunits*2
    idlist()=List(dccmsel) and Obj Solute,Format='MOLNAME RESName RESNUM' 
    if textwidth<4500
      textobj1=MakeTextObj Units,Width=(textwidth),Height=(textwidth)
      Font Arial,Height=(pointwidth*0.6),Color=Yellow,Depth=0.5,DepthCol=Red
      for i=1 to dccmunits
        PosText X=(textwidth*0.5+pointwidth*-0.5*(dccmunits+8)),
                Y=(textwidth*0.5+pointwidth*(0.5*dccmunits-i)),justify=left
        Print (idlist(i))
      # Duplicate the labels at the bottom
      textobj2=DuplicateObj (textobj1)
      RotateObj (textobj2),Z=(90*coordsys)
    # Change the target of Print back from the text object to the console
    PrintCon
    # Calculate the projected height of the DCCM and move it to a spot where it fits the screen exactly
    r=RadiusObj (dccmobj1)
    r=r+5
    _,scrsizey,scrscale=ScreenSize
    s=PixToA*(scrsizey/scrscale-50)*0.5/sqrt (r*r/2)
    PosObj all,Z=(EyeDis/s-EyeDis)
    DelObj not (dccmobj1) (dccmobj3) Units Solute
    if not count block
      # Make a ray-traced picture of visualized DCCM
      SwitchObj Solute,off
      SwitchObj ZeroLevel,off
      SaveScreenshot 'dccm1','DCCM'
      SwitchObj Solute,on
      SwitchObj ZeroLevel,on
    # Save the matrix
    SaveTab DCCM,(resultbase)_dccm(id),Format=Text,Columns=(dccmunits),NumFormat=6.3f,
            'Dynamic Cross-Correlation Matrix for (dccmunits) selected units'
    # Save the matrix to be printed in the report
    MakeTab DCCM2,Dimensions=2,Columns=(dccmunits+1)
    # Top left corner
    Tabulate "DCCM"
    # Top row with all residue IDs
    Tabulate List(dccmsel) and Obj Solute,Format='MOLNAME\\n RESName\\n RESNUM' 
    # Fill the other rows
    idlist()=List(dccmsel) and Obj Solute,Format='MOLNAME``RESName``RESNUM'
    for i=1 to dccmunits
      # First the residue ID
      Tabulate '(idlist(i))'
      # Tabulate row i of table DCCM
      Tabulate Tab DCCM,Row=(i)
    # Then show ribbon and make sure DCCM atoms are visible so that arrows are visible
    Style Ribbon
    ShowAtom (dccmsel)
    # Get a list of the first atom in each unit in dccmsel
    idlist()=List(dccmsel) and Obj Solute
    # Assign the DCCM to a list for quick access
    dccmlist()=DCCM
    # Show and style solute object
    ShowSolute
    Style Ribbon,Stick
    HideRes Protein
    Show(dccmsel)
    # Show blue and red lines between strongly anti- and correlated pairs
    # If the selcetion consists of residues or molecules, the first visible atom is selected for drawing arrows
    for i=1 to dccmunits
      for j=i+1 to dccmunits
        corr=dccmlist((i-1)*dccmunits+j)
        if corr<-dccmcut
          col='Blue'
        elif corr>dccmcut
          col='Red'
        else
          continue
        izeroed=zeroed dccmunits+i
        jzeroed=zeroed dccmunits+j
        if not count atm(izeroed)
          dccmatm(izeroed)=ListAtom (idlist(i)) Visible
        if not count dccmatm(jzeroed)
          dccmatm(jzeroed)=ListAtom (idlist(j)) Visible
        corrlist(count corrlist+1)='(-abs corr) (dccmatm(izeroed)) (dccmatm(jzeroed)) (col)'
    correlations=count corrlist
    if correlations
      if correlations>10000
        # Show at max 10000 correlations
        correlations=10000
      # Sort by correlation
      corrlist()=sort corrlist
      for i=1 to correlations
        _,atm1,atm2,col=split corrlist(i)
        # Draw correlations
        ShowArrow Start=AtAtom,(atm1),End=AtAtom,(atm2),Radius=0.1,Heads=0,Color=(col)
    # Save the solute with arrows
    SaveYOb Solute,(resultbase)_dccm(id)
    if not count block
      SaveScreenshot 'dccm2','motion correlations'
    SwitchAll On
    DelObj Solute
    # Save the visualized matrix
    NumberObj all
    PosObj all,Z=(EyeDis/s-EyeDis)
    SaveSce (resultbase)_dccm(id)
    # Write to report
    WriteReport Heading,2,'Dynamic Cross-Correlation Matrix'
    WriteReport Paragraph,
      'The dynamic cross-correlation matrix [DCCM] is a square matrix, whose rows and columns '
      'match the selected units `(dccmsel)`. To change this selection, edit the `dccmsel` variable '
      'at the beginning of this macro. The DCCM shows how the movements of all selected pairs correlate. '
      'The values in the DCCM range from -1 [perfectly anti-correlated] to +1 [perfectly correlated]. '
      'The values along the diagonal are always +1 [because the motion of an atom is perfectly correlated to itself]. '
      'The DCCM element for units i and j is obtained with the following formula: '
    Writereport Image,Filename=(YASARADir)/doc/DCCM1.png,Style=Center,Name="formula_dccm"
    WriteReport Paragraph,
      'Here `d` is the displacement between the current position and the average position of the '
      'selected unit, and the angle brackets indicate the average over all samples. '
      'The highest correlations off the diagonal can often be found for bridged cysteines.'
    if not count block
      WriteReport Paragraph,'The image below shows the correlation directly in the solute object:'
      WriteReport Image,(resultbase)_dccm2.png,Style=Figure,
                  'Blue and red lines are shown between (correlations) strongly anti- and correlated residue pairs. '
                  'To change the threshold value for the correlation lines edit the `dccmcut` variable at the beginning of this macro. '
                  'To look at this structure interactively, open the file "file://(basename resultbase)_dccm.yob" in YASARA.',Delete=yes
      WriteReport Paragraph,
        'In the image below, the DCCM is visualized with colors ranging from '
        '(dccmcol(1)) [-1, fully anti-correlated] to (dccmcol(2)) [+1, fully correlated]. '
      WriteReport Image,(resultbase)_dccm1.png,Style=Figure,
        'Visualization of the dynamic cross-correlation matrix. Open the file "file://(basename resultbase)_dccm.sce" in YASARA '
        'to look at this matrix visualization interactively. '
        'In the scene file, the zero level [0, not correlated] is indicated with a wire-frame grid. ',Delete=Yes
    else
      WriteReport Paragraph,
      'Open the file "file://(basename resultbase)_dccm(id).sce" in YASARA to look at the dynamic cross-correlation matrix visualization interactively. '
      'In the scene file, the DCCM is visualized with colors ranging from (dccmcol(1)) [-1, fully anti-correlated] to (dccmcol(2)) [+1, fully correlated] and '
      'the zero level [0, not correlated] is indicated with a wire-frame grid. '
      'To view the correlation directly in the solute object open the file "file://(basename resultbase)_dccm(id).yob" in YASARA. '
      'Blue and red lines are shown between strongly anti- and correlated residue pairs. '
      'To change the threshold value for the correlation lines edit the `dccmcut` variable at the beginning of this macro.'
    captiontext='Dynamic cross-correlation matrix. '
                'The full table is also available in text format, you need a proper text editor '
                'without line wrapping to look at this file: "file://(basename resultbase)_dccm(id).tab". '
    if tabrowsmax
      captiontext=captiontext+'Note: At most (tabrowsmax) rows of the DCCM are shown above. '
                              'Change the `tabrowsmax` variable in the macro to adjust this number. '
    WriteReport Table,DCCM2,'(captiontext)',RowsMax=(tabrowsmax)
    
# Additional files block
Writereport Heading,2,'Additional files'
WriteReport Paragraph,'The following additional files have been created:'
Writereport Heading,3,'The main data table'
WriteReport Paragraph,
  'The main table contains all collected data in a single file. The column names match the '
  'names used above for graphs in plots and columns in tables. You can find a more detailed '
  'explanation of this table in the user manual at Recipes > Run a molecular dynamics simulation > '
  'Analyzing a trajectory. If you parse this file automatically, keep in mind that the number of '
  'columns can change any time, so you have to use the names in the first table row to find the '
  'columns of interest: "file://(basename resultbase)_analysis(id).tab"'
if (tables>1)
  Writereport Heading,3,'Extra data tables'
  WriteReport Paragraph,'User defined extra data tables:'
  for i=2 to tables
    WriteReport Paragraph,'"file://(basename resultbase)_analysis_(tabnamelist(i))(id).tab"'
if (plotmaptabs>1)
  Writereport Heading,3,'Per-atom and per-residue data tables'
  WriteReport Paragraph,'Data of the per-atom and per-residue plots:'
  for i=1 to plotmaptabs
    WriteReport Paragraph,'"file://(basename plotmaptablist(i)).tab"'
Writereport Heading,3,'The structures'
paragraphtext=''
if count block
  paragraphtext=' of the current block'
WriteReport Paragraph,
  'The `time averaged structure`(paragraphtext) in PDB format: "file://(basename resultbase)_average(id).pdb"'
WriteReport Paragraph,
  'The `snapshot with the minimum solute energy`. Either just the solute in PDB format '
  '"file://(basename resultbase)_energymin(id).pdb", or the complete system including '
  'solvent as a YASARA scene "file://(basename resultbase)_energymin(id).sce".'
if last
  WriteReport Paragraph,
    'The `last snapshot` of the simulation. Either just the solute in PDB format '
    '"file://(basename resultbase)_last.pdb", or the complete system including solvent '
    'as a YASARA scene "file://(basename resultbase)_last.sce"'
if rmsdmin
  WriteReport Paragraph,'One `representative structure for each cluster` in YOb format:'
  for i=1 to clustermembers
    WriteReport Paragraph,'"file://(basename clusterfilenamelist(i)).yob"'
Writereport Heading,3,'The RMSF tables'
WriteReport Paragraph,
  'A table that lists the Root Mean Square Fluctuations [RMSFs] of all atoms in [A] is available '
  'here: "file://(basename resultbase)_rmsf(id).tab". The RMSFs have also been converted '
  'to B-factors and stored in the B-factor field of the time-average structure above. '
WriteReport Paragraph,
  'A table with average atom RMSFs per residue can be found here: "file://(basename resultbase)_rmsfres(id).tab".'
if hiresplotted and (last or not count block)
  Writereport Heading,3,'High resolution plots'
  WriteReport Paragraph,
    'To facilitate publication, high resolution versions of the plots above have been '
    'created with a 4:3 aspect ratio suited for printing in a single column of a typical '
    'journal article. Just look at the figure number above to find the right file:'
  for i=1 to 9999
    size=FileSize (resultbase)_report_figure(i)_hires.png
    if size
      WriteReport Paragraph,'"file://(basename resultbase)_report_figure(i)_hires.png"'

if last or not count block
  # End report
  WriteReport End
  HideMessage
  if runWithMacro and ConsoleMode
    # Exit YASARA if this macro was provided as command line argument in console mode
    Exit
  elif not ConsoleMode
    # Open report in browser if this macro was not run in console mode
    ShowURL file://(resultbase)_report.html

# ADD A PLOT OR TABLE TO THE HTML REPORT
# ======================================
# 'xvalue' is the number of the value in valuelist that will be used as X-value for scatter plots
# 'valuelist' is either a list of values to plot or tabulate, or a command to run that returns such a list of values
# 'title' is the title of the plot or table
# 'xlabel' is the label of the X-axis for scatter plots
# 'ylabel' is the label of the Y-axis
# 'graphnamestr'/'datacolumnstr' is a string with the names of the graphs/columns to create in the plot/table, joined with ' '
# 'task' is a global variable that tells us what to do: 'Tabulate' to collect data, 'AddHeader'
# to add the graph/column names to global variable 'header', and 'WriteReport' to write the actual report
# (in the latter case the global variable 'ycolumn' is the starting column to plot or tabulate).
def Plot valuelist(),title,ylabel,graphnamestr
  ShowData 'Plot',0,valuelist,title,"None",ylabel,graphnamestr

def ScatterPlot valuelist(),xvalue,title,xlabel,ylabel,graphnamestr
  ShowData 'Plot',xvalue,valuelist,title,xlabel,ylabel,graphnamestr

def WriteTable valuelist(),title,datacolumnstr
  ShowData 'Main',0,valuelist,title,"None","None",datacolumnstr

def WriteExtraTable tablename,valuelist(),datacolumnstr
  ShowData tablename,0,valuelist,"None","None","None",datacolumnstr

def ShowData resulttype,xvalue,valuelist(),title,xlabel,ylabel,graphnamestr
  global resultbase,task,hiresplotted,header,xcolumn,ycolumn,plottimestring,tabtimestring,figurewidth,tabrowsmax,snapshots,tabnamelist,tabheaderlist,tabcolslist,tables,id
  graphsmax=30
  graphcolorlist()=''
  command=''
  if type valuelist1=='StrongString'
    # User provided a command to run
    if count valuelist!=1
      RaiseError 'When plotting the output of a YASARA command as "(title)", exactly one command must be provided, not (count valuelist)'
    # Run the command and collect the output as 'valuelist'
    command=''+valuelist1
    valuelist()=(command)
  values=count valuelist
  # Check if valuelist is empty
  if !values
    RaiseError 'The list of values for (resulttype) "(title)" is empty. Please check the selection'
  if values>graphsmax and resulttype=='Plot'
    RaiseError 'Plot "(title)" has (values) graphs. The maximum allowed number of graphs per plot is (graphsmax)'
  # Check the graph names
  graphnamelist()=split graphnamestr
  graphnames=count graphnamelist
  if values<graphnames
    # Too many graphnames, delete them
    for i=graphnames to values+1 step -1
      DelVar graphnamelist(i)
  elif graphnamelist(graphnames)=='Auto'
    # Automatic naming, add graph names if necessary
    for i=graphnames to values
      graphnamelist(i)='Graph(i)'
  elif values>graphnames
    if resulttype=='Plot'
      valuename='graph'
    else
      valuename='data column'
    RaiseError 'You tried to add (values) (valuename)s to the (resulttype) "(title)", but provided only (graphnames) (valuename) names "(graphnamestr)". '
                'Choose "Auto" as the last or only one (valuename) name for automatic naming'
  # Check which task to perform
  if task=='Tabulate'
    SelectTab Main
    if resulttype!='Plot' and resulttype!='Main'
      SelectTab '(resulttype)'
    for i=1 to count valuelist
      # A WeakString can contain commas so it must be protected with quotes
      if type valuelist(i)=='WeakString'
        Tabulate '(valuelist(i))'
      else
        Tabulate (valuelist(i))
  elif task=='AddHeader'
    # Add graph names to table header strings, which will be the first line of the saved tables
    table=1
    if !tables or (resulttype!='Plot' and resulttype!='Main')
      tables=tables+1
      table=tables
      tabnamelist(table)=resulttype
      if resulttype=='Plot'
        tabnamelist(table)='Main'
      tabcolslist(table)=2
      tabheaderlist(table)="    Time[ps]     Time[ns]"
    for i=1 to values
      graphname=graphnamelist(i)
      graphnamelen=strlen graphname
      if graphnamelen>12
        RaiseError 'Graph or table column name "(graphname)" is (graphnamelen) characters long, but at most 12 characters are allowed'
      tabheaderlist(table)=tabheaderlist(table)+" "*(13-graphnamelen)+graphname
    tabcolslist(table)=tabcolslist(table)+values
  elif task=='WriteReport' and (resulttype=='Plot' or resulttype=='Main')
    # Write the actual report
    WriteReport Heading,3,'(title)'
    charlist()=crack graphnamelist1
    code=ord '(upper charlist1)'
    if code>=ord 'A' and code<=ord 'Z' and ' ' not in charlist and '.' not in charlist
      # Variables must start with a letter and not contain spaces or dots
      global (graphnamelist1)Text
      if isfunction Print(graphnamelist1)Info
        # Macro contains a function that writes infos to the report
        Print(graphnamelist1)Info
      elif count (graphnamelist1)Text
        # Macro contains a variable with infos to write
        WriteReport Paragraph,'((graphnamelist1)Text)'
    if resulttype=='Plot'
      # Create plot
      if xvalue
        caption=title+' [vertical axis] as a function of (graphnamelist(xvalue)) [horizontal axis]'
        DelVar graphnamelist(xvalue)
        plotvalues=values-1
        plottype='Scatter'
        plotxcolumn=ycolumn-1+xvalue
      else
        caption=title+' [vertical axis] as a function of simulation time [horizontal axis]'
        xlabel='Simulation time in (plottimestring)'
        plotvalues=values
        plottype='Line'
        plotxcolumn=xcolumn
      if command!=''
        caption=caption+', obtained with the command "(command)"'
      caption=caption+'.'
      # Check if one graph covers another or graph is all zero
      if values>1 and !xvalue
        text=''
        textlist()=' Graph',', graph',' graph',' and graph',' has',' have'
        idx=0
        # First, check if there are all zero graphs
        for i=1 to values
          # Initilaize donelist for the covering graph part
          donelist(i)=0
          # Get all values of one column
          columnlist()=Tab Main,(ycolumn+i-1)
          # Add sum of all values to sumlist
          sumlist(i)=Sum columnlist
          if !sumlist(i)
            # If the sum of all values is zero then we asume that all values are zero
            if idx<2
              idx=idx+1
            # Add text
            text=text+textlist(idx)+' `(graphnamelist(i))`'
        if idx
          # Add ending text
          text=text+textlist(idx+4)+' all zero values.'
        # Second, check if one graph covers another
        idx=0
        for i=values to 2 step -1
          if !donelist(i)
          # Check if graph covers another only if that graph itself is not covered and already flaged in donelist
            for j=i-1 to 1 step -1
              if sumlist(j)==sumlist(i) and !donelist(j)
                # If the sum of all values between two graphs is equal we asume that the graphs are identical
                donelist(j)=1
                if idx<2
                  # Add starting text
                  text=text+textlist(idx+1)+' `(graphnamelist(i))` completely covers '
                  idx=3
                # Add text for every graph that is covered
                text=text+textlist(idx)+' `(graphnamelist(j))`'
                idx=4
            if idx
              # Reset index
              idx=1
        if idx
          # Add ending text
          text=text+', they share the same values.'
        if text!=''
          # Add text to caption
          caption=caption+' Note:'+text
      if graphnamelist1=='TotalEnergy' and snapshots>1
        tempval1=Tab Main,(ycolumn),1
        tempval2=Tab Main,(ycolumn),2
        if tempval1<tempval2
          caption=caption+' Note: The first value of the plot [(0.00+tempval1)], coming from the energy minimized starting structure, '
                          'has been replaced with the second value of the plot [(0.00+tempval2)] to show this plot with a smaller energy range '
                          'and thus a higher resolution. '
          Tab Main,(ycolumn),1,Set=(tempval2)
      elif graphnamestr=='CellLengthX CellLengthY CellLengthZ'
        # Set cell length graph colors to match cell colors
        graphcolorlist()='Red','Green','Blue'
      elif graphnamestr=='SelHBonds RecHBonds SlvHBonds TotHBonds'
        # Set colors to match the description
        graphcolorlist()='Blue','Green','Red','Grey'
      elif graphnamestr=='Helix Sheet Turn Coil Helix310 HelixPi'
        # Set graph color for secondary structure content plot
        secstrnamelist()='Helix','Sheet','Turn','Coil','Helix310','HelixPi'
        for i=1 to count secstrnamelist
          graphcolorlist(i)=ColorPar SecStr,(secstrnamelist(i))
      WriteReport Plot,Caption='(caption)',Main,Width=(figurewidth),Height=480,Title='(title)',
                  Type=(plottype),XColumn=(plotxcolumn),YColumn=(ycolumn),YColumns=(plotvalues),XLabel=(xlabel),
                  YLabel=(ylabel),LegendPos='Outside',Graphname=(graphnamelist),Graphcol=(graphcolorlist)
      if hiresplotted
        SavePlot Filename="LastReportPlot_hires",Main,Width=1600,Height=1200,Title='(title)',
                 Type=(plottype),XColumn=(plotxcolumn),YColumn=(ycolumn),YColumns=(plotvalues),XLabel=(xlabel),
                 YLabel=(ylabel),Graphname=(graphnamelist),Graphcol=(graphcolorlist)
      if graphnamelist1=='TotalEnergy' and snapshots>1
        Tab Main,(ycolumn),1,Set=(tempval1)
    else
      # Create table
      caption=title+' as a function of simulation time [first column]'
      if command!=''
        caption=caption+', obtained with the command "(command)"'
      caption=caption+'.'
      if tabrowsmax
        caption=caption+' Note: At most (tabrowsmax) table rows are shown. '
                        'Change the `tabrowsmax` variable in the macro to adjust the number of shown table rows. '
                        'The full table can be found in "file://(basename resultbase)_analysis(id).tab".'
      WriteReport Table,Main,Caption='(caption)',NumFormat=.2f,RowsMax=(tabrowsmax),
                  InfoColumn=(xcolumn),DataColumn=(ycolumn),DataColumns=(values),
                  InfoColumnName='(tabtimestring)',DataColumnName=(graphnamelist)
    ycolumn=ycolumn+values

# PLOT A RESIDUE/TIME HEATMAP
# ============================
# A heatmap plot with time as X- and residue number as Y-coordinate will be added to the report and saved to disk.
def PlotRes name,sellist(),valuelist(),valmin,valmax,title,graphnamestr,graphcolstr
  PlotMap 'Res',name,sellist(),valuelist(),valmin,valmax,title,graphnamestr,graphcolstr

# PLOT AN ATOM/TIME HEATMAP
# =========================
# A heatmap plot with time as X- and atom number as Y-coordinate will be added to the report and saved to disk.
def PlotAtom name,sellist(),valuelist(),valmin,valmax,title,graphnamestr,graphcolstr
  PlotMap 'Atom',name,sellist(),valuelist(),valmin,valmax,title,graphnamestr,graphcolstr

# PLOT A HEATMAP
# ==============
def PlotMap type,name,sellist(),valuelist(),valmin,valmax,title,graphnamestr,graphcolstr
  global resultbase,task,figurewidth,hiresplotted,plottimestring,simtime,xcolumn,plotmaptabs,plotmaptablist(),(name)Text,id
  
  selformat=''
  if type=='Res'
    selformat='RESNUMWIC'
  numlist()=List(type) (join sellist),(selformat)
  mollist()=ListMol (type) (join sellist)
  start=1
  for i=1 to count mollist
    selections=Count(type) (mollist(i)) and (join sellist)
    if selections
      if task=='AddHeader'
        # Create tables
        MakeTab (name)Mol(i),2,(selections+2)
        Tabulate 'Time_ps','ns\(type)Num'
        # Tabulate the selection numbers
        for j=start to start+selections-1
          Tabulate (numlist(j))
        start=start+selections
      elif task=='Tabulate'
        # Tabulate data
        SelectTab (name)Mol(i)
        Tabulate (simtime/1000),(simtime/1000000)
        # Tabulate values
        for j=start to (start+selections-1)
          Tabulate (valuelist(j))
        start=start+selections
        SelectTab Main
      elif task=='WriteReport'
        if i==1
          # Write text to report only once
          WriteReport Heading,3,'(title)'
          if isfunction Print(name)Info
            # Macro contains a function that writes infos to the report
            Print(name)Info
          elif count (name)Text
            # Macro contains a variable with infos to write
            WriteReport Paragraph,'((name)Text)'
        molname=ListMol (mollist(i)),Format='MOLNAME'
        moltitle='(title) of molecule (molname)'
        plotmaptabs=plotmaptabs+1
        typename=lower type
        plotmaptablist(plotmaptabs)='(resultbase)_plot(typename)_(name)Mol(molname)(id)'
        if graphnamestr==''
          # Generate legend if none provided
          vals=valmax-valmin+1
          for j=1 to vals
            graphnamelist(j)='(j+valmin-1)'
          graphnames=count graphnamelist
          legendcaption=''
        else
          graphnamelist()=split graphnamestr
          graphnames=count graphnamelist
          legendcaption='. Values 1-(graphnames) in the table correspond to the (graphnames) labels in the plot legend.'
        if (count graphnamelist==2 and graphnamelist1=='Min' and graphnamelist2=='Max')
          # Continuous heatmap
          used=1
          meanplotted=1
          legendcaption=''
          # If valmin/valmax is provided then values will be capped for the plot but not in the saved data table
          graphnamelist1='Min (valmin)'
          graphnamelist2='Max (valmax)'
          caption='data'
          # Save data table using floats
          SaveTab (name)Mol(i),(plotmaptablist(plotmaptabs)),NumFormat=12.3f,(moltitle)
        else
          # For the plot, cap values with valmin and valmax, shift all values so they start with 1
          used=0
          meanplotted=0
          cap()='',''
          shift=1-valmin
          tablist()=Tab (name)Mol(i)
          cols=selections+2
          rows=count tablist/cols
          # List the percentages of occurrence of each value for each selection
          for j=1 to graphnames*cols
            percentlist(j)=0
          for j=2 to rows
            for k=1 to cols
              if k<3
                # Round time stamps down so that the micro seconds listed in the table change only
                # whenever a new micro second of simulation time is reached,
                # as all values saved in the table file will be integers.
                time=Tab (name)Mol(i),(k),(j)
                time=time//1
                Tab (name)Mol(i),(k),(j),Set=(time)
              else
                idx=(j-1)*cols+k
                if tablist(idx)<valmin
                  tablist(idx)=valmin
                  cap1='<='
                elif tablist(idx)>valmax
                  tablist(idx)=valmax
                  cap2='>='
                tablist(idx)=tablist(idx)+shift
                if tablist(idx)>1
                  used=1
                num=k+(tablist(idx)-1)*cols
                percentlist(num)=1.+percentlist(num)
          # Adapt capped graph names
          if graphnamelist1=='(valmin)'
            graphnamelist1='(cap1)(valmin)'
          if graphnamelist(graphnames)=='(valmax)'
            graphnamelist(graphnames)='(cap2)(valmax)'
          for j=1 to graphnames
            idx=(j-1)*cols
            percentlist(idx+1)=""+'(graphnamelist(j))'
            percentlist(idx+2)="%"
            for k=3 to cols
              # Convert count to percentage and round so sum will always be 100
              percentlist(idx+k)=round (percentlist(idx+k)*100./(rows-1))
          # Add percentages to table
          SelectTab (name)Mol(i)
          Tabulate (percentlist)
          caption='data including percentages'
          # Save data table using integers only
          SaveTab (name)Mol(i),(plotmaptablist(plotmaptabs)),NumFormat=10.0f,(moltitle)
          # Update table with shift and caps
          DelTab (name)Mol(i)
          MakeTab (name)Mol(i),2,(cols)
          tablist1=""+'(tablist1)'
          tablist2=""+'(tablist2)'
          Tabulate (tablist)
          SelectTab Main
        if used
          # Make a plot for each molecule
          WriteReport Plot,'(title) as a function of simulation time [horizontal axis] for each (type) number [vertical axis]. '
                           'A table with the raw (caption) is available here: "file://(basename plotmaptablist(plotmaptabs)).tab"(legendcaption)',
                      (name)Mol(i),Width=(figurewidth),Height=480,Title=(moltitle),Type=Heatmap,
                      XColumn=(xcolumn),YColumn=3,YColumns=(selections),XLabel='Simulation time in (plottimestring)',
                      YLabel='(type) number',LegendPos='Outside',Graphname=(graphnamelist),Graphcol=(split graphcolstr)
          if hiresplotted
            SavePlot Filename="LastReportPlot_hires",(name)Mol(i),Width=1600,Height=1200,Title=(moltitle),Type=Heatmap,
                     XColumn=(xcolumn),YColumn=3,YColumns=(selections),XLabel='Simulation time in (plottimestring)',
                     YLabel='(type) number',Graphname=(graphnamelist),Graphcol=(split graphcolstr)
          if meanplotted
            # For continuous heatmaps plot mean values as function of selection number
            MakeTab (name)Mol(i)_mean,2,2
            for j=1 to count numlist
              vallist()=Tab (name)Mol(i),Column=(2+j)
              # First value is the selection number
              Tabulate (vallist(1))
              DelVar vallist(1)
              # Tabulate mean
              Tabulate (mean vallist)
            SaveTab (name)Mol(i)_mean,(plotmaptablist(plotmaptabs))_mean,NumFormat=12.3f,'Mean (moltitle)'
            WriteReport Plot,'Mean (title) [vertical axis] as a function of (type) number [horizontal axis]. '
                             'A table with the raw data is available here: "file://(basename plotmaptablist(plotmaptabs))_mean.tab".',
                        (name)Mol(i)_mean,Width=(figurewidth),Height=480,Title='Mean (moltitle)',
                        Type=Line,XColumn=1,YColumn=2,YColumns=1,XLabel='(type) number',
                        YLabel='Mean (title)',LegendPos='Outside',Graphname='Mean'
            if hiresplotted
              SavePlot Filename="LastReportPlot_hires",(name)Mol(i)_mean,Width=1600,Height=1200,Title='Mean (moltitle)',
                       Type=Line,XColumn=1,YColumn=2,YColumns=1,XLabel='(type) number',YLabel='Mean (title)'

# SHOW THE ENTIRE SIMULATED SYSTEM
# ================================
def ShowSystem
  global ligandsel 
  SwitchAll On
  alpha,beta,gamma = OriObj SimCell
  RotateAll Y=(-beta)
  RotateAll Z=(-gamma)
  RotateAll X=(-alpha)
  Style BallStick
  Style Ribbon,Stick
  ColorBonds Order
  if ligandsel!=''
    BallRes (ligandsel)
  ZoomAll Steps=0

# SHOW THE SOLUTE OBJECT
# ======================
def ShowSolute
  global ligandsel 
  # Show only the solute and orient along the major axes
  SwitchAll Off
  SwitchObj Solute,On
  NiceOriObj Solute
  # First show hetgroups as balls&sticks, but ions and the ligand as balls
  Style BallStick
  BallAtom all with 0 bonds to all
  if ligandsel!=''
    BallRes (ligandsel)
    # Make sure that the ligand faces to the front
    CenterObj Solute
    TransformObj Solute
    _,_,cenz = GroupCenter (ligandsel) Obj Solute
    if cenz>0
      RotateObj Solute,Y=180
  ZoomObj Solute,Steps=0
  
# SAVE A RAY-TRACED SCREENSHOT
# ============================  
def SaveScreenshot fileid,description
  global figurewidth,resultbase
  ShowMessage 'Creating ray-traced picture of the (description)...'
  Wait 1
  RayTrace Filename=(resultbase)_(fileid).png,X=(figurewidth),Zoom=1,LabelShadow=No,Display=Off,Outline=0,Background=On

# SHORTEN A LIST
# ==============
# Returns a list with those items in longdatalist, whose corresponding item in longidlist is in shortidlist.
# longdatalist and longidlist must have the same length, all items in shortidlist must be in longidlist and in the same order.
def ShortList longdatalist(),longidlist(),shortidlist()
  ids=count longidlist
  idx=1
  for i=1 to count shortidlist
    while longidlist(idx)!=shortidlist(i)
      idx=idx+1
    shortdatalist(i)=longdatalist(idx)
  return (shortdatalist)
