# Parameterizing modified nucleic acids for molecular simulations in the AMBER MD software environment [Article v2.0]

Rodrigo Galindo-Murillo<sup>1</sup>, Akanksha Manghrani<sup>2</sup>, Daniel Roe<sup>3</sup>, Olivia Love<sup>4</sup>, Pablo
D. Dans<sup>5</sup>, Thomas E. Cheatham III<sup>4</sup>, Christina Bergonzo<sup>2*</sup>

<sup>1</sup>Department of Medicinal Chemistry, Ionis Pharmaceuticals, 2855 Gazelle Court,
Carlsbad, California 92010, United States; <sup>2</sup>Institute for Bioscience and Biotechnology
Research, National Institute of Standards and Technology and the University of
Maryland, 9600 Gudelsky Way, Rockville, MD, 20850, United States; <sup>3</sup>Laboratory of
Computational Biology, National Heart Lung and Blood Institute, National Institutes of
Health, Bethesda, MD 20892, United States; <sup>4</sup>Department of Medicinal Chemistry,
College of Pharmacy, University of Utah, 2000 East 30 South Skaggs 306, Salt Lake City,
Utah 84112, United States; <sup>5</sup>Computational Biophysics Group, Department of Biological
Sciences, CENUR Litoral Norte, Universidad de la República, 50000 Salto, Uruguay.
Bioinformatics Unit. Institute Pasteur of Montevideo, Uruguay

***For correspondence:**
christina.bergonzo@nist.gov (CB)

## Abstract

Parameterizing modified nucleic acids is a difficult but necessary task for expanding the simulated space of oligonucleotides, including both naturally occurring structures and those with pharmaceutical relevance. In lieu of expensive and difficult chemical synthesis in the laboratory, computer simulations are often performed to make predictions for sequence and structure effects, as well as downstream critical quality attributes. To enable these simulations, modifications have to be parameterized to faithfully represent their effect on nucleotides. This is a non-trivial process, complicated by the fact that it may be the first thing researchers must figure out before they can build their structures and start their initial simulations. To enable these research projects, we created modXNA, a code that assembles pre-parameterized modules of the base, backbone, and sugar, to create bespoke combinations of modifications. In the following tutorial, we provide background on force field parameterization in the Amber software ecosystem and detail the steps necessary to perform parameterization of modified nucleic acids using modXNA.



## 1. Introduction

Molecular dynamics (MD) force fields rely on additive parameters, including bonded and non-bonded terms, to describe a system using physics-based potentials. What follows is a brief overview of Amber force field parameterization, to establish a baseline understanding and provide context for the modXNA [1] parameterization process. The original Amber force field for proteins and nucleic acids addressed parameterization of torsion and angle parameters, deriving these from experimental data [2]. This force field was a precursor to the current equation used to describe the potential energy, introduced in Cornell et al. 1995 [3] (Figure 1). The Cornell et al. Amber ff94 reparametrized the electrostatics using the restrained electrostatic potential (RESP) fitting protocol [3–5]. Subsequent force field modifications undertook refitting nucleic acid sugar pucker and chi torsional parameters [6]. The most recent nucleic acid parameters for RNA combine the refitting of alpha and gamma torsion dihedrals from parmbsc0 [7] with new chi parameters from the OL3 corrections [8], and new chi, epsilon, and zeta torsions in bsc1 [9] for DNA, as well as revisions to beta, alpha, and gamma torsions in the OL15 [9] and OL21 [10] DNA force fields, respectively. In addition to bonded terms, there are non-bonded terms which depend on representative point charges associated with each atom as well as the combination of their atom type specific van der Waals radii. In order to incorporate modifications into the nucleic acid MD, we will have to create parameters for the force field. The following review contains a comprehensive overview of parameters for modified nucleotides and their parameterization methodology [11].

When developing parameters for new residues, charges are fit to each atom using RESP fitting [3–5]. RESP fitting [5] begins with calculating the molecular electrostatic potential. The Hartree-Fock (HF) level of theory with the 6-31G* basis set is used due to opportune cancellation of error between the overestimation of the polarity of molecules and popular water models’ (TIP3P and SPC) enhanced dipole over gas-phase values. A fixed grid of points is created in the solvent accessible space around the molecule, outside the van der Waals radius of the molecule. Quantum mechanics (QM) is used to calculate the electrostatic potential at each point. Least squares fitting assigns partial charges to each atom in the molecule, with the constraint that the sum
of all charges must be equal to the overall charge of the molecule. An added penalty function (the “restrained” part of RESP) eliminates the overestimation of polar atom’s
calculated electrostatics by scaling down the magnitude of their charges [5]. 

Nucleic acids were originally RESP fit using a two-stage method, which involved splitting each nucleotide into two representative groups, the deoxyribonucleoside base and
dimethylphosphate, representing the backbone phosphate group - only the B-form conformation was considered. In the first stage of fitting, hyperbolic restraints were applied to all heavy atoms. Additional charge constraints were applied to the “joint” regions of each of these functional groups in order to correctly assign charges for the 3’ and 5’ termini, and to enforce an integer charge for the entire system. The sugar atoms, with the exception of C1’ and H1’, were intermolecularly equivalenced, meaning that their charges were standardized in all 4 deoxyribonucleosides. The OP1and OP2 atoms were equivalenced, as well as the hydrogen atoms of the NH2 groups. In the second stage of charge fitting, restraints were increased, freezing heavy atoms and
assigning charges for hydrogens of CH3 and CH2 groups. Equivalencing of C2’H2 and C5’H2 groups establishes the same charges for all hydrogens in the group.

## 2. Scope of the tutorial

Parameterization of modifications, including the designation of new atom types with their corresponding hybridization, must be developed for use in MD simulations. Based on the current form of the force field (Figure 1), we need to develop parameters for bonds, angles, dihedrals, and non-bonded interactions, including van der Waals and electrostatics.

![](C:\Users\cnb11\Desktop\Work\siRNA\figures\figure1.png)

**Figure 1. **AMBER force field equation for molecular dynamics simulations.
Bold, blue font denotes variables in the equation that must
be parameterized by the user.

The updated procedure for parameterization, introduced in modXNA [1], preserves some of the early parameterization methods – specifically, we wanted to make sure that in the absence of a modification, the underlying nucleotide parameters are identical to the current force field, that intermolecular equivalences were preserved, and that joint regions were treated fairly (Figure 2). The split from two components
(a base containing a deoxyribose and a phosphate group) to three components (base, sugar, and backbone) was undertaken in order to modify the smallest possible unit (Figure 3). Three modules are chosen by the user: one backbone, one sugar, and one base. The Catalog hosted online (at [https://modxna.chpc.utah.edu]()) and on the GitHub
repository for this publication outlines the fragments that have already been parameterized. The three-letter “Fragment IDs” for each fragment are designated in the input file in order: [backbone] [sugar] [base].

![](C:\Users\cnb11\Desktop\Work\siRNA\figures\figure2.png)

**Figure 2. **A schematic of the workflow used in modXNA.

![](C:\Users\cnb11\Desktop\Work\siRNA\figures\figure3.png)**Figure 3. **Modules which make up the final nucleotide. All three parts are required for modXNA. The red groups denote “capping groups” that were present during QM calculations to act as placeholders for the other modules that make up a complete nucleotide. These capping groups are removed when modXNA assembles the residue.

In the following tutorial, we outline a few common applications where modified  nucleotides are required: first, incorporating a pseudouridine base in an RNA duplex;
second, incorporating 2’-O-methylated nucleotides in an RNA duplex, including modifying the 5’ and 3’ termini; and third, building a heavily modified antisense oligonucleotide. The first two examples focus on building a parameter/topology file and coordinate file pair from a deposited PDB file with modifications, and the third example builds a parameter/topology file and coordinate file pair for a ssRNA from sequence. In each case, the user will produce a parameter/topology file and coordinate file to use to run molecular dynamics in AmberMD. Information about using the force fields in other software packages can be found at the following link: https://ambermd.org/AmberModels.php.

### 2.1 Using modXNA.sh, prerequisites

ModXNA uses several programs in AmberTools, and the best way to obtain the up-to-date accessory programs it uses is to download and compile AmberTools25 (or later versions) from source, or make use of the conda package for binary installation of AmberTools25 (or later versions). ModXNA can be downloaded from [modxna.chpc.utah.edu](). After choosing the three fragments that make up the backbone, sugar and base, the user runs the modXNA.sh script. The script utilizes Amber’s LEaP (tleap), CPPTRAJ [12] and sander programs. Header lines in each fragment’s mol2 file are used to set flags which strip capping groups, allowing the combination of fragments into a residue. After joint corrections are applied, and charge equivalencing occurs, the full nucleotide residue is built in tleap, and minimized using sander. A final step sees the nucleotide residue loaded into tleap to create an Amber "lib” (library) file. The final output of the program, necessary to build this modified residue in tleap, is the XXX.lib file, where “XXX” is the random 3 letter name set as
this specific modification’s residue name. Additional temporary outputs include topology, coordinates, and mol2 or pdb files from each step of the stitching together of the modules, which is useful to ensure the correct modification is being added to the correct location. More information about Amber file formats can be found at the following link: [Amber File Formats](https://ambermd.org/FileFormats.php)

In addition to charge assignment for atoms in each module, the other parameters of the nucleic acid force field must be assigned. Bonds, angles, and dihedrals are assigned
based on similarity to current parameters or through the General Amber Force Field (GAFF) [13]. Van der Waals radii are similarly assigned through comparison. A flowchart
describing a general procedure for using modXNA is shown in Figure 4.

![](C:\Users\cnb11\Desktop\Work\siRNA\figures\ModXNA%20general%20pipeline.png)**Figure 4. ** ModXNA overall pipeline flowchart. Part 1 - Users will need to define their residue by comparison to the modXNA catalog. Part 2 - Users will need to build their modified nucleotides with modXNA.sh. Part 3 - To be done concurrently with Part 2, users should visually inspect the modified nucleotide. Part 4- Use the library files built in Part 2 along with frcmod.modxna to build a system parameter and coordinate
file pair using LEaP. 

## 3. Building modified oligonucleotides with modXNA

### 3.1 Incorporating pseudouridine in an RNA duplex

A simple initial example is running MD simulations of an RNA duplex that contains a pseudo-uridine (PSU) modification. We will base this example on the Nuclear Magnetic
Resonance (NMR) determined structure deposited with PDB ID 6I1W [14] that contains the PSU residue at position 5 of the sense strand. After downloading the legacy PDB Format file (6I1W.pdb) some cleaning is required to have a working PDB file that we can then use for building parameter and input coordinate files for MD simulations. There are multiple ways to do this. In this case, we will use the prepareforleap command available in CPPTRAJ [15]. Briefly, this command removes information that the PDB file includes that is either not needed by or would cause issues with LEaP, removes
water molecules, and performs a brief check of the residues to report problems that LEaP might encounter. In a text file (prepareforleap.cpptraj), add the commands (Listing 1):

> Listing 1. prepareforleap.cpptraj
> 
> parm 6I1W.pdb
> loadcrd 6I1W.pdb name edit
> prepareforleap crdset edit name FromPrepForLeap \
> out FromPrepForLeap−tleap.in leapunitname x \
> pdbout FromPrepForLeap.pdb nowat noh
> go

we execute the command using:

```
$ cpptraj -i prepareforleap.cpptraj
```

CPPTRAJ will open the 6I1W.pdb file, delete the water molecules, strip hydrogens and generate a new PDB file (from-prepareforleap.pdb). The CPPTRAJ output will print a
warning referring to the PSU residue:

*Potential problem : PSU_5_A is an unrecognized name and may not have parameters .*

This is because PSU is not a residue with parameters in the standard AMBER force fields. We will now use modXNA to create the parameters for the PSU residue. Even though modXNA is available in AmberTools25, it is recommended that the user download the latest version of modXNA from the website (https://modxna.chpc.utah.edu/ or this article’s github repository (https://github.com/ManghraniA/ModXNALiveComms). As commented in the introduction, modXNA requires an input file that declares the backbone, the sugar, and the base. Looking at the catalog of available fragments on the modXNA website, we need the RNA phosphate backbone (modXNA code RPO), a ribose sugar (modXNA code RC3), and the pseudo uridine base (modXNA code PUU). The input file (PSU.in.modxna, Listing 2) is:

> Listing 2. PSU.in.modxna
> 
> \##modXNA input file
> 
> RPO RC3 PUU

We run the script with:

`$ modXNA.sh -i PSU.in.modxna -m PSU`

The -i flag indicates the input file, and the -m flag indicates the name that we want for the generated nucleotide. Without the -m flag, modXNA will generate a random 3-letter
residue name. The script will generate multiple temporary files (with the prefix “tmp.”) that are useful to troubleshoot if something went wrong. If the script is completed successfully, we will have a PSU.lib file with the nucleotide, created from the selected fragments in the modXNA input script. We can check the structure of the generated nucleotide using the tmp.opt.pdb file in any molecular visualizer (Figure 5).

![](C:\Users\cnb11\Desktop\Work\siRNA\figures\figure4.png)

**Figure 5. ** Visualization of the tmp.opt.pdb structure for new residue name PSU.

We can now build our Amber topology and coordinates files using LEaP. We use the following tleap script (build.tleap, Listing 3):

> Listing 3. build.tleap
> 
> \###########################
> source leaprc.DNA.OL15
> source leaprc.RNA.OL3
> source leaprc.water.opc
> \##########################
> 
> \## Load the LIB file created with modXNA
> 
> loadoff PSU.lib
> 
> \## Load the force mod file included in modXNA
> 
> loadamberparams /home/user/modXNA/dat/frcmod.modxna
> \######################################
> x = loadpdb FromPrepForLeap.pdb
> solvateoct x OPCBOX 9.0
> addionsrand x Na+ 0
> addionsrand x Cl− 0
> saveamberparm x complex.parm7 complex.crd
> charge x
> quit

We can run the above script with the following comand:

`$ tleap -s -f build.tleap`

LEaP will load the declared AMBER force fields, load the LIB file and the frcmod files, then load the PDB structure generated from the prepareforleap command. The files complex.parm7 and complex.crd will be successfully generated. Note that atom names for residues built by modXNA may not match your particular PDB file. In this case, it is usually acceptable to keep commonly named atoms, delete atoms with naming differences, and have LEaP build in the remainder of the nucleotide based on the coordinates and connectivity in the lib file.

### 3.2 Incorporating 2'-O-methylated (O-methyl) nucleotides in an RNA duplex

In the second example we will build the PDB 1NAO [16], which consists of a DNA-RNA hybrid where the DNA contains 2’-O-methyl ribose sugar nucleotides in the anti-sense
strand. As before, we download the PDB file and use a similar CPPTRAJ script to obtain a pruned ‘clean’ PDB (prepareforleap.cpptraj, Listing 4):

> Listing 4. prepareforleap.cpptraj
> 
> parm 1nao.pdb
> loadcrd 1nao.pdb name edit
> prepareforleap crdset edit name FromPrepForLeap \
> out FromPrepForLeap−tleap.in leapunitname x \
> pdbout FromPrepForLeap.pdb nowat noh
> go

Run the script with:

```
$ cpptraj -i prepareforleap.cpptraj
```

and we get several warnings:

*Potential problem: OMG_10_B is an unrecognized name and may not have parameters.
Potential problem: OMU_11_B is an unrecognized name and may not have parameters.
Potential problem: OMC_12_B is an unrecognized name and may not have parameters.
Potential problem: OMC_17_B is an unrecognized name and may not have parameters.
Potential problem: OMC_18_B is an unrecognized name and may not have parameters.*

These are the C2’-O-methyl modified residues that are not included in the AMBER force fields for standard nucleotides. We will use modXNA to create the following new residues, outlined in Table 1: 5MG, OMU, OMC, 3MC.



**Table 1.** Description of input file modules for building four modified residues in PDBID 1NAO.

| Residue Name | Backbone fragment code | Sugar fragment code | Base fragment code | Notes       |
| ------------ |:----------------------:|:-------------------:|:------------------:| ----------- |
| 5MG          | 5PO                    | OME                 | DGG                | 5'-terminal |
| OMU          | DPO                    | OME                 | RUU                |             |
| OMC          | DPO                    | OME                 | DCC                |             |
| 3MC          | DPO                    | OME                 | DCC                | 3'-terminal |

It is important to note that 5’ and 3’ terminal nucleotides are unique from the central nucleic acid fragments that are expected to have both a 5’ and 3’ bond. For this system, we require a 5’-terminal modified residue (guanine with C2’-Omethyl sugar) and a 3’-terminal modified residue (cytosine with C2’-O-methyl sugar). Terminal nucleotide parameterization is only possible with modXNA version 1.8 or later, and the only termini parameterized are hydroxyl capping groups. We will first focus on the central residues: a C2’-O-methyl uridine and a C2’-O-methyl cytosine are required to match the
sequence from the experimental structure. 

Create a new directory with the name OMU. Inside the directory, create the modXNA input file (OMU.in.modxna, Listing
5):

> Listing 5. OMU.in.modxna
> 
> \### C3'-OMe-uridine
> 
> DPO OME RUU

Then, run modXNA.sh with:

```
$ modXNA.sh -i OMU.in.modxna -m OMU
```

Create a new directory with the name OMC. Inside the directory, create the modXNA input file (OMC.in.modxna, Listing 6):

> Listing 6. OMC.in.modxna
> 
> \## C3'-OMe-cytosine
> 
> DPO OME DCC

Then, run the modXNA.sh scrip tinside the OMD directory:

```
$ modXNA.sh -i OMC.in.modxna -m OMC
```

We can open the files OMU/tmp.opt.pdb and OMC/tmp.opt.pdb and confirm that both nucleotides have been created correctly (Figure 6).

![](C:\Users\cnb11\Desktop\Work\siRNA\figures\figure5.png)**Figure 6.**Left: OMU, right: OMC. Both structures are named tmp.opt.pdb, created using modXNA.

These residues are positions OMU11, OMC12 and OMC17 of the anti-sense strand (B chain) of the original 1NAO structure. To generate the 5’-terminal OMG residue, we can use the modXNA script with the option --5cap. We are also going to use the 5PO backbone module to specify the 5’ cap. We are going to name the 5’-terminal residue 5MG, to not confuse with a central fragment. Create a new folder with the name 5MG, as well as the modXNA input file (5MG.in.modxna, Listing 7):

> Listing 7. 5MG.in.modxna
> 
> \### 5'-C2'-OMe-guanine
> 
> 5PO OME DGG

Run the modXNA.sh script inside the 5MG directory:

```
$ modXNA.sh -i 5MG.in.modxna -m 5MG --5cap
```

Notice the extra flag --5cap and the 5PO backbone module are both specified to indicate to modXNA that we want a 5’-terminal residue. A 5’ residue must have a
companion 3’ residue for the charge on the entire system to be an integer. For the 3’-terminal residue, create a new folder with the name 3MC and create the modXNA input file (3MC.in.modxna, Listing 8):

> Listing 8. 3MC.in.modxna
> 
> \### C2'-OMe-cytosine
> 
> DPO OME DCC

And run modXNA:

```
$ modXNA.sh -i 3MC.in.modxna -m 3MC --3cap
```

Notice the extra flag --3cap to indicate to modXNA that we want a 3’-terminal residue. The resulting 5’-terminal and 3’-terminal are shown in Figure 7. Note that a 5’terminal requires a 3’terminal counterpart, otherwise the charge on the molecule will not be an integer. We now have all the required library files (5MG.lib, OMU.lib, OMC.lib and 3MC.lib) to build the RNA-DNA hybrid that matches the structure.

![](C:\Users\cnb11\Desktop\Work\siRNA\figures\figure6-new.png)**Figure 7.**Left: 5' terminal 5MG residue. Right: 3' terminal 3MC residue. Both are the tmp.opt.pdb files created by modXNA.

To build the RNA-DNA hybrid 1NAO, we will first combine the different lib files into a single file. To do this, create the text file create-1nao-lib.tleap (Listing 9):

> Listing 9. create-1nao-lib.tleap
> 
> loadoff 3MC/3MC.lib
> saveoff 3MC 1NAO.lib
> loadoff 5MG/5MG.lib
> saveoff 5MG 1NAO.lib
> loadoff OMC/OMC.lib
> saveoff OMC 1NAO.lib
> loadoff OMU/OMU.lib
> saveoff OMU 1NAO.lib
> quit

Run the file with:

<div>
$ tleap -s -f create-1nao-lib.tleap
</div>

The second thing we must do is edit the FromPrepFor-Leap.pdb file to denote 5’ and 3’ terminal residues. Copy FromPrepForLeap.pdb to EditTerminiFromPrepForLeap.pdb
and change the following residue names: OMG 10 to 5MG, OMC 18 to 3MC. 

Once these two steps are completed, we can create input for tleap and build the parmtop and crd files, as shown below in Listing 10:

> Listing 10. Build1nao.tleap
> 
> \###########################
> source leaprc.DNA.OL15
> source leaprc.RNA.OL3
> source leaprc.water.opc
> \###########################
> loadoff 1NAO.lib
> loadamberparams /home/user/modXNA/dat/frcmod.modxna
> \###########################
> x = loadpdb EditTerminiFromPrepForLeap.pdb
> solvateoct x OPCBOX 9.0
> addionsrand x Na+ 0
> saveamberparm x 1nao.parm7 1nao.rst7
> quit

### 3.3 Building a new modified sequence

In this next example, we will use modXNA to build a single strand antisense oligonucleotide (ASO) that complements with an RNA target. The visual representation of the ASO is depicted in Figure 8.

![](C:\Users\cnb11\Desktop\Work\siRNA\figures\figure7.jpg)

**Figure 8.**Sketch representation of the ASO. Sugar modifications: Blue is cET (k), red is C2’-O-methyl (m) and orange is C2’-Omethoxyethyl (e). Backbone modifications: red circle between the bases is phosphate linkage and grey circle with an ‘s’ is Phosphorothioate linkage. The <sup>m</sup>C is 5 methyl cytosine.

ASOs are typically heavily modified and contain a hybrid DNA-RNA-DNA strand. In the example here, the modifications include a 5 methyl cytosine base (M5C) modification, a
gapmer (CET) sugar modification, a 2’-O-methoxyethyl (MOE) sugar modification, a 2’-O-methyl (OME) sugar modification, and a phosphorothioate (PS1) backbone modification. To build the ASO, we need the above modifications combined into the following new nucleotides: 5CG, CEG, MEC, MOG, MOA, POG, POA, MSC, OMU, 3EC, detailed in Table 2.

**Table 2.** Description of input file modules for building ten modified residues in model ASO.

| Residue name | Backbone fragment code | Sugar fragment code | Base fragment code | Details of fragments                                    | Notes           |
| ------------ | ---------------------- | ------------------- | ------------------ | ------------------------------------------------------- | --------------- |
| 5CG          | 5PO                    | CET                 | RGG                | 5' hydroxyl cap backbone + cET + guanine                | 5'-terminal     |
| CEG          | RPO                    | CET                 | RGG                | Phosphate backbone + cET + guanine                      | Central residue |
| MEC          | RPO                    | CET                 | M5C                | Phosphate backbone + cET + 5 methyl cytosine            | Central residue |
| MOG          | RPO                    | MOE                 | RGG                | Phosphate backbone + MOE + guanine                      | Central residue |
| MOA          | RPO                    | MOE                 | RAA                | Phosphate backbone + MOE + adenine                      | Central residue |
| POG          | PS1                    | MOE                 | RGG                | Phosphorothioate (S) backbone + MOE + guanine           | Central residue |
| POA          | PS1                    | MOE                 | RAA                | Phosphorothioate (S) backbone + MOE + adenine           | Central residue |
| MSC          | PS1                    | CET                 | M5C                | Phosphorothioate (S) backbone + cET + 5 methyl cytosine | Central residue |
| OMU          | PS1                    | MOE                 | RUU                | Phosphorothioate (S) backbone + O-methyl + uracil       | Central residue |
| 3EC          | PS1                    | CET                 | M5C                | Phosphorothioate (S) backbone + cET + 5 methyl
cytosine | 3'-terminal     |

All the new residues follow the same procedure:

1. Make a new folder (to keep things organized) with the new residue name (i.e. 5CG)

2. In the new folder, make a modXNA input file with the fragment code of the backbone, sugar and base (5CG.in.modxna with the codes: 5PO CET DGG, for the
   residue 5CG)

3. Run modXNA (remember to use the --5cap for the 5CG residue and --3cap flag for the 3EC residue): \`$ modXNA .sh -i 5CG.in. modxna -m 5CG --5cap`

4. Check that the new residues have been built successfully by opening the tmp.opt.pdb structure.

The resulting structures are shown in Figure 9. 

![](C:\Users\cnb11\Desktop\Work\siRNA\figures\figure8.png)**Figure 9.**Visualizations of all the tmp.opt.pdb files generated by modXNA for each of the new nucleotides.



Before we can use these 10 library files, we need to combine all the different LIB files into a single file. This will make it easier to load just one bigger LIB file into LEaP instead of loading 10 independent LIB files. To do this, create the text file (create-ASO-lib.tleap, Listing 11):

> Listing 11. create-ASO-lib.tleap

> loadoff 3EC/3EC.lib
> saveoff 3EC ASO.lib
> loadoff 5CG/5CG.lib
> saveoff 5CG ASO.lib
> 
> loadoff CEG/CEG.lib
> saveoff CEG ASO.lib
> loadoff MEC/MEC.lib
> saveoff MEC ASO.lib
> loadoff MOA/MOA.lib
> saveoff MOA ASO.lib
> loadoff MOG/MOG.lib
> saveoff MOG ASO.lib
> loadoff MSC/MSC.lib
> saveoff MSC ASO.lib
> loadoff OMU/OMU.lib
> saveoff OMU ASO.lib
> loadoff POA/POA.lib
> saveoff POA ASO.lib
> loadoff POG/POG.lib
> saveoff POG ASO.lib
> quit

Run the file with:

`$ tleap -s -f create-ASO-lib.tleap`

This will load all the LIB files we created with modXNA and save the information into the same ‘ASO.lib’ file. If you inspect the ASO.lib file generated by tleap, you will see that all the new residues are in one single file. 

Now that we have all the required nucleotides, we need to build the ASO. There are several methodologies that we can use to proceed. The simplest one is generating a single strand with the correct sequence of nucleotide analogs using the ‘sequence’ command available in LEaP, using the following script (Aso-single-strand.tleap, Listing 12):

> Listing 12. Aso-single-strand.tleap
> 
> \###########################
> source leaprc.DNA.OL24
> source leaprc.RNA.OL3
> source leaprc.water.opc
> \###########################
> loadoff ASO.lib
> loadamberparams /home/user/modXNA/dat/frcmod.modxna
> \###########################
> ASO=sequence{5CG CEG MEC MOG MOA MOG POG POA MSC OMU 3EC}
> solvateoct ASO OPCBOX 9.0
> addionsrand ASO Na+ 0
> addionsrand ASO Cl− 0
> saveamberparm ASO ASO−ss.parm7 ASO−ss.crd
> quit

LEaP will create a linearized structure with the correct sequence of residues. It will have some structural overlaps, however, and if we run a quick structural minimization using
AMBER/sander (using the 9-step protocol from Amber-MDPrep [15]), we obtain a structure with no clashes, shown in Figure 10.

![](C:\Users\cnb11\Desktop\Work\siRNA\figures\figure9.png)**Figure 10.**Top: PDB structure generated directly from the ASO-single-strand.parm7 and ASO-single-strand.crd. Bottom: structure obtained after running minimization/equilibration using AmberMD-Prep.

## 4. Conclusions

In the above tutorial, we have outlined the use of modXNA.sh, a script and database of pre-parameterized modules, in order to create bespoke modified nucleotides. We have shown examples of building structures pulled from the PDB, which contain modified residues (1NAO and 6I1W), as well as using the program to create a heavily modified single stranded RNA ASO. We have outlined the procedure as well as identified some common pitfalls. While we know that building new residues for MD simulations can be difficult, we hope this work streamlines the process and decreases the time- to simulating for users. Users should consider this tutorial as a starting point for simulations, and proceed to validate the parameters from modXNA in a way that best addresses the problem they are investigating, using these citations [16–20] as examples.

## 5. Content and Links

All scripts used in this tutorial as well as the modXNA script are available at:
https://github.com/ManghraniA/ModXNALiveComs


modXNA is available at:
https://modxna.chpc.utah.edu


AmberTools can be installed through conda or by downloading the source code and compiling with CMake. Instructions for installing AmberTools are available at:
https://ambermd.org/Installation.php


For questions, to report an issue, or to recommend a newfragment, please contact us at: modxnainfo@gmail.com

## Author Contributions

CB and RGM wrote the manuscript; RGM performed simulations; DRR, PDD, TEC3, OL, and AM contributed to proofreading and testing; AM formatted for publication.


## Potentially Conflicting Interests

 The authors declare no competing financial interest.

## Funding Information

 The authors gratefully acknowledge the computing and data resources from the Center for High-Performance Computing at the University of Utah and the University of Maryland at IBBR. Funding for this research came from the National Institutes of Health R-01 GM-081411. D.R.R. was funded through NHLBI. P.D.D. was funded through CSIC I+D 22520220100376UD.

## Disclaimer

This article was funded by the National Institute of Standards and Technology and is not subject to U.S. Copyright. Certain commercial equipment, instruments, software, or materials are identified in this paper to foster understanding. Such identification does not imply recommendation or endorsement by the National Institute of Standards and Technology, nor does it imply that the materials or equipment identified are necessarily the best available for the purpose.

## Author Information

### ORCID:

 Rodrigo Galindo-Murillo: orcid.org/0000-0001-5847-4143 

Akanksha Manghrani: orcid.org/0000-0002-6498-7305 

Daniel R. Roe: orcid.org/0000-0002-5834-2447 

Olivia Love: orcid.org/0000-0002-7743-8495 

Pablo D. Dans: orcid.org/0000-0002-5927-372X 

Thomas E. Cheatham III: orcid.org/0000-0003-0298-3904 

Christina Bergonzo: orcid.org/0000-0003-1990-2912



## References

[1] Love O, Galindo-Murillo R, Roe DR, Dans PD, Cheatham III TE, Bergonzo C. modXNA: A Modular Approach to Parametrization of Modified Nucleic Acids for Use with Amber Force Fields. Journal of Chemical Theory and Computation. 2024; 20:9354–9363.

https://doi.org/10.1021/acs.jctc.4c01164


[2] Weiner SJ, Kollman PA, Nguyen DT, Case DA. An all atom force field for simulations of proteins and nucleic acids. Journal of Computational Chemistry. 1986; 7(2):230–252.

https://doi.org/10.1002/jcc.540070216


[3] Cornell WD, Cieplak P, Bayly CI, Kollman PA. Application of RESP Charges To Calculate Conformational Energies, Hydrogen Bond Energies, and Free Energies of Solvation. Journal of the American Chemical Society. 1993; 115:9620–9631.

https://doi.org/10.1021/ja00074a030


[4] Cieplak P, Cornell WD, Bayly C, Kollman PA. Application of the multimolecule and multiconformational RESP methodology to biopolymers: Charge derivation for DNA, RNA, and proteins. Journal of Computational Chemistry. 1995; 16:1357–1377.

https://doi.org/10.1002/jcc.540161106


[5] Bayly CI, Cieplak P, Cornell WD, Kollman PA. A well behaved electrostatic potential based method using charge restraints for deriving atomic charges: The RESP model.
Journal of Physical Chemistry. 1993; 97:10269–10280.

https://doi.org/10.1021/j100142a004

[6] Cheatham III TE, Cieplak P, Kollman PA. A modified version of the Cornell et al. force field with improved sugar pucker phases and helical repeat. Journal of Biomolecular Structure and Dynamics. 1999; 16:845–862.

https://doi.org/10.1080/07391102.1999.10508297


[7] Pérez A, Marchán I, Svozil D, Šponer J, Cheatham III TE, Laughton CA, Orozco M. Refinement of the AMBER force field for nucleic acids: Improving the description of α/γ
conformers. Biophysical Journal. 2007; 92:3817–3829.

https://doi.org/10.1529/biophysj.106.097782


[8] Zgarbová M, Otyepka M, Šponer J, Mladek A, Banas P, Cheatham III TE, Jurečka P. Refinement of the Cornell et al. nucleic acids force field based on reference quantum
chemical calculations of glycosidic torsion profiles. Journal of Chemical Theory and Computation. 2011; 7:2886–2902.

https://doi.org/10.1021/ct200162x


[9] Ivani I, Dans PD, Noy A, Perez A, Faustino I, Hospital A, Walther J, Andrio P, Goni R, Balaceanu A, Portella G, Battistini F, Gelpi JL, Gonzalez C, Vendruscolo M, Laughton CA,
Harris SA, Case DA, Orozco M. Parmbsc1: A refined force field for DNA simulations. Nature Methods. 2015; 13:55–58.

https://doi.org/10.1038/nmeth.3658


[10] Zgarbová M, Šponer J, Jurečka P. Z-DNA as a Touchstone for Additive Empirical Force Fields and a Refinement of the Alpha/Gamma DNA Torsions for AMBER. Journal
of Chemical Theory and Computation. 2021; 17:6292–6301.

https://doi.org/10.1021/acs.jctc.1c00697


[11] Roe DR, Cheatham III TE. PTRAJ and CPPTRAJ: Software for processing and analysis of molecular dynamics trajectory data. Journal of Chemical Theory and Computation. 2013; 9(7):3084–3095.

https://doi.org/10.1021/ct400341p


[12] Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA. Development and testing of a general Amber force field. Journal of Computational Chemistry. 2004; 25:1157–1174.

https://doi.org/10.1002/jcc.20035


[13] Deb I, Popenda L, Sarzynska J, Lahiri A, Walczak A, Biala E, Gdaniec Z, Kierzek R, Banerjee R. Computational and NMR studies of RNA duplexes with an internal pseudouridine-adenosine base pair. Scientific Reports. 2019; 9. https://doi.org/10.1038/s41598-019-52637-0


[14] Roe DR, Bergonzo C. prepareforleap: An automated tool for fast PDB-to-parameter generation. Journal of Computational Chemistry. 2022; 43(13):930–935. https://doi.org/10.1002/jcc.26847


[15] Roe DR, Brooks BR. A protocol for preparing explicitly solvated systems for stable molecular dynamics simulations. Journal of Chemical Physics. 2020; 153(5):054123.

https://doi.org/10.1063/5.0013849


[16] Linzer JT, Aminov E, Abdullah AS, Kirkup CE, Diaz Ventura RI, Bijoor VR, Jung J, Huang S, Tse CG, Álvarez Toucet E, Onghai HP, Ghosh AP, Grodzki AC, Haines ER, Iyer AS, Khalil MK, Leong AP, Neuhaus MA, Park J, Shahid A, et al. Accurately Modeling
RNA Stem-Loops in an Implicit Solvent Environment. Journal of Chemical Information and Modeling. 2024; 64(15):6092–6104.

https://doi.org/10.1021/acs.jcim.4C00756

[17] Bergonzo C, Grishaev A. Critical Assessment of RNA and DNA Structure Predictions via Artificial Intelligence: the Imitation Game. Journal of Chemical Information and Modeling. 2025; 65:3544–3554. 

https://doi.org/10.1021/acs.jcim.5c00245


[18] Winkler L, Cheatham III TE. Benchmarking the Drude Polarizable Force Field Using the r(GACC) Tetranucleotide. Journal of Chemical Information and Modeling. 2023; 63(8):2505–2511. 

https://doi.org/10.1021/acs.jcim.3c00250


[19] Love O, Galindo-Murillo R, Zgarbová M, Šponer J, Jurečka P, Cheatham III TE. Assessing the Current State of Amber Force Field Modifications for DNA–2023 Edition. Journal of Chemical Theory and Computation. 2023; 19:4299–4307.

https://doi.org/10.1021/acs.jctc.3c00233


[20] Bergonzo C, Grishaev A. Maximizing accuracy of RNA structure in refinement against residual dipolar couplings. Journal of Biomolecular NMR. 2019; 73:117–139.

https://doi.org/10.1007/s10858-019-00236-6







