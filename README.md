# Modelling and sampling protease substrates

## Package to model peptide substrates bound to annotated protease structures and run promiscuity analysis

* From publication *"Modelling peptide substrates bound to proteases: insights to study binding promiscuity"*
* BMC Bioinformatics, 2020
* Authors: Rodrigo Ochoa, Mikhail Magnitov, Roman A. Laskowski, Pilar Cossio, Janet M. Thornton

## Purpose

The goal of these scripts is to provide tools for modelling natural substrates bound to proteases with structural information available. The first script is focused on the modelling, which includes protocols for the modification of non-natural amino acids by natural counterparts, and the modelling of missing regions within the 8-mer peptide region that binds the protease binding site. In addition, a second script allows the model of any peptide based on the templates obtained from the first scripts, in order to run dynamic analysis and calculate averages of the structures observables reported in the paper, including the accessible surface area (ASA) and the interface interaction energy.

## Third-party tools required:

- BioPython: https://biopython.org/wiki/Download
- RDKit: https://github.com/rdkit/rdkit/releases
- Modeller: https://salilab.org/modeller/download_installation.html
- DSSP: https://github.com/cmbi/hssp/archive/3.1.4.tar.gz
- Rosetta Commons: https://www.rosettacommons.org/software/license-and-download

The BioPython and RDKit modules can be installed directly from package repositories. Modeller can be installed freely after obtaining an academic license. DSSP require to be compiled using the source code of the latest version. For the Rosetta functionalities the recommended is to follow the installation instructions, and take note of the Rosetta version that will be provided in the script.

## Extra files and databases required:

To run the modelling script, a file and a database are required to perform the analysis. These are:

1. components.cif (file with information of compounds and modified ligands in the PDB). The file can be downloaded from: ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif (272 MB)

**NOTE:** Put the file into the auxiliar folder of the code project

2. MEROPS mySQL database (mySQL schema containing the information required to identify enzyme substrates). The database can be downloaded from: ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/meropsweb121.tar.gz (727 MB)

**NOTE:** To install the mysql database you can follow these instructions:

- Enter your mySQL installed engine and run the command:

`mysql> CREATE DATABASE merops;`

- Then from the command line, import the schema as:

`mysql -u username -p merops < merops.sql`
You can put here your user name and the password after prompted in the terminal.
**The user, password and MEROPS database name should be provided in the script to run the analysis succesfully**

## How to run the modelling script

The basic command line to run the script is:

`model_substrate_protease.py [-h] [-s STRUCTURE] [-f FAMILY]
                                   [-n MAX_SUBSTRATES] [-m MODE]
                                   [-t SIM_THRESHOLD] [-r ROSETTA] -u USER -p
                                   PASSWORD -d DATABASE`
                                       
where the arguments are:

```
optional arguments:
  -h, --help         show this help message and exit
  -s STRUCTURE       Add a protease structure PDB id to see if it is possible
                     to model substrates
  -f FAMILY          Select the family of proteases used as reference to model
                     the substrates from 1) serine, 2) cysteine
  -n MAX_SUBSTRATES  Define the maximum number of substrates that can be
                     modelled per structure available
  -m MODE            Choose a mode to run the script from two options: 1)
                     ready, 2) complete. Ready are the peptide templates with
                     natural amino acids, and Complete are the peptides
                     containing at least one amino acid
  -t SIM_THRESHOLD   Similarity threshold to decide which amino acid will
                     replace the NNAA present in the template
  -r ROSETTA         Version of Rosetta that will be implemented
  -u USER            User to access the local installation of the MEROPS
                     database
  -p PASSWORD        Password to access the local installation of the MEROPS
                     database
  -d DATABASE        Name of the local installation of the MEROPS database
 ```
 
The required arguments are the information required to access the mySQL MEROP database. For the other parameters the script has default values, including a PDB of reference in case is part of the annotated proteases, or just run the models for all the families and structures present in the dataset. **However, please be aware of changing the Rosetta version through the flag or directly in the script by default.**

## Examples

### Model a peptide based on a substrate composed of natural amino acids

The modelling of substrates bound to proteases has two main modes of running. One involves the selection of proteases bound to peptides composed of natural amino acids (called *ready* in the script). The second mode is when the substrate has one non-natural amino acids (NNAA) to modify. The following is a command example for the first case:

`python3.5 model_substrate_protease.py -f serine -n 1 -m ready -u user -p password -d merops`

In this case, we are allowing the modelling of only 1 substrate per structure available in the dataset bound to peptides with natural amino acids. For this script two families are available: serine and cysteine proteases. Here we selected serine as the reference family. Finally the MEROPS credentials are provided to look for reported substrates

After that, the models will be stored in **models/model_ready** with the name `[structure]\_[substrate sequence].pdb`, where the structure is the PDB id and the substrate sequence is the peptide that was fully modelled using the protocol. A report of the model is provided in the following form:
```
pdb,pepTemplate,pepModel,chain,merops,old_aa,new_aa,pos_aa,uniprot
1smf,-CAKSI--,SCAKSIIG,I,S01.151,THR,ALA,3,O60256
...
```
Here the fragment CAKSI was modelled into the 8-mer peptide SCAKSIIG, which is part of the substrate protein with UniProtKB id O60256. The model was performed based on the protease structure with PDB is 1smf that belongs to the MEROPS family S01.151 (trypsin). In addition, an amino acid from the original structure had to be changed by another one reported in the substrate, in this case a threonine by an alanine in the position 3 of the peptide.

### Model a peptide based on a substrate containing one non-natural amino acid (NNAA)

For the second case, we require the modification of a NNAA for a natural amino acid based on a similarity threshold defined by the user (called *complete* in the script). It means that we can model a substrate depending on how similar we want the template and new amino acids should be. The following is an example based on a particular PDB structure available in the dataset:

`python3.5 model_substrate_protease.py -s 1tps -f serine -n 1 -m complete -t 0.4 -u user -p password -d merops`

Here we are modelling a substrate based on the template present in the structure with PDB is 1tps. The natural amino acid that will replace the present NNAA require to be at least 40\% similar or more based on Tanimoto comparison of the molecular representations of the amino acid side chains. We selected serine as the reference family, and the MEROPS credentials are also provided to look for reported substrates.

After that, the models will be stored in **models/model_complete** with the name `[structure]\_[substrate sequence]\_modelled.pdb`, where the structure is the PDB id and the substrate sequence is the peptide that was fully modelled using the protocol. A report of the model is provided in the following form:

```
pdb,pepTemplate,pepModel,chain,merops,old_aa,new_aa,pos_aa,uniprot,sim
1tps,-LTREL--,FLTRELAE,B,S01.151,DLE,LEU,247,P23396,1.0
```

Here the fragment LTREL was modelled into the 8-mer peptide FLTRELAE, which is part of the substrate protein with UniProtKB id P23396. The model was performed based on the protease structure with PDB is 1tps that belongs to the MEROPS family S01.151 (trypsin). In addition, a NNAA from the original structure had to be changed by another one reported in the substrate, in this case the residue DLE by an L-Leucine, with a similarity of 100\%.

### How to run the sampling script

After having modelled 8-mer peptides in reference structures, it is possible to call the second script for modelling any peptide of interest, run a simulation of the system using the backrub method from Rosetta and calculate structural descriptors from the trajectory. The basic command line to run the script is:

`run_dynamic_proteases.py [-h] -p PATH -s SEQUENCE -c CHAIN [-r ROSETTA]`
                                       
where the arguments are:

```
optional arguments:
  -h, --help   show this help message and exit
  -p PATH      Path of the structure that will be used to run the analysis
  -s SEQUENCE  Sequence of the peptide that will be modelled and sampled
  -c CHAIN     Chain identifier of the peptide in the structure
  -r ROSETTA   Version of Rosetta that will be implemented
 ```

The required arguments are the path of the model that we want to use as template, the sequence of the peptide that will be modelled, and the chain of the peptide in the structure of reference. **Please be aware of changing the Rosetta version through the flag or directly in the script by default.**

### Run the sampling of a modelled peptide using Rosetta

To run the dynamic analysis of a peptide of reference, we can call the script in the followin form:

`python run_dynamic_proteases.py -p models/model_ready/1smf_SCAKSIIG_modelled.pdb -s TGYHKLPR -c B`

Here we provide the path of the modelled structure that will be used as reference, the sequence TGYHKLPR of the new peptide, and the chain identifiier of the peptide, which is B in this case. The modelled peptide is subjected to Backrub for a defined number of Monte Carlo steps (this can be changed in the script), and the snapshots are used to calculate per peptide amino acid two structural observables: the accessible surface area using DSSP, and the interaction energy using an interface scoring of Rosetta. The new modelled peptide is stored in the route **dynamic/models**, and the observables per amino acid are stored in the routes **dynamic/observables/asa_[peptide sequence]** and **dynamic/observables/energy_[peptide sequence]**.

The averages per peptide are reported in the file `final_averages_[peptide sequence]`, with the following format for peptide TGYHKLPR:

```
Amino_acid Position Average_ASA Average_Energy
THR	1	0.72	-2.57
GLY	2	0.45	-2.17
TYR	3	0.28	-6.23
HIS	4	0.014	-9.59
LYS	5	0.31	-5.67
LEU	6	0.28	-4.31
PRO	7	0.52	-2.12
ARG	8	1.0	-1.17
```

The third column represents the average ASA for each amino acid in the peptide, and the four column the average energy. The script can be embedded in a loop to calculate the observable for a set of peptides of interest, as the example of the random libraries mentioned in the publication.

## Support

In case the protocol is useful for other research projects and require some advice, please contact us to the email: rodrigo.ochoa@udea.edu.co
