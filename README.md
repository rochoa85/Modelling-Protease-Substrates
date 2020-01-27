# Modelling and sample protease substrates

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
- Rosetta Commons: https://www.rosettacommons.org/software/license-and-download

The BioPython and RDKit modules can be installed directly from package repositories. Modeller can be installed freely after obtaining an academic license. For the Rosetta functionalities the recommended is to follow the installation instructions, and take note of the Rosetta version that will be provided in the script.

## Extra files and databases required:

To run the modelling script, a file and a database are required to perform the analysis. These are:

1. components.cif (file with information of compounds and modified ligands in the PDB)
The file can be downloaded from: ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif (272 MB)
**NOTE:** Put the file into the auxiliar folder of the code project

2. MEROPS mySQL database (mySQL schema containing the information required to identify enzyme substrates)
The database can be downloaded from: ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/meropsweb121.tar.gz (727 MB)
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

After that, the models will be stored in **models/model_ready** with the name [structure]\_[substrate sequence].pdb, where the structure is the PDB id and the substrate sequence is the peptides that was fully modelled using the protocol. A report of the model is provided in the following form:
```
pdb,pepTemplate,pepModel,chain,merops,old_aa,new_aa,pos_aa,uniprot
1smf,-CAKSI--,SCAKSIIG,I,S01.151,THR,ALA,3,O60256
...
```
Here the fragment CAKSI was modelled into the 8-mer peptide SCAKSIIG, which is part of the substrate protein with UniProtKB id O60256. The model was done on the protease structure with PDB is 1smf that belongs to the MEROPS family S01.151 (trypsin). In addition, an amino acid from the original structure had to be changes by another one reported in the substrate, in this case a threonine by an alanine in the position 3 of the peptide.

### Model a peptide based on a substrate containing one non-natural amino acid (NNAA)

To model peptides of interest in an allele with crystal structure available, the user can select from four alleles that have available crystals in the local folder provided. Based on the selection, the script reads the file `list_MHC_crystals.txt` and select the corresponding structure. One way to run the script is as follows:

`python3.5 model_peptides_MHC_complexes.py -l list_peptides.txt -m model -a 0301`

The file `<sequence>_modelled.pdb` will be generated. In addition, a report with the statistics of the modelling will be created. By default the name of the file is `stats_peptide_models.txt`- However, the name can be changed with the `-o` flag. An example of the report is the following:

```
Crystal	Allele	Peptide_template	Peptide_model	Core_template	Core_model	Alignment_score	Identities	Identical_pockets
1t5x	DRB1*01:01	AAYSDQATPLLLSPR	IPTAFSIGKTYKPEE	FSIGKTYKP	YSDQATPLL	0	2	1
The model was completed
```
In addition to the prediction of the cores, the peptides are aligned to evaluate how similar they are, and identities are mapped in the full sequence and in the core region. This is useful to understand how different is the sequence content of the peptides from the available peptide template.

### Run additional sampling using Rosetta Backrub

After modelling the peptides, it is possible to run a simulation of the system using the backrub method from Rosetta. For that purpose, just run the script with the following options:

`python3.5 model_peptides_MHC_complexes.py -l list_peptides.txt -m backrub -a 0101`

In that case, the model with the additional simulation will be run in the local computer. The output of the trajectory can be used to calculate any structural descriptor of interest using available packages to study interactions, or calculate other variables such as secondary structures or accessible surface areas, among others.

## Support

In case the protocol is useful for other research projects and require some advice, please contact us to the email: rodrigo.ochoa@udea.edu.co
