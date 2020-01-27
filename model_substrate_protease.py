#!/usr/bin/python3

"""
Package to model peptide substrates bound to annotated protease structures
NOTE: The protocol requires of auxiliary programs and Unix system commands - Tested on Ubuntu 16.04

From publication "Modelling peptide substrates bound to proteases: insights to study binding promiscuity"
Journal: BMC Bioinformatics
Authors: Rodrigo Ochoa, Mikhail Magnitov, Roman Laskowski, Pilar Cossio, Janet Thornton
Year: 2020

Third-party tools required:

BioPython: https://biopython.org/wiki/Download - Ubuntu package: python3-biopython
RDKit: https://github.com/rdkit/rdkit/releases
Modeller: https://salilab.org/modeller/download_installation.html
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__credits__ = ["Rodrigo Ochoa","Mikhail Magnitov","Roman Laskowski","Pilar Cossio","Janet Thornton"]
__license__ = "GPL"
__version__ = "1.0"
__email__ = "rodrigo.ochoa@udea.edu.co"

########################################################################################
# Modules to import
########################################################################################

import csv
import ast
import os
import subprocess
import argparse
import re
import MySQLdb as mdb

# Third-party modules: Biopython, RDKit and Modeller
from Bio.PDB import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import modeller
from modeller.automodel import *

########################################################################################
# Extra files and databases to install
########################################################################################

# 1. components.cif (file with information of compounds and modified ligands in the PDB)
# The file can be downloaded from: ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif (272 MB)
# NOTE: Put the file into the auxiliar folder of the code project

# 2. MEROPS mySQL database (mySQL schema containing the information required to identify enzyme substrates)
# The database can be downloaded from: ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/meropsweb121.tar.gz (727 MB)
# NOTE: To install the mysql database you can follow these instructions:
# - Enter your mySQL installed engine and run the command:
#
# mysql> CREATE DATABASE merops;
#
# - Then from the command line, import the schema as:
#
# mysql -u username -p merops < merops.sql (you can put here your user name and the password after prompted in the terminal)
#
# The user, password and MEROPS database name should be provided in the script to run the analysis succesfully

########################################################################################
# Functions
########################################################################################

def map_peptides(list_enzymes):
    """
    Read a csv file with the annotated proteases structures to generate dictionaries for the two categories included: peptides with natural amino acids and NNAAs
    
    Argument:
    list_enzymes -- csv file with the annotated structures
    
    Return:
    ready -- dictionary with information of peptide complexed to peptide fragments composed of natural amino acids
    model -- dictionary with information of peptide complexed to peptide fragments composed of max one non-natural amino acid (NNAA)
    """
    
    # List of amino acids
    aminoacids=["ALA","ASP","GLU","PHE","HIS","ILE","LYS","LEU","MET","GLY","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR","CYS"]
    
    # Start the dictionaries with the required information
    ready={}
    model={}
    
    # Read the CSV file containing the structures selected for each class
    with open(list_enzymes) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                line_count += 1
            else:
                # Store the basic information required from the csv file
                pdb=row[0]
                ligand_chain=row[3]
                merops=row[11]
                aa=row[15].split("-")
                flag_nAA=0
                
                # Check if any of the amino acids is non natural to activate a flag
                for a in aa:
                    if a not in aminoacids: flag_nAA=1
                
                # Save the PDBs with natural amino acids
                if flag_nAA==0:
                    pos=ast.literal_eval(row[17])
                    numbering=ast.literal_eval(row[16])
                    ready[pdb]=(aa,pos,merops,numbering,ligand_chain)
                        
                # If there are any NNAA, check that is only one
                if flag_nAA==1:
                    pos=ast.literal_eval(row[17])
                    numbering=ast.literal_eval(row[16])
                    model[pdb]=(aa,pos,merops,numbering,ligand_chain)    
                                        
                line_count += 1
    
    # Return the two dictionaries with the information requird for modelling
    return ready,model

##################################################################################

def count_merops_id(list_struct):
    """
    Count from the dictionary how many Merops Ids are in the list and count them
    
    Arguments:
    list_struct -- Dictionary with the information extracted from the csv file
    
    Return:
    total_mers -- Dictionary containing the counts of the Merops IDs
    """
    total_mers={}
    
    # Iterate over all the complexes
    for ele in list_struct:
        if list_struct[ele][2] not in total_mers: total_mers[list_struct[ele][2]]=1
        else: total_mers[list_struct[ele][2]]+=1
    
    # Return the dictionary with the counts
    return total_mers

##################################################################################

def modelling_fragment(pdb,pep_chain,pepT,pepM,rosetta_path,folder_path,new_aa=None):
    """
    Function to model the fragment of the desired peptide
    
    Arguments:
    pdb -- Code of the structure used to model the peptide
    pep_chain -- Chain of the peptide on the PDB
    pepT -- Template peptide sequence
    pepM -- Peptide sequence to model
    rosetta_path -- Path of the Rosetta installation in the local computer
    folder_path -- Path of the folder where the models will be saved
    new_aa -- In case the peptide was subjected to a mutation before, it will contain the modification. By default is None
    
    Return
    PDB structure of the peptide after growing flanking amino acids - Structure ready for promiscuity analysis
    """
    
    # Check if an amino acid was mutated to select the input structure accordingly
    if new_aa:
        os.system('cp models/{}/{}_{}_modelled.pdb .'.format(folder_path,pdb,new_aa))
        code = '{}_{}_modelled'.format(pdb,new_aa)
    else:
        os.system('cp models/{}/{}_relaxed.pdb .'.format(folder_path,pdb))
        code = '{}_relaxed'.format(pdb)
    
    # Start Modeller environment and write the sequence file 
    e = modeller.environ()
    m = modeller.model(e, file=code)
    aln = modeller.alignment(e)
    aln.append_model(m, align_codes=code)
    aln.write(file=code+'.seq')
    
    # Obtain the protein sequence with the peptide
    infoSeq=[x.strip() for x in open(code+'.seq')]
    header=[]
    sequenceLine=''
    for info in infoSeq:     
        if ">" not in info and ":" not in info:
            if info:
                sequenceLine+=info
        else:
            if info: header.append(info)
    
    # Split the sequence information as required by Modeller
    last_line=sequenceLine.split("/")
    sequenceTemp=last_line[0]+"/"+pepT+"*"
    sequenceMod=last_line[0]+"/"+pepM+"*"
    
    seqTempList = [sequenceTemp[i:i+75] for i in range(0, len(sequenceTemp), 75)]
    seqModList = [sequenceMod[i:i+75] for i in range(0, len(sequenceMod), 75)]
    
    # Create the alignment file
    alignmentFile=open("alignment.ali","w")
    for h in header: alignmentFile.write(h+"\n")
    for s in seqTempList: alignmentFile.write(s+"\n")
    alignmentFile.write("\n>P1;{}_fill\nsequence:::::::::\n".format(code))
    for s in seqModList: alignmentFile.write(s+"\n")
    alignmentFile.close()
    
    # Directories for input atom files
    e.io.atom_files_directory = ['.', '../atom_files']
    a = automodel(e, alnfile='alignment.ali', knowns='{}'.format(code), sequence='{}_fill'.format(code))
    a.starting_model= 1
    a.ending_model  = 1
    a.make()
    
    # Create resfile
    rosetta_config_file=open("resfile.config","wb")
    rosetta_config_file.write("NATRO\n")
    rosetta_config_file.write("start\n")
    rosetta_config_file.close()

    # Renumbering of the chain residues
    os.system("{}/main/source/bin/fixbb.default.linuxgccrelease -in:file:s {}_fill.B99990001.pdb -resfile resfile.config -renumber_pdb -per_chain_renumbering".format(rosetta_path,code))
    parser = PDBParser()
    structure = parser.get_structure('REF',"{}_fill.B99990001_0001.pdb".format(code))
    io = PDBIO()
    io.set_structure(structure)
    io.save("pre-post-modelled.pdb")
    
    # Relaxing of the final complex
    os.system("{}/main/source/bin/relax.default.linuxgccrelease -database {}/main/database -in:file:s pre-post-modelled.pdb -relax:thorough -relax:bb_move false".format(rosetta_path,rosetta_path))
    parser = PDBParser()
    structure = parser.get_structure('REF',"pre-post-modelled_0001.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save("post-modelled.pdb")
    
    # Deleter temporal files and save final model in the route of interest
    os.system("rm {}* score.sc resfile.config alignment.ali pre-post-modelled*".format(code))
    os.system("mv post-modelled.pdb models/{}/{}_{}_modelled.pdb".format(folder_path,pdb,pepM))
    
    # Delete temporal files
    if new_aa:
        os.system('rm models/{}/{}_{}_modelled.pdb .'.format(folder_path,pdb,new_aa))
        os.system('rm models/{}/{}_relaxed.pdb .'.format(folder_path,pdb))
    else:
        os.system('rm models/{}/{}_relaxed.pdb .'.format(folder_path,pdb))

##################################################################################

def replace_amino_acid_natural(pep_pdb,pep_chain,old_aa,new_aa,pep_position):
    """
    Left only backbone atoms for the replacement of a natural amino acid in a position of the pdb chain
        
    Arguments:
    pep_pdb -- Biopython object with the pdb structure of interest
    pep_chain -- Chain containing the sequence where the replacement will be done
    old_aa -- Amino acid that will be replaced - natural one
    new_aa -- Amino acid that will be mutated
    pep_position -- Position in the chain where the mutation will be done
    
    Output:
    A PDB structure called "pre-mutated.pdb" with only the backbone atoms
    """
    
    # Read the PDB file
    residues=pep_pdb.get_residues()
    chain=pep_pdb[0][pep_chain]
    
    # Report the mutation made
    message="The residue {} in chain {} and position {} will be changed by {}".format(old_aa,pep_chain,pep_position,new_aa)
    print(message)
    
    # Rename the residue
    chain[pep_position].resname=new_aa
    
    # Delete the other atoms leaving only the atoms of the backbone
    ids=[]
    for a in chain[pep_position]:
        atomId=a.id
        if atomId not in ("N","CA","O","C"): ids.append(atomId)
    for i in ids: chain[pep_position].detach_child(i)
        
    # Saving the new structure
    io = PDBIO()
    io.set_structure(pep_pdb)
    io.save("pre-mutated.pdb")
    
##################################################################################

def replace_amino_acid(pep_pdb,pep_chain,old_aa,new_aa,pep_position):
    """
    Left only backbone atoms for the replacement of a non-natural amino acid (NNAA) in a position of the pdb chain
        
    Arguments:
    pep_pdb -- Biopython object with the pdb structure of interest
    pep_chain -- Chain containing the sequence where the replacement will be done
    old_aa -- Amino acid that will be replaced - non-natural one (NNAA)
    new_aa -- Amino acid that will be mutated
    pep_position -- Position in the chain where the mutation will be done
    
    Output:
    A PDB structure called "pre-mutated.pdb" with only the backbone atoms
    """
    
    # Dictionaries with non-natural amino acids that will be modelled
    # NOTE: Here should be included new NNAA in case is required
    
    mod_amino_backbone={'DLE':("N","CA","O","C"),'RNG':("CH","CA","O","C"),'DPN':("N","CA","O","C"),
                        'AR7':("N","CA","O","CF","C"),'VAI':("N","CA","O","C"),'ASJ':("N","CA","O","C"),'LYK':("N","CA","O","C"),
                        'PHQ':("C2","O2","O1","C1"),'AKZ':("N","CA","O","C"),'ASA':("N","CA","O","C"),'4H0':("N","CA","CB","CG"),
                        'AA1':("N3","N6","O5","C4"),'MX3':("N20","N2","O6","C11"),'MX4':("N16","N9","O8","C7"),'MY1':("N4","N2","O3","C1"),
                        'MX5':("N20","N2","O6","C11"),'FPR':("N","CA","O","C"),'MY5':("N2","N1","O5","C4"),'1U8':("N","CA","O","C"),
                        'MY3':("N20","N2","O6","C11"),'MKE':("N","CA","O","C")}
    
    mod_amino_changes={'RNG':[("CH","N","N")],'AR7':[("CF","C","C")],'PHQ':[("C1","C","C"),("O1","O","O"),("O2","CA","C"),("C2","N","N")],
                       '4H0':[("CB","C","C"),("CG","O","O")],'AA1':[("N3","N","N"),("N6","CA","C"),("C4","C","C"),("O5","O","O")],
                       'MX3':[("N20","N","N"),("N2","CA","C"),("C11","C","C"),("O6","O","O")],'MX4':[("N16","N","N"),("N9","CA","C"),("C7","C","C"),("O8","O","O")],
                       'MY1':[("N4","N","N"),("N2","CA","C"),("C1","C","C"),("O3","O","O")],'MX5':[("N20","N","N"),("N2","CA","C"),("C11","C","C"),("O6","O","O")],
                       'MY5':[("N2","N","N"),("N1","CA","C"),("C4","C","C"),("O5","O","O")],'MY3':[("N20","N","N"),("N2","CA","C"),("C11","C","C"),("O6","O","O")]}
    
    # Read the PDB file
    chain=pep_pdb[0][pep_chain]
    for residue in chain:
        resPos=int(residue.get_full_id()[3][1])
        if resPos==pep_position:
            # Put the new name
            residue.resname=new_aa
            print(residue.get_full_id())
        
    # Delete the other atoms leaving only the atoms of the backbone
    ids=[]
    atom_numbers=[]
    atom_reference=0
    # Special code for manipulating non-natural amino acids in the PDB
    for c,a in enumerate(chain[("H_{}".format(old_aa),pep_position," ")]):
        atomId=a.id
        if c==0:
            atom_reference=a.serial_number
        if atomId not in mod_amino_backbone[old_aa]: ids.append(atomId)
    for i in ids: chain[("H_{}".format(old_aa),pep_position," ")].detach_child(i)
     
    # Append the corresponding atoms
    for a in chain[("H_{}".format(old_aa),pep_position," ")]:
        if old_aa in mod_amino_changes:
            for ch in mod_amino_changes[old_aa]:
                if ch[0]==a.id:
                    a.fullname=' {}'.format(ch[1])
                    a.name=' {}'.format(ch[1])
                    a.element='{}'.format(ch[2])
        atom_numbers.append(a.serial_number)
    
    
    # Saving the new structure
    io = PDBIO()
    io.set_structure(pep_pdb)
    io.save("pre-mutated.pdb")
    
    # Replace in the file the HETATM flags by ATOM records
    for i,num in enumerate(atom_numbers):
        if len(str(num))==4:
            os.system("sed -i 's#HETATM {atom}#ATOM   {atom}#g' pre-mutated.pdb".format(atom=atom_reference+i))
        if len(str(num))==3:
            os.system("sed -i 's#HETATM  {atom}#ATOM    {atom}#g' pre-mutated.pdb".format(atom=atom_reference+i))
        if len(str(num))==2:
            os.system("sed -i 's#HETATM   {atom}#ATOM     {atom}#g' pre-mutated.pdb".format(atom=atom_reference+i))
        if len(str(num))==1:
            os.system("sed -i 's#HETATM    {atom}#ATOM      {atom}#g' pre-mutated.pdb".format(atom=atom_reference+i))
    
    # Save the pre-mutated file
    os.system("grep -v HETATM pre-mutated.pdb > temp; mv temp pre-mutated.pdb")
    
##################################################################################

def mutateRosetta(pdb,pep_chain,new_aa,pep_position,rosetta_path,folder_path):
    """
    Prediction of the mutated side chain using Rosetta functionalities
    
    Arguments:
    pdb -- Code of the structure used to model the peptide
    pep_chain -- Chain containing the sequence where the replacement will be done
    new_aa -- Amino acid that will be mutated
    pep_position -- Position in the chain where the mutation will be done
    rosetta_path -- Path of the Rosetta installation in the local computer
    folder_path -- Path of the folder where the models will be saved
    
    Output:
    A PDB structure with the new predicted side chain atoms
    """
    
    # Dictionary of amino acids
    aminoacids={"ALA":"A","ASP":"D","GLU":"E","PHE":"F","HIS":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","GLY":"G",
                "ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y","CYS":"C"}
    
    # Create Rosetta configuration file
    rosetta_config_file=open("resfile.config","w")
    rosetta_config_file.write("NATRO\n")
    rosetta_config_file.write("start\n")
    rosetta_config_file.write("\t{} {} PIKAA {} EX 1 EX 2 EX 3 EX 4 EX_CUTOFF 0\n".format(pep_position,pep_chain,aminoacids[new_aa]))
    rosetta_config_file.close()
    
    # Run the fixbb tool to mutate the amino acid
    os.system("{}/main/source/bin/fixbb.default.linuxgccrelease -in:file:s pre-mutated.pdb -resfile resfile.config".format(rosetta_path))
    parser = PDBParser()
    structure = parser.get_structure('REF',"pre-mutated_0001.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save("post-pre-mutated.pdb")
    
    # Relax the generated structure
    os.system("{rosetta}/main/source/bin/relax.default.linuxgccrelease -database {rosetta}/main/database -in:file:s post-pre-mutated.pdb -relax:thorough -relax:bb_move false".format(rosetta=rosetta_path))
    parser = PDBParser()
    structure = parser.get_structure('REF',"post-pre-mutated_0001.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save("post-mutated.pdb")
    
    # Delete temporal files and save the model with NNAA replace by a natural amino acid
    os.system("rm pre-mutated.pdb pre-mutated_0001.pdb post-pre-mutated.pdb post-pre-mutated_0001.pdb score.sc resfile.config")
    os.system("mv post-mutated.pdb models/{}/{}_{}_modelled.pdb".format(folder_path,pdb,new_aa))

################################################################################

def relax_ready(pdb,rosetta_path,folder_path):
    """
    Only relax a PDB structure containing a peptide fragment with natural amino acids
    
    Arguments:
    pdb -- Code of the structure used to model the peptide
    rosetta_path -- Path of the Rosetta installation in the local computer
    folder_path -- Path of the folder where the models will be saved
    
    Output:
    Relaxed structure in the folder provided by the user
    """
    
    # Copy the structure of reference
    os.system("grep -v HETATM auxiliar/prepared_structures/{}_fixed.pdb > post-pre-mutated.pdb".format(pdb))
    
    # Run the relax protocol
    os.system("{rosetta}/main/source/bin/relax.default.linuxgccrelease -database {rosetta}/main/database -in:file:s post-pre-mutated.pdb -relax:thorough -relax:bb_move false".format(rosetta=rosetta_path))
    parser = PDBParser()
    structure = parser.get_structure('REF',"post-pre-mutated_0001.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save("post-mutated.pdb")
    
    # Delete temporal files and save the relaxed structure
    os.system("rm post-pre-mutated.pdb post-pre-mutated_0001.pdb score.sc resfile.config")
    os.system("mv post-mutated.pdb models/{}/{}_relaxed.pdb".format(folder_path,pdb))

################################################################################

def similarity(smiles,aaComp):
    """
    Calculate the similarity between a NNAA and a natural amino acid
    
    Arguments:
    smiles -- SMILES representation of the NNAA
    aaComp -- Return the similarity with the amino acid of interest
    
    Return
    similarities[aaComp] -- Similarity values calculated with the RDKit functionalities
    """
    
    # Dictionary of amino acid SMILES
    aminoSmiles = {'GLY':'NCC(=O)O','ALA':'N[C@@]([H])(C)C(=O)O','ARG':'N[C@@]([H])(CCCNC(=N)N)C(=O)O','ASN':'N[C@@]([H])(CC(=O)N)C(=O)O','ASP':'N[C@@]([H])(CC(=O)O)C(=O)O',
                   'CYS':'N[C@@]([H])(CS)C(=O)O','GLU':'N[C@@]([H])(CCC(=O)O)C(=O)O','GLN':'N[C@@]([H])(CCC(=O)N)C(=O)O','HIS':'N[C@@]([H])(CC1=CN=C-N1)C(=O)O',
                   'ILE':'N[C@@]([H])(C(CC)C)C(=O)O','LEU':'N[C@@]([H])(CC(C)C)C(=O)O','LYS':'N[C@@]([H])(CCCCN)C(=O)O','MET':'N[C@@]([H])(CCSC)C(=O)O',
                   'PHE':'N[C@@]([H])(Cc1ccccc1)C(=O)O','PRO':'N1[C@@]([H])(CCC1)C(=O)O','SER':'N[C@@]([H])(CO)C(=O)O','THR':'N[C@@]([H])(C(O)C)C(=O)O',
                   'TRP':'N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)O','TYR':'N[C@@]([H])(Cc1ccc(O)cc1)C(=O)O','VAL':'N[C@@]([H])(C(C)C)C(=O)O'}
    
    # Dictionary storing al the similarities
    similarities={}
    for aa in aminoSmiles:
        mol1=Chem.MolFromSmiles(smiles)
        mol2=Chem.MolFromSmiles(aminoSmiles[aa])
        
        # Generation of fingerprints
        fp1=AllChem.GetMorganFingerprintAsBitVect(mol1,2,2048)
        fp2=AllChem.GetMorganFingerprintAsBitVect(mol2,2,2048)
        similarity=DataStructs.TanimotoSimilarity(fp1,fp2)
        
        # Store the calculated similarities
        similarities[aa]=similarity
    
    # Return the similarity value of interest
    return "{}".format(similarities[aaComp])

####################################################################################

def modelling(list_struct,user,password,database):
    """
    Function to model the complex based on a template containing natural amino acids
    NOTE: To run this function successfully is required a local installation of the MEROPS MySQL instance
    
    Arguments:
    list_struct -- Dictionary with the information about the structures that will be modelled
    user -- MySQL local user information
    password -- MySQL local password information
    database -- MySQL local database information
    
    Return:
    list_to_model_final -- List with the structures that can be modelled based on the assigned criteria
    """
    
    # Dictionaries of amino acids and the conversion between name representations
    aminoacids=["ALA","ASP","GLU","PHE","HIS","ILE","LYS","LEU","MET","GLY","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR","CYS"]
    aminoacid_change={"ALA":"A","ASP":"D","GLU":"E","PHE":"F","HIS":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","GLY":"G",
                "ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y","CYS":"C"}
    
    # General Connection to the MySQL database
    connection={'host':'localhost','user':user,'pwd':password,'database':database} 
    con = mdb.connect(connection['host'],connection['user'],connection['pwd'],connection['database']);
    substrates_chosen={}
    
    # Creation of the list that will store the information for modelling
    list_to_model_final=[]
    
    # Iterate over the proteins from the list
    for protein in list_struct:
        # Store some initial information
        mer=list_struct[protein][2]
        ligand_chain=list_struct[protein][4]
        if mer not in substrates_chosen: substrates_chosen[mer]=[]
        
        # Run the MySQL query
        query="select * from Substrate_search where code='{}'".format(mer)  
        cur = con.cursor()
        cur.execute(query)
        results=cur.fetchall()
        
        counter=0
        
        # Store additional information from the structure
        amino_subs=list_struct[protein][0]
        position_subs=list_struct[protein][1]
        number_subs=list_struct[protein][3]
        
        # Run the analysis if results were obtained from the query
        if results:
            for r in results:
                # Store the amino acids of the substrate at each position
                incomplete=0
                positions={}
                positions['S4']=r[4].upper()
                positions['S3']=r[5].upper()
                positions['S2']=r[6].upper()
                positions['S1']=r[7].upper()
                positions["S1'"]=r[8].upper()
                positions["S2'"]=r[9].upper()
                positions["S3'"]=r[10].upper()
                positions["S4'"]=r[11].upper()
                
                # Check that the substrate has store the full sequence
                for key in positions:
                    if positions[key]=="-": incomplete=1
                
                # Store additional information from the query
                uniprot=r[13]
                organism=r[15]

                # Check how many amino acids are the same between the template and the sequence to model
                matches=0
                for i,realPos in enumerate(position_subs):
                    if list_struct[protein][0][i]==positions[realPos]:
                        matches+=1
                
                # Criterum to move on with the modelling - Check for human substrates with just one amino acid difference and not repeated full substrate
                if matches>=1 and len(amino_subs)<=matches+1 and organism=="Homo sapiens" and uniprot not in substrates_chosen[mer] and incomplete==0:
                    substrates_chosen[mer].append(uniprot)
                    
                    # By default put as None the replacement of amino acids
                    old_aa=None
                    new_aa=None
                    aa_position=None
                    
                    # Check if there is a mutation to store the information about the change
                    for amino in amino_subs:
                        aa1=amino
                        aa2=positions[position_subs[amino_subs.index(amino)]]
                        if aa1!=aa2:
                            old_aa=aa1
                            new_aa=aa2
                            aa_position=number_subs[amino_subs.index(amino)]
                    
                    # Annotate the peptide sequence that will be used as template and the one that will be modelled                
                    pepTemplate=''
                    pepModel=''
                    for pos_pep in ['S4','S3','S2','S1',"S1'","S2'","S3'","S4'"]:
                        if pos_pep in position_subs:
                            pepTemplate+=aminoacid_change[positions[pos_pep]]
                            pepModel+=aminoacid_change[positions[pos_pep]]
                        else:
                            pepTemplate+='-'
                            pepModel+=aminoacid_change[positions[pos_pep]]
                    
                    # Store the information in the final list and print it
                    array_info=[protein,mer,pepTemplate,pepModel,old_aa,new_aa,aa_position,ligand_chain,uniprot]
                    if array_info not in list_to_model_final:
                        list_to_model_final.append(array_info)
                        print(array_info)
    
    # Return the list with the structures to model
    return list_to_model_final

####################################################################################

def modelling_complete(list_struct,sim_threshold,user,password,database):
    """
    Function to model the complex based on a template containing a non-natural amino acid (NNAA)
    NOTE: To run this function successfully is required a local installation of the MEROPS MySQL instance
    
    Arguments:
    list_struct -- Dictionary with the information about the structures that will be modelled
    sim_threshold -- Similarity threshold to decide wiich natural amino acids will replace the NNAA
    user -- MySQL local user information
    password -- MySQL local password information
    database -- MySQL local database information
    
    Return:
    list_to_model_final -- List with the structures that can be modelled based on the assigned criteria
    """
    
    # Dictionaries of amino acids and the conversion between name representations
    aminoacids=["ALA","ASP","GLU","PHE","HIS","ILE","LYS","LEU","MET","GLY","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR","CYS"]
    aminoacid_change={"ALA":"A","ASP":"D","GLU":"E","PHE":"F","HIS":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","GLY":"G",
                "ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y","CYS":"C"}
    
    # General Connection to the MySQL database
    connection={'host':'localhost','user':user,'pwd':password,'database':database} 
    con = mdb.connect(connection['host'],connection['user'],connection['pwd'],connection['database']);
    substrates_chosen={}
    
    # Creation of the list that will store the information for modelling
    list_to_model_final=[]
    
    # Iterate over the proteins from the list
    for protein in list_struct:
        # Store some initial information
        mer=list_struct[protein][2]
        ligand_chain=list_struct[protein][4]
        if mer not in substrates_chosen: substrates_chosen[mer]=[]
        
        # Run the MySQL query
        query="select * from Substrate_search where code='{}'".format(mer)
        cur = con.cursor()
        cur.execute(query)
        results=cur.fetchall()
        
        counter=0
        
        # Store additional information from the structure
        amino_subs=list_struct[protein][0]
        position_subs=list_struct[protein][1]
        number_subs=list_struct[protein][3]
        
        # Pre-calculate the similarity between the NNAA and all the natural amino acids
        # NOTE: The smiles from the NNAA is obtained from the components.cif file from the PDB
        sim_associations=[]
        for amino in amino_subs:
            if amino not in aminoacids:
                old_aa=amino
                database=open("auxiliar/components.cif","r").read()
                index = re.findall('{} SMILES_CANONICAL CACTVS.*$'.format(old_aa), database,re.MULTILINE)
                smiles=index[0].strip().split()[-1].strip("\"")
                for aa in aminoacids:
                    # Run the similarity function
                    sim=similarity(smiles,aa)
                    sim_associations.append((old_aa,aa,sim))
        
        # Run the analysis if results were obtained from the query
        if results:
            for r in results:
                # Store the amino acids of the substrate at each position
                incomplete=0
                positions={}
                positions['S4']=r[4].upper()
                positions['S3']=r[5].upper()
                positions['S2']=r[6].upper()
                positions['S1']=r[7].upper()
                positions["S1'"]=r[8].upper()
                positions["S2'"]=r[9].upper()
                positions["S3'"]=r[10].upper()
                positions["S4'"]=r[11].upper()
                
                # Check that the substrate has store the full sequence
                for key in positions:
                    if positions[key]=="-": incomplete=1
                
                # Store additional information from the query
                uniprot=r[13]
                organism=r[15]
                
                # Check how many amino acids are the same between the template and the sequence to model
                matches=0
                for i,realPos in enumerate(position_subs):
                    if list_struct[protein][0][i]==positions[realPos]:
                        matches+=1
                
                # Criterum to move on with the modelling - Check for human substrates with just one amino acid difference and not repeated full substrate
                if matches>=1 and len(amino_subs)<=matches+1 and organism=="Homo sapiens" and uniprot not in substrates_chosen[mer] and incomplete==0:
                    substrates_chosen[mer].append(uniprot)
                    
                    # Obtain the similarity between the amino acids based on the pre-calculated values
                    for amino in amino_subs:
                        if amino not in aminoacids:
                            old_aa=amino
                            new_aa=positions[position_subs[amino_subs.index(amino)]]
                            aa_position=number_subs[amino_subs.index(amino)]
                            for sim_pairs in sim_associations:
                                if old_aa==sim_pairs[0] and new_aa==sim_pairs[1]:
                                    sim=float(sim_pairs[2])
                    
                    # Annotate the peptide sequence that will be used as template and the one that will be modelled
                    pepTemplate=''
                    pepModel=''
                    for pos_pep in ['S4','S3','S2','S1',"S1'","S2'","S3'","S4'"]:
                        if pos_pep in position_subs:
                            aa_ref=amino_subs[position_subs.index(pos_pep)]
                            if aa_ref not in aminoacids: aa_ref=new_aa
                            pepTemplate+=aminoacid_change[aa_ref]
                            pepModel+=aminoacid_change[aa_ref]
                        else:
                            pepTemplate+='-'
                            pepModel+=aminoacid_change[positions[pos_pep]]
                    
                    # Store the information in the final list and print it
                    array_info=[protein,mer,pepTemplate,pepModel,sim,old_aa,new_aa,aa_position,ligand_chain,uniprot]
                    if array_info not in list_to_model_final:
                        # Check the similarity threshold to store the data
                        if sim >= sim_threshold:
                            list_to_model_final.append(array_info)
                            print(array_info)
    
    # Return the list with the structures to model
    return list_to_model_final

####################################################################################

def run_model(model,rosetta_path,folder_path):
    """
    Function to run the complete model protocol for proteases bound to peptides with natural amino acids
    
    Arguments:
    model -- Element from a list with the required information to run the modelling
    rosetta_path -- Path of the Rosetta installation in the local computer
    folder_path -- Path of the folder where the models will be saved
    
    Output:
    Run all the modelling functions, save the models and create a file with the report
    """
    
    # Information required to generate the csv report of the modelling
    columns=["pdb","pepTemplate","pepModel","chain","merops","old_aa","new_aa","pos_aa","uniprot"]
    to_model_data=[]
    
    # Information extracted from the list element
    pdb=model[0]
    mer=model[1]
    pepT=model[2]
    pepM=model[3]
    old_aa=model[4]
    new_aa=model[5]
    pos_aa=model[6]
    chain_aa=model[7]
    uniprot=model[8]
    
    # Save information for the csv report
    to_model_data.append({"pdb":pdb,"pepTemplate":pepT,"pepModel":pepM,"chain":chain_aa,"merops":mer,"old_aa":old_aa,"new_aa":new_aa,"pos_aa":pos_aa,"uniprot":uniprot})
    
    # Run the relaxation of the structure
    relax_ready(pdb,rosetta_path,folder_path)
    
    # If an amino acid require to be mutated, run the mutation protocol
    if new_aa:   
        parser = PDBParser()
        reference = parser.get_structure('REF',"models/{}/{}_relaxed.pdb".format(folder_path,pdb))
        print(reference,chain_aa,old_aa,new_aa,pos_aa)
        replace_amino_acid_natural(reference,chain_aa,old_aa,new_aa,pos_aa)
        mutateRosetta(pdb,chain_aa,new_aa,pos_aa,rosetta_path,folder_path)
        
    # If missing amino acids should be modelled, run the protocol
    if '-' in pepT:
        modelling_fragment(pdb,chain_aa,pepT,pepM,rosetta_path,folder_path,new_aa)
    
    # Open a report file with the information used for modelling
    with open("report_model_ready.csv", 'w') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=columns)
        writer.writeheader()
        for data in to_model_data: writer.writerow(data)

####################################################################################

def run_model_complete(model,rosetta_path,folder_path):
    """
    Function to run the complete model protocol for proteases bound to peptides with non-natural amino acids (NNAA)
    
    Arguments:
    model -- Element from a list with the required information to run the modelling
    rosetta_path -- Path of the Rosetta installation in the local computer
    folder_path -- Path of the folder where the models will be saved
    
    Output:
    Run all the modelling functions, save the models and create a file with the report
    """
    
    # Information required to generate the csv report of the modelling
    columns=["pdb","pepTemplate","pepModel","chain","merops","old_aa","new_aa","pos_aa","uniprot","sim"]
    to_model_data=[]
    
    # Information extracted from the list element
    pdb=model[0]
    mer=model[1]
    pepT=model[2]
    pepM=model[3]
    sim=model[4]
    old_aa=model[5]
    new_aa=model[6]
    pos_aa=model[7]
    chain_aa=model[8]
    uniprot=model[9]
    
    # Save information for the csv report
    to_model_data.append({"pdb":pdb,"pepTemplate":pepT,"pepModel":pepM,"chain":chain_aa,"merops":mer,"old_aa":old_aa,"new_aa":new_aa,"pos_aa":pos_aa,"uniprot":uniprot,"sim":sim})
    
    # Run the mutation protocol
    parser = PDBParser()
    reference = parser.get_structure('REF',"auxiliar/prepared_structures/{}_fixed.pdb".format(pdb))
    print(reference,chain_aa,old_aa,new_aa,pos_aa)
    replace_amino_acid(reference,chain_aa,old_aa,new_aa,pos_aa)
    mutateRosetta(pdb,chain_aa,new_aa,pos_aa,rosetta_path,folder_path)
    
    # If missing amino acids should be modelled, run the protocol
    if '-' in pepT:
        modelling_fragment(pdb,chain_aa,pepT,pepM,rosetta_path,folder_path,new_aa)
        
    # Open a report file with the information used for modelling
    with open("report_model_complete.csv", 'w') as csvfile:
         writer = csv.DictWriter(csvfile, fieldnames=columns)
         writer.writeheader()
         for data in to_model_data: writer.writerow(data)

########################################################################################
########################################################################################
########################################################################################
# Main execution
########################################################################################
########################################################################################
########################################################################################
if __name__ == '__main__':
    
    # Script arguments
    parser = argparse.ArgumentParser(description='Package to model peptide substrates bound to annotated protease structures')
    parser.add_argument('-s', dest='structure', action='store',default='none',
                        help='Add a protease structure PDB id to see if it is possible to model substrates')
    parser.add_argument('-f', dest='family', action='store',default='serine',
                        help='Select the family of proteases used as reference to model the substrates from 1) serine, 2) cysteine')
    parser.add_argument('-n', dest='max_substrates', action='store',default=1,
                        help='Define the maximum number of substrates that can be modelled per structure available')
    parser.add_argument('-m', dest='mode', action='store', default="ready",
                        help='Choose a mode to run the script from two options: 1) ready, 2) complete. Ready are the peptide templates with natural amino acids, and Complete are the peptides containing at least one amino acid')
    parser.add_argument('-t', dest='sim_threshold', action='store',default=0.4,
                        help='Similarity threshold to decide which amino acid will replace the NNAA present in the template')
    parser.add_argument('-r', dest='rosetta', action='store', default="rosetta_src_2016.32.58837_bundle",
                        help='Version of Rosetta that will be implemented')
    parser.add_argument('-u', dest='user', action='store', required=True,
                        help='User to access the local installation of the MEROPS database')
    parser.add_argument('-p', dest='password', action='store', required=True,
                        help='Password to access the local installation of the MEROPS database')
    parser.add_argument('-d', dest='database', action='store', required=True,
                        help='Name of the local installation of the MEROPS database')
    args = parser.parse_args()
    
    ####################################################################################
    # Assignment of parameters
    
    # Rosetta version installed in the computer
    rosetta_version=args.rosetta
    bash = "locate -b {} | head -n1".format(rosetta_version)
    rosetta_path = subprocess.check_output(['bash','-c', bash]).strip()
    
    # MySQL data. NOTE: it depends on the local installation of the MEROPS database. See EXTRA FILES section at the beginning of the script
    user=args.user
    password=args.password
    database=args.database
    
    # Obtain the dictionaries from the csv file of interest
    family=args.family
    if family=="serine":
        ready,model=map_peptides('auxiliar/serine_proteases_pockets.csv')
    if family=="cysteine":
        ready,model=map_peptides('auxiliar/cysteine_proteases_pockets.csv')
    
    #print ready
    mode_modelling=args.mode
    struct=args.structure.lower()
    sim_threshold=args.sim_threshold
    max_subs=int(args.max_substrates)
    
    # Create folders in case they do not exist
    os.system("mkdir models")
    os.system("mkdir models/model_ready")
    os.system("mkdir models/model_complete")
    
    # Check if the mode of modelling is for structures having all natural amino acids
    if mode_modelling=="ready":
        # Check the structure 
        if struct!="none":
            if struct in ready:
                print("The structure {} is present in the annotated set".format(struct))
                # Run the modelling process for the complexes with natural amino acids
                list_to_model=modelling(ready,user,password,database)
                print("######################")
                # Check if the model is in the list of structures that can be modelled
                flag_stop=1
                for m in list_to_model:
                    protein=m[0]
                    if struct==protein: flag_stop=0
                
                # Check the flag and report the message
                if flag_stop==1:
                    print("The structure {} does not have substrate available to do the modelling".format(struct))
                    print("Exiting ...")
                    exit()
                else:
                    print("The structure {} has substrates available to do the modelling".format(struct))
            # The structure is not in the original dataset
            else:
                print("The structure {} is not available in the annotated set and the modelling will stop".format(struct))
                print("Exiting ...")
                exit()
        else:
            # Run the modelling process for the complexes with natural amino acids
            list_to_model=modelling(ready,user,password,database)
            print("######################")
            print(list_to_model)
                
        # Move to model the substrate if there is information available
        # Case where the pdb structure is directly provided in the script
        if struct!="none":
            counter=1
            for m in list_to_model:
                protein=m[0]
                if struct==protein:
                    # Check the maximum number of structures required per protein
                    if counter <= max_subs:
                        run_model(m,rosetta_path,"model_ready")
                        counter+=1
                    else:
                        print("The number of models required for {} has been completed".format(struct))
                        print("Exiting ...")
                        exit()
        else:
            # General case where an n number of models is generated per structure
            counter_proteins={}
            prot_ref=""
            for m in list_to_model:
                protein=m[0]
                if protein!=prot_ref: counter=1
                prot_ref=protein
                counter_proteins[protein]=counter
                # Check the maximum number of structures required per protein
                if counter_proteins[protein] <= max_subs:
                    run_model(m,rosetta_path,"model_ready")
                    counter+=1
                else:
                    print("The number of models required for {} has been completed".format(struct))
                    print("Exiting ...")
                    exit()

    # Check if the mode of modelling is for structures having one NNAA to replace by    
    if mode_modelling=="complete":
        #print("Chao")
        # Check the structure 
        if struct!="none":
            if struct in model:
                print("The structure {} is present in the annotated set".format(struct))
                # Run the modelling process for the complexes with NNAAs
                list_to_model=modelling_complete(model,sim_threshold,user,password,database)
                print("######################")
                # Check if the model is in the list of structures that can be modelled
                flag_stop=1
                for m in list_to_model:
                    protein=m[0]
                    if struct==protein: flag_stop=0
                
                # Check the flag and report the message
                if flag_stop==1:
                    print("The structure {} does not have substrate available to do the modelling".format(struct))
                    print("Exiting ...")
                    exit()
                else:
                    print("The structure {} has substrates available to do the modelling".format(struct))
            # The structure is not in the original dataset
            else:
                print("The structure {} is not available in the annotated set and the modelling will stop".format(struct))
                print("Exiting ...")
                exit()
        else:
            # Run the modelling process for the complexes with NNAAs
            list_to_model=modelling_complete(model,sim_threshold,user,password,database) # PENDING READY
            print("######################")
            print(list_to_model)
        
        # Move to model the substrate if there is information available
        # Case where the pdb structure is directly provided in the script
        if struct!="none":
            counter=1
            for m in list_to_model:
                protein=m[0]
                if struct==protein:
                    # Check the maximum number of structures required per protein
                    if counter <= max_subs:
                        run_model_complete(m,rosetta_path,"model_complete")
                        counter+=1
                    else:
                        print("The number of models required for {} has been completed".format(struct))
                        print("Exiting ...")
                        exit()
        else:
            # General case where an n-number of models is generated per structure
            counter_proteins={}
            prot_ref=""
            for m in list_to_model:
                protein=m[0]
                if protein!=prot_ref: counter=1
                prot_ref=protein
                counter_proteins[protein]=counter
                # Check the maximum number of structures required per protein
                if counter_proteins[protein] <= max_subs:
                    run_model_complete(m,rosetta_path,"model_complete")
                    counter+=1
                else:
                    print("The number of models required for {} has been completed".format(struct))
                    print("Exiting ...")
                    exit()