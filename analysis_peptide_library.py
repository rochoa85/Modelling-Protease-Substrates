#!/usr/bin/python

"""
Protocol to calculate the structural observables for a library of peptides bound to a modelled protease structure
NOTE: The protocol requires of auxiliary programs and Unix system commands - Tested on Ubuntu 16.04

From publication "Modelling peptide substrates bound to proteases: insights to study binding promiscuity"
Journal: PLoS Computational Biology
Authors: Rodrigo Ochoa, Mikhail Magnitov, Roman Laskowski, Pilar Cossio, Janet Thornton
Year: 2020
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

import os
import subprocess
import sys
import numpy as np
from statistics import mean

########################################################################################
# Main protocol
########################################################################################

# Dictionary with amino acids    
aminoacid_change={"ALA":"A","ASP":"D","GLU":"E","PHE":"F","HIS":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","GLY":"G",
                  "ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y","CYS":"C"}
aminoacids={"A":"ALA","D":"ASP","E":"GLU","F":"PHE","H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","G":"GLY",
            "N":"ASN","P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR","C":"CYS"}

# Dictionary to normalize the calculated energy per amino acid
norm_energy={"A":88.6,"D":111.1,"E":138.4,"F":189.9,"H":153.2,"I":166.7,"K":168.6,"L":166.7,"M":162.9,"G":60.1,
            "N":114.1,"P":112.7,"Q":143.8,"R":173.4,"S":89.0,"T":116.1,"V":140.0,"W":227.8,"Y":193.6,"C":108.5}

# Dictionaries storing the average values per peptide position
P4_avg_asa={}; P4_avg_energy={}
P3_avg_asa={}; P3_avg_energy={}
P2_avg_asa={}; P2_avg_energy={}
P1_avg_asa={}; P1_avg_energy={}
P1p_avg_asa={}; P1p_avg_energy={}
P2p_avg_asa={}; P2p_avg_energy={}
P3p_avg_asa={}; P3p_avg_energy={}
P4p_avg_asa={}; P4p_avg_energy={}

# Read the list of peptide from the library
peptides=[x.strip() for x in open("list_peptide_to_model.txt")]

# Iterate over the peptides
for count,pep in enumerate(peptides):

    print "Processing peptide %s ..." %pep
    # Iterate over the peptide amino acids
    for position in range(1,9):
        aa_single=pep[position-1]
        aa=aminoacids[pep[position-1]]
        chain='B'
        
        ######################################3
        # CREATE FOLDER, copy files, move there and run the script
        ######################################3
        
        # Read the ASA values
        asa_prev=[float(x.strip()) for x in open("%s/asa/asa_BR_%s%d_%s.txt" %(pep,aa,position,pep))]
        averageAsa=mean(asa_prev)
        
        # Read the energy values
        energy_prev=[float(x.strip().split()[1]) for x in open("%s/energy/energy_BR_%s%d_%s.txt" %(pep,aa,position,pep))]
        averageEnergy=mean(energy_prev)
    
        # Store the values per position in the peptide
        if position==1:
            if aa_single not in P4_avg_asa:
                P4_avg_asa[aa_single]=[averageAsa]
                P4_avg_energy[aa_single]=[averageEnergy]
            else:
                P4_avg_asa[aa_single].append(averageAsa)
                P4_avg_energy[aa_single].append(averageEnergy)
        if position==2:
            if aa_single not in P3_avg_asa:
                P3_avg_asa[aa_single]=[averageAsa]
                P3_avg_energy[aa_single]=[averageEnergy]
            else:
                P3_avg_asa[aa_single].append(averageAsa)
                P3_avg_energy[aa_single].append(averageEnergy)
        if position==3:
            if aa_single not in P2_avg_asa:
                P2_avg_asa[aa_single]=[averageAsa]
                P2_avg_energy[aa_single]=[averageEnergy]
            else:
                P2_avg_asa[aa_single].append(averageAsa)
                P2_avg_energy[aa_single].append(averageEnergy)
        if position==4:
            if aa_single not in P1_avg_asa:
                P1_avg_asa[aa_single]=[averageAsa]
                P1_avg_energy[aa_single]=[averageEnergy]
            else:
                P1_avg_asa[aa_single].append(averageAsa)
                P1_avg_energy[aa_single].append(averageEnergy)
        if position==5:
            if aa_single not in P1p_avg_asa:
                P1p_avg_asa[aa_single]=[averageAsa]
                P1p_avg_energy[aa_single]=[averageEnergy]
            else:
                P1p_avg_asa[aa_single].append(averageAsa)
                P1p_avg_energy[aa_single].append(averageEnergy)
        if position==6:
            if aa_single not in P2p_avg_asa:
                P2p_avg_asa[aa_single]=[averageAsa]
                P2p_avg_energy[aa_single]=[averageEnergy]
            else:
                P2p_avg_asa[aa_single].append(averageAsa)
                P2p_avg_energy[aa_single].append(averageEnergy)
        if position==7:
            if aa_single not in P3p_avg_asa:
                P3p_avg_asa[aa_single]=[averageAsa]
                P3p_avg_energy[aa_single]=[averageEnergy]
            else:
                P3p_avg_asa[aa_single].append(averageAsa)
                P3p_avg_energy[aa_single].append(averageEnergy)
        if position==8:
            if aa_single not in P4p_avg_asa:
                P4p_avg_asa[aa_single]=[averageAsa]
                P4p_avg_energy[aa_single]=[averageEnergy]
            else:
                P4p_avg_asa[aa_single].append(averageAsa)
                P4p_avg_energy[aa_single].append(averageEnergy)

# Create folders in the peptide workspace
os.system("mkdir matrices")

# Calculate the mean values of ASA among all the peptides per position
means_asa={'P4':{},'P3':{},'P2':{},'P1':{},'P1p':{},'P2p':{},'P3p':{},'P4p':{}}
for aa in aminoacids:
    if aa in P4_avg_asa: means_asa['P4'][aa]=mean(P4_avg_asa[aa])
    else: means_asa['P4'][aa]= 0.0
    
    if aa in P3_avg_asa: means_asa['P3'][aa]=mean(P3_avg_asa[aa])
    else: means_asa['P3'][aa]= 0.0
    
    if aa in P2_avg_asa: means_asa['P2'][aa]=mean(P2_avg_asa[aa])
    else: means_asa['P2'][aa]= 0.0
    
    if aa in P1_avg_asa: means_asa['P1'][aa]=mean(P1_avg_asa[aa])
    else: means_asa['P1'][aa]= 0.0
    
    if aa in P1p_avg_asa: means_asa['P1p'][aa]=mean(P1p_avg_asa[aa])
    else: means_asa['P1p'][aa]= 0.0
    
    if aa in P2p_avg_asa: means_asa['P2p'][aa]=mean(P2p_avg_asa[aa])
    else: means_asa['P2p'][aa]= 0.0
    
    if aa in P3p_avg_asa: means_asa['P3p'][aa]=mean(P3p_avg_asa[aa])
    else: means_asa['P3p'][aa]= 0.0
    
    if aa in P4p_avg_asa: means_asa['P4p'][aa]=mean(P4p_avg_asa[aa])
    else: means_asa['P4p'][aa]= 0.0

# Generating the files with the complete average values in the library
print "Creating means asa ..."
means_asa_file=open("matrices/means_full_asa.txt","wb")
means_asa_file.write("AA\tP4\tP3\tP2\tP1\tP1p\tP2p\tP3p\tP4p\n")
for aa in P4_avg_asa:
    means_asa_file.write("%s\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n" %(aa,means_asa['P4'][aa],
                             means_asa['P3'][aa],means_asa['P2'][aa],means_asa['P1'][aa],means_asa['P1p'][aa],
                             means_asa['P2p'][aa],means_asa['P3p'][aa],means_asa['P4p'][aa]))
means_asa_file.close()

# Calculate the mean values of the energy among all the peptides per position
means_energy={'P4':{},'P3':{},'P2':{},'P1':{},'P1p':{},'P2p':{},'P3p':{},'P4p':{}}
for aa in aminoacids:
    if aa in P4_avg_energy: means_energy['P4'][aa]=mean(P4_avg_energy[aa])
    else: means_energy['P4'][aa]= 0.0
    
    if aa in P3_avg_energy: means_energy['P3'][aa]=mean(P3_avg_energy[aa])
    else: means_energy['P3'][aa]= 0.0
    
    if aa in P2_avg_energy: means_energy['P2'][aa]=mean(P2_avg_energy[aa])
    else: means_energy['P2'][aa]= 0.0
    
    if aa in P1_avg_energy: means_energy['P1'][aa]=mean(P1_avg_energy[aa])
    else: means_energy['P1'][aa]= 0.0
    
    if aa in P1p_avg_energy: means_energy['P1p'][aa]=mean(P1p_avg_energy[aa])
    else: means_energy['P1p'][aa]= 0.0
    
    if aa in P2p_avg_energy: means_energy['P2p'][aa]=mean(P2p_avg_energy[aa])
    else: means_energy['P2p'][aa]= 0.0
    
    if aa in P3p_avg_energy: means_energy['P3p'][aa]=mean(P3p_avg_energy[aa])
    else: means_energy['P3p'][aa]= 0.0
    
    if aa in P4p_avg_energy: means_energy['P4p'][aa]=mean(P4p_avg_energy[aa])
    else: means_energy['P4p'][aa]= 0.0

# Generating the files with the complete average values in the library
print "Creating means energy ..."
means_energy_file=open("matrices/means_full_energy.txt","w")
means_energy_file.write("AA\tP4\tP3\tP2\tP1\tP1p\tP2p\tP3p\tP4p\n")
for aa in P4_avg_energy:
    means_energy_file.write("%s\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n" %(aa,means_energy['P4'][aa],
                             means_energy['P3'][aa],means_energy['P2'][aa],means_energy['P1'][aa],means_energy['P1p'][aa],
                             means_energy['P2p'][aa],means_energy['P3p'][aa],means_energy['P4p'][aa]))
means_energy_file.close()