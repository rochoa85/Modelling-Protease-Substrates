#!/usr/bin/python

"""
Package to run dynamic analysis of peptide substrates bound to annotated protease structures
NOTE: The protocol requires of auxiliary programs and Unix system commands - Tested on Ubuntu 16.04

From publication "Modelling peptide substrates bound to proteases: insights to study binding promiscuity"
Journal: BMC Bioinformatics
Authors: Rodrigo Ochoa, Mikhail Magnitov, Roman Laskowski, Pilar Cossio, Janet Thornton
Year: 2020

Third-party tools required:

BioPython: https://biopython.org/wiki/Download
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

import subprocess
import os
import sys
import argparse
import numpy as np
from random import shuffle
from statistics import mean

# Third-party modules: Biopython
from Bio.PDB import *

########################################################################################
# Functions
########################################################################################

def replace_amino_acid(pep_pdb,pep_chain,old_aa,new_aa,pep_position):
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
    
    aminoacids={"A":"ALA","D":"ASP","E":"GLU","F":"PHE","H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","G":"GLY",
                "N":"ASN","P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR","C":"CYS"}

    # Read the PDB file
    residues=pep_pdb.get_residues()
    chain=pep_pdb[0][pep_chain]
    
    # Report the mutation made
    message="The residue {} in chain {} and position {} will be changed by {}".format(old_aa,pep_chain,pep_position,new_aa)
    print(message)
    
    # Rename the residue
    chain[pep_position].resname=aminoacids[new_aa]
    
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

def generateMutations(pep1,pep2):
    """
    Based on two sequences, generate the pairs required to perform single point mutations on the peptide
    
    Arguments:
    pep1 -- Sequence of reference
    pep2 -- Sequence that will be modelled
    
    Return:
    mutation_list -- A list with tuples of the mutation that will be performed
    """
    
    # List that will contain the mutations to be performed
    mutation_list=[]
    
    # Check if the positions contains X characters, and store the pairs if the old and new AA are different
    if pep1[0]=="X":
        counter=0
        for i,aa in enumerate(pep1):
            if aa=="X": pass
            else:
                counter+=1
                if aa!=pep2[i]:
                    mutation_list.append((aa,pep2[i],counter))
    else:
        for i,aa in enumerate(pep1):
            if aa!=pep2[i]:
                mutation_list.append((aa,pep2[i],i+1))
                
    # Return the list
    return mutation_list

##################################################################################

def generatePairs(pep1,pep2):
    """
    Based on two sequences, generate the all the possible pairs position by position
    
    Arguments:
    pep1 -- Sequence of reference
    pep2 -- Sequence that will be modelled
    
    Return:
    mutation_list -- A list with tuples all the assignments position by position
    """
    
    # List that will contain all the associations
    mutation_list=[]
    for i,aa in enumerate(pep1):
        mutation_list.append((aa,pep2[i],i+1))
        
    # Return the list
    return mutation_list

##################################################################################

def mutateRosetta(pep_chain,new_aa,pep_position,rosetta_path):
    """
    Prediction of the mutated side chain using Rosetta functionalities
    
    Arguments:
    pep_chain -- Chain containing the sequence where the replacement will be done
    new_aa -- Amino acid that will be mutated
    pep_position -- Position in the chain where the mutation will be done
    rosetta_path -- Path of the Rosetta installation in the local computer
    
    Output:
    A PDB structure called "pos-mutated.pdb" with the new predicted side chain atoms
    """
    
    # Creation of the Rosetta configuration file
    rosetta_config_file=open("resfile.config","w")
    rosetta_config_file.write("NATRO\n")
    rosetta_config_file.write("start\n")
    rosetta_config_file.write("\t{} {} PIKAA {} EX 1 EX 2 EX 3 EX 4 EX_CUTOFF 0\n".format(pep_position,pep_chain,new_aa))
    rosetta_config_file.close()
    
    # Perform the prediction of the side chain
    os.system("{}/main/source/bin/fixbb.default.linuxgccrelease -in:file:s pre-mutated.pdb -resfile resfile.config".format(rosetta_path))
    parser = PDBParser()
    structure = parser.get_structure('REF',"pre-mutated_0001.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save("post-pre-mutated.pdb")
    
    # Relax the system
    os.system("{rosetta}/main/source/bin/relax.default.linuxgccrelease -database {rosetta}/main/database -in:file:s post-pre-mutated.pdb -relax:thorough -relax:bb_move false".format(rosetta=rosetta_path))
    parser = PDBParser()
    structure = parser.get_structure('REF',"post-pre-mutated_0001.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save("post-mutated.pdb")

##################################################################################

def modelling(mutList,model,pep_to_model,pepChain,rosetta_path):
    """
    Modelling of the new peptide by iterative single-point mutations
        
    Arguments:
    mutList -- List containing all the mutations that will be performed
    model -- Route of the model PDB file that will be used as template
    pep_to_model -- Original sequence of the peptide to model
    pepChain -- Chain of the peptide in the structure
    rosetta_path -- Path of the Rosetta installation in the local computer
        
    Output:
    PDB file of the peptide mutated (XXX_mutated_final.pdb)
    """
    
    # Iterate over the list of pairs that will be mutated
    counter_step=0
    for j,mut in enumerate(mutList):
        old_aa=mut[0]
        new_aa=mut[1]
        position=mut[2]
        
        # Run the mutations only if there are not amino acids to insert. In this case this happen always
        if old_aa != "X" and new_aa != "X":
            print(mut)
            if counter_step==0:
                parser = PDBParser()
                reference = parser.get_structure('REF',"{}".format(model))
                counter_step+=1
            else:
                parser = PDBParser()
                reference = parser.get_structure('REF',"post-mutated.pdb")
            
            # Replace the amino acid and run the prediction of the new side chain rotamer
            replace_amino_acid(reference,pepChain,old_aa,new_aa,position)
            mutateRosetta(pepChain,new_aa,position,rosetta_path)
            os.system("rm pre-mutated.pdb pre-mutated_0001.pdb post-pre-mutated.pdb post-pre-mutated_0001.pdb resfile.config score.sc")
    
    # Delete previous files
    os.system("mv post-mutated.pdb {}_pre_mutation.pdb".format(pep_to_model))
    
    # Relaxing of the final complex
    os.system("{rosetta}/main/source/bin/relax.default.linuxgccrelease -database {rosetta}/main/database -in:file:s {pep}_pre_mutation.pdb -relax:thorough -relax:bb_move false".format(rosetta=rosetta_path,pep=pep_to_model))
    parser = PDBParser()
    structure = parser.get_structure('REF',"{}_pre_mutation_0001.pdb".format(pep_to_model))
    io = PDBIO()
    io.set_structure(structure)
    io.save("post-modelled.pdb")
    
    # Delete temporal files and rename the final structure
    os.system("rm score.sc resfile.config {pep}_pre_mutation_0001.pdb {pep}_pre_mutation.pdb".format(pep=pep_to_model))
    os.system("mv post-modelled.pdb {}_mutated_final.pdb".format(pep_to_model))
    
##################################################################################

def run_sampling(rosetta_path,pep_to_model):
    """
    Function to run the Backrub analysis and refinement of the modelled protease-peptide complex
    
    Arguments:
    rosetta_path -- Path of the Rosetta installation in the local computer
    pep_to_model -- Original sequence of the peptide to model
    
    Output:
    Trajectory calculated using Backrub with parameters defined by the user
    """
    ################
    # Pendiente cuadrar el refinement con flexpepdock. Eso esta listo en servidor .9
    ################
    
    # Runnning Backrub
    os.system("{rosetta}/main/source/bin/backrub.default.linuxgccrelease -database {rosetta}/main/database -s {pep}_mutated_final.pdb -ignore_unrecognized_res -ex1 -ex2 -extrachi_cutoff 0 -backrub:ntrials 10000 -mc_kt 1.2 -ignore_zero_occupancy=false -trajectory=true -initial_pack".format(rosetta=rosetta_path,pep=pep_to_model))
    
    # Delete additional files
    os.system("rm {pep}_mutated_final_0001_last.pdb".format(pep=pep_to_model))
    os.system("rm {pep}_mutated_final_0001_low.pdb".format(pep=pep_to_model))
    os.system("rm {pep}_mutated_final_0001.pdb".format(pep=pep_to_model))
    os.system("rm score.sc")

##################################################################################
    
def calculate_observables(pep_to_model,pepChain,rosetta_path):
    """
    Function to calculate the observables from the Backrub trajectory
    
    Arguments:
    pep_to_model -- Original sequence of the peptide to model
    pepChain -- Chain of the peptide in the structure
    rosetta_path -- Path of the Rosetta installation in the local computer
        
    Return:
    asaTotal -- dictionary containing all the asa values calculated from the snapshots
    energyTotal -- dictionary containing all the interface energy scores calculated from the snapshots
    """
    aminoacids={"A":"ALA","D":"ASP","E":"GLU","F":"PHE","H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","G":"GLY",
                "N":"ASN","P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR","C":"CYS"}
    
    # Input trajectory information
    traj="{}_mutated_final_0001_traj.pdb".format(pep_to_model)
    compareChain='A' # First chain of the pdb
    position_reference=1
    
    # Upload the the text from the trajectory in a list
    PDB_text=[x.strip() for x in open(traj)]
    
    # Start the numbering of the models
    model_number = 1
    new_file_text = ""
    
    # Dictionaries containing the calculated observables
    asaTotal={} # Dictionary for ASA
    energyTotal={} # Dictionary rosetta energy residue
    
    # Assign to each residue an empty list
    for i,aa in enumerate(pep_to_model):
        residue=aminoacids[aa]
        position=position_reference+i
        asaTotal[residue+str(position)]=[]
        energyTotal[residue+str(position)]=[]
        
    # Create folders in the peptide workspace
    os.system("mkdir asa energy")
    
    # Loop over the PDB lines to detect when a model is defined
    for line in PDB_text: 
        if line == "ENDMDL":
            # Save file with file number in name
            output_file = open("model_" + str(model_number) + ".pdb", "w")
            output_file.write(new_file_text.rstrip('\r\n')) #rstrip to remove trailing newline
            output_file.close()
    
            # Modify some nomenclatures to avoid issues with Rosetta
            os.system("sed -i 's/CD  ILE {}/CD1 ILE {}/g' model_{}.pdb".format(pepChain,pepChain,str(model_number)))
            os.system("sed -i 's/ASH/ASP/g' model_{}.pdb".format(str(model_number)))
            os.system("sed -i 's/GLH/GLU/g' model_{}.pdb".format(str(model_number)))
            os.system("sed -i 's/HID/HIS/g' model_{}.pdb".format(str(model_number)))
            os.system("sed -i 's/HIE/HIS/g' model_{}.pdb".format(str(model_number)))
            os.system("sed -i 's/HIP/HIS/g' model_{}.pdb".format(str(model_number)))
            
            # Running score interface
            os.system("{rosetta}/main/source/bin/InterfaceAnalyzer.default.linuxgccrelease -database {rosetta}/main/database/ -in:file:s model_{mod}.pdb -out:file:scorefile score_interface.sc -add_regular_scores_to_scorefile true -pack_input true -pack_separated true -interface A_B".format(rosetta=rosetta_path,mod=str(model_number)))
            
            # Get the residue number
            parser = PDBParser()
            structure = parser.get_structure('REF',"model_{mod}.pdb".format(mod=str(model_number)))
            model = structure[0]
            out_structure=model["A"]
            
            final_residue=0
            # Store the numbers of the residues that should be changed
            for residue in out_structure:
                final_residue=int(residue.get_full_id()[3][1])
            
            # Run over the output and store the energies for each residue
            for i,aa in enumerate(pep_to_model):
                residue=aminoacids[aa]
                position=position_reference+i
                bash = "grep _{} model_{}_0001.pdb".format(position+final_residue,str(model_number))
                output = subprocess.check_output(['bash','-c', bash])
                energyTotal[residue+str(position)].append(output.strip())
            
            # Delete temporal files
            os.system("rm model_{}_0001.pdb score_interface.sc".format(str(model_number)))
    
            # Calculate ASA
            p = PDBParser()
            structure = p.get_structure("MODEL", "model_{}.pdb".format(str(model_number)))
            model = structure[0]
            dssp = DSSP(model, 'model_{}.pdb'.format(str(model_number)), dssp='mkdssp')
            
            # Store the information for each amino acid
            for keys in list(dssp.keys()):
                if keys[0]==pepChain:
                    asaTotal[aminoacids[dssp[keys][1]]+str(keys[1][1])].append(str(dssp[keys][3]))
            
            # Save the model
            model_number += 1
            new_file_text = ""
        elif not line.startswith("MODEL"):
            data=line.split()
            if data[0]=="ATOM":
                if data[4] != compareChain:
                    new_file_text += '\n'
                    compareChain=data[4]    
            if data[0]=="ATOM" or data[0]=="TER":
                new_file_text += line + '\n'
    
    # Delete at the end all the split models
    os.system("rm model*.pdb {}".format(traj))
    
    # Return the dictionaries
    return asaTotal,energyTotal

##################################################################################

def obtain_averages(asaTotal,energyTotal,pep_to_model):
    """
    Based on the calculation of the observables calculate final averages and store the values in files
    
    Arguments:
    asaTotal -- Dictionary with information of the calculated ASA values
    energyTotal -- Dictionary with information of the calculated energy values
    pep_to_model -- Original sequence of the peptide to model
    
    Output:
    Files containing the single values, as well as a final report with the average values per peptide amino acid
    """
    aminoacids={"A":"ALA","D":"ASP","E":"GLU","F":"PHE","H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","G":"GLY",
                "N":"ASN","P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR","C":"CYS"}
    
    # File with the final averages
    final_averages=open("final_averages_{}.txt".format(pep_to_model),"w")
    final_averages.write("Amino_acid\tPosition\tAverage_ASA\tAverage_Energy\n")
    position_reference=1
    
    # Iterate over the peptide amino acids
    for i,aa in enumerate(pep_to_model):
        residue=aminoacids[aa]
        position=position_reference+i
        
        # Generate the files to store the calculated observables
        fileAsa=open("asa/asa_BR_{}{}_{}.txt".format(residue,position,pep_to_model),"w")
        fileEnergy=open("energy/energy_BR_{}{}_{}.txt".format(residue,position,pep_to_model),"w")
        
        # Write the values
        for value in asaTotal[residue+str(position)]:
            fileAsa.write(str(value)+"\n")
              
        for value in energyTotal[residue+str(position)]:
            fileEnergy.write(str(value)+"\n")
        
        # Close the files
        fileAsa.close()
        fileEnergy.close()
    
        # Open the ASA file to calculate the mean value
        asa_prev=[float(x.strip()) for x in open("asa/asa_BR_{}{}_{}.txt".format(residue,position,pep_to_model))]
        averageAsa=mean(asa_prev)
        
        # Open the energy file to calculate the mean value
        energy_prev=[float(x.strip().split()[1]) for x in open("energy/energy_BR_{}{}_{}.txt".format(residue,position,pep_to_model))]
        averageEnergy=mean(energy_prev)
        
        # Save the line
        final_averages.write("{}\t{}\t{}\t{}\n".format(residue,position,averageAsa,averageEnergy))
    
    final_averages.close()

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
    parser.add_argument('-p', dest='path', action='store',required=True,
                        help='Path of the structure that will be used to run the analysis')
    parser.add_argument('-s', dest='sequence', action='store', required=True,
                        help='Sequence of the peptide that will be modelled and sampled')
    parser.add_argument('-c', dest='chain', action='store', required=True,
                        help='Chain identifier of the peptide in the structure')
    parser.add_argument('-r', dest='rosetta', action='store', default="rosetta_src_2016.32.58837_bundle",
                        help='Version of Rosetta that will be implemented')
    args = parser.parse_args()
    
    ####################################################################################
    # Assignment of parameters
    model=args.path
    pep2_final=args.sequence
    pepChain=args.chain
    rosetta_version=args.rosetta

    # Get the sequence of the peptide used as template
    aminoacids={"ALA":"A","ASP":"D","GLU":"E","PHE":"F","HIS":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","GLY":"G",
                "ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y","CYS":"C"}
    pep1_final=""
    parser = PDBParser()
    reference = parser.get_structure('REF',model)
    for chain in reference[0]:
        if chain.get_id()==pepChain:
            for residue in chain:
                seq=aminoacids[residue.get_resname()]
                pep1_final=pep1_final+str(seq)
    
    # Get the rosetta path
    bash = "locate -b {} | head -n1".format(rosetta_version)
    rosetta_path = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
    
    # Create folders in case they do not exist
    os.system("mkdir dynamic")
    os.system("mkdir dynamic/models")
    os.system("mkdir dynamic/observables")
    
    print("Starting modelling of peptide {} ...".format(pep2_final))
    
    # Generate the list of mutations and shuffle the order
    mutList=generateMutations(pep1_final,pep2_final)
    shuffle(mutList)
    
    # Run the protocol to model the new peptide
    modelling(mutList,model,pep2_final,pepChain,rosetta_path)
    
    # Sample and refine the complex
    run_sampling(rosetta_path,pep2_final)
    
    # Calculate the observables required for the analysis
    asa_data,energy_data=calculate_observables(pep2_final,pepChain,rosetta_path)
    
    # Obtain the final averages and the files with the data
    obtain_averages(asa_data,energy_data,pep2_final)
    
    # Move files
    os.system("mv {}_* dynamic/models".format(pep2_final))
    os.system("mv asa dynamic/observables/asa_{}".format(pep2_final))
    os.system("mv energy dynamic/observables/energy_{}".format(pep2_final))
