#!/bin/bash

# Brandon L
#
# Submission script for generating protomers for smiles.
# DIR_OF_SMILES - replace with path to directory containing all smiles files 
# that you would like to protonate
#
# Written to accomodate parallelization of protomer generaiton on login04.
#

export DIR_OF_SMILES=/lustre/fs6/lyu_lab/scratch/blam/work/vacht/chemspace_library/chloroacetamides/protonate

bash /lustre/fs6/lyu_lab/scratch/blam/soft/docking_tools/build_tools/protomer_generaition_on_login04/submitGenProtomers.sh
