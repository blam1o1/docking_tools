#! /bin/bash

#Brandon L.
# Repair Incomplete OUTDOCK written by broken DOCK3/dock64 binary 
# dirlist_repair - paths to jobs subdirectories of output directory

repair=./repair_outdock.py
WORKDIR=$PWD



source /lustre/fs6/lyu_lab/scratch/blam/soft/miniconda3/etc/profile.d/conda.sh
conda activate analysis_env

#Iterate over job paths:
for dir in `cat dirlist_repair`;do
    echo "$dir"
    cd $dir
    mv OUTDOCK OUTDOCK_orig
    input_file_path=$dir/OUTDOCK_orig
    save_file_path=$dir/OUTDOCK
    python "$repair" "$input_file_path" "$save_file_path"
    cd "$WORKDIR"
done
    
