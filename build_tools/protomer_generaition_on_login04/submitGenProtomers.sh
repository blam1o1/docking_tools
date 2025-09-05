#!/bin/bash

export SBASE=/lustre/fs6/lyu_lab/scratch/blam/soft/docking_tools/build_tools/protomer_generaition_on_login04

function exists {
        env_name=$1
        desc=$2
        if [ -z "${!env_name}" ]; then
                echo "expected input arg: $env_name"
                echo "arg description: $desc" 
                failed=1
        fi
}
function exists_warning {
        env_name=$1
        desc=$2
        default=$3
        if [ -z "${!env_name}" ]; then
                echo "optional env arg missing: $env_name"
                echo "arg description: $desc"
                echo "defaulting to $default"
                export $env_name="$default"
        fi
}


exists DIR_OF_SMILES "directory containing smile files"
exists_warning PH "pH to protonate protomers at" 7.4

export WORKDIR=$PWD

echo "Starting protomer gen process.. creating protomer_gen directory"
mkdir $WORKDIR/protomer_gen
cd $WORKDIR/protomer_gen
export PROT_DIR=$PWD

#Split smiles files into batches corresponding to num available nodes:
num_nodes=$(sinfo -N 2>/dev/null | grep lyu | grep -e idle -e mix -e alloc | wc -l)
echo "Total nodes available: $num_nodes"
total_lines=$(cat $DIR_OF_SMILES/*smi | wc -l)
#total_lines=$(wc -l $DIR_OF_SMILES/*smi)
lines_per=$((total_lines / num_nodes +1 )) #+1 ensures no remainder that will create another job that there is no node for
echo "Splitting smiles into $lines_per line  batches"
cat $DIR_OF_SMILES/*smi | split -l $lines_per - genProt__
ls $PROT_DIR/genProt__* > glist


#Simplest option would be to just do it as whole batches per node. No need to worry about node availablility changing with one long job

for batch in $(cat glist);do

        cd $PROT_DIR

        unset NODE_NUMBER
        unset INPUT_SMILES
        unset BASENAME

        basename=$(basename "$batch")
        echo "$basename"
        mkdir "dir_$basename"
        mv $basename  "dir_$basename/$basename.smi"

        nodes_in_use=$(awk '{print $NF}' <(squeue | grep lyu | grep protP | grep -v -e None -e Priority) <(grep node $PROT_DIR/dir*/*.sh) | sort -u)
        echo 

        echo -e "Nodes in use:\n$nodes_in_use\n"
        
        if [ -z "$nodes_in_use" ]; then
                nodes_available=$(sinfo -N | grep lyu | grep -v -e comp -e down -e drain -e drng | awk '{print $1}' | sort -u)
        else 
                nodes_available=$(sinfo -N | grep lyu | grep -v -e comp -e down -e drain -e drng | grep -v -f <(echo "$nodes_in_use") |  awk '{print $1}' | sort -u)
        fi

        echo -e "Available Nodes:\n$nodes_available\n "        

        node_submit=$(head -1 <(echo "$nodes_available"))
        echo "Submitting on $node_submit"

        export NODE_NUMBER=$node_submit
        export INPUT_SMILES=$PROT_DIR/dir_$basename/$basename.smi
        export BASENAME=$basename

        #Generating submission script and submitting protomer generation job
        cd "dir_$basename"
        bash $SBASE/create-submit.sh
        echo -e "Generated submission script: genProtomerOnNode.sh"
        echo -e "If job fails you are able to resubmit using this script...\n"

        echo -e "Submitting job and then sleeping for 5s...\n"
        sbatch genProtomerOnNode.sh
        sleep 5

        echo -e "----------------------------------------------------\n"

done










