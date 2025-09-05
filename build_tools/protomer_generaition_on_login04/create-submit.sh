#!/bin/bash

# Brandon L
# Setup for subdirectories where protonation jobs will run

node=$NODE_NUMBER
ph=$PH
smiles=$INPUT_SMILES
dir=$BASENAME
partition=$(scontrol show node $node | grep Partitions | awk -F'=' '{print $2}')
cat <<EOF > genProtomerOnNode.sh
#!/bin/bash
#SBATCH --job-name=protP
#SBATCH -p $partition
#SBATCH -w $node
#SBATCH --output=GenProt.out
#SBATCH --error=GenProt.err

ph=$PH

export WORKDIR=\$PWD
source /lustre/fs6/lyu_lab/scratch/lyu_soft/jchem/current/env.sh
source /lustre/fs6/lyu_lab/scratch/lyu_soft/java11/env.sh
source /lustre/fs6/lyu_lab/scratch/lyu_soft/dock/versions/dock37/DOCK-3.7-trunk/env.sh

split -l 50000 "$INPUT_SMILES" protomerGen_fold_

ls protomerGen_fold_* > plist


for file in \$(cat plist) ;do
    echo -e "Processing smiles for \$file\n"
    echo -e "----------------------------------------------------\n"
    cd \$WORKDIR
    mkdir "dir_\$file"
    cp \$file "dir_\$file/\${file}.smi"
    
    cd "dir_\$file"
    echo \$PWD


    bash /lustre/fs6/lyu_lab/scratch/lyu_soft/dock/versions/dock37/DOCK-3.7-trunk/analysis/new_DUDE_SCRIPTS/generate_protomers_from_smiles.sh -H \$ph \$PWD/\${file}.smi

done

EOF
