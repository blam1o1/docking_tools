#!/bin/bash

# Brandon L
# Combine results of parallelized protomer generation

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <protomer_gen_dir>"
    exit 1
fi

function dominantState {
    input="$1"-protomers-expanded.ism
    return="$1"_protomers.smi
    awk '$3 == 0 {print $1 " " $2}' $input > $return
}

function rename_stereoCenters {
    input="$1"_protomers.smi
    return="$1"_protomers_extended.smi
    awk '{if ($2 in id_array) {id_array[$2]++; $2 = $2 "." id_array[$2]} else {id_array[$2] = 0; $2 = $2 ".0"}} 1'  $input > $return

}


WORKDIR=$PWD
DIR_TO_COMBINE="$WORKDIR/$1"

for dir in $(ls -d "$DIR_TO_COMBINE"/dir*/dir*); do
    echo -e "Combining smiles in $(basename "$dir"):"
    name=$(basename "$dir"/*smi | awk -F'.' '{print $1}')
    dominantState "$dir"/working/protonate/"$name"
    rename_stereoCenters "$dir"/working/protonate/"$name"
    
    cat "$dir"/working/protonate/"$name"_protomers_extended.smi >> $WORKDIR/protomers_extended.smi

done





