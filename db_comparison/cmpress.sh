!#/bin/bash

infernalpress="/home/ama/bin/infernal-1.1.4/src/cmpress"
files="*cm"

for i in $files;
do
    echo "For file: $i"
    $infernalpress $i

done;
