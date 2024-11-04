#!/bin/bash

## run consecutive scans so that we can control the number of threads used

analyses="Pfam SUPERFAMILY Gene3D ProSiteProfiles PRINTS PANTHER"

for appl in $analyses
do
    base=${1/%.fa/}_$appl
#    echo "file $appl => $base"
    nice interproscan.sh -i $1 -dp -f tsv,gff3 -appl $appl -goterms -b $base
#    nice interproscan.sh -i $1 -dp -f tsv,gff3 -appl $appl -goterms -b $base
done

## note the following timing information (for a 100 sequences)
## Pfam:             10:29:04:897 -> 10:30:13:522   1:09
## SUPERFAMILY:      10:30:15:453 -> 10:31:36:761   1:21
## Gene3D:           10:31:38:718 -> 10:34:04:109   2:30
## ProSiteProfiles:  10:34:06:099 -> 10:35:08:762   1:02
## PRINTS-42.0:      10:35:10:776 -> 10:36:05:022   0:55
##
## but in each case it takes close to 16 seconds to load the
## java application. I should probably make the fasta files
## a bit bigger. Increase to 500.
##
## Also note that most of the analyses use 1 thread by default
## but that this can be changed in interproscan.properties
##
## I've updated this so that a single thread should be used by
## all (prositeprofiles was using 4)
## But the thread usage is still rather inconsistent.
