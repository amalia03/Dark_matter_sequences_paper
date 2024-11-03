#!/bin/bash
#OrfPredictor.pl fastaFile BLASTX bFlag strand E-value output

filenames="$(ls *lengths.fa | perl -pe 's/([0-9]+[A-Z]+).*/$1/')"

#echo $filenames

for file in $filenames;
do
#echo $file

fa_file="$(ls $file*lengths.fa)"

echo "Starting ORF analysis for sample $fa_file"

mkdir $file
cd $file

perl /home/ama/bin/OrfPredictor/OrfPredictor.pl ../$fa_file blank.txt 0 both 1 $file

#echo $file
cd ..

done
