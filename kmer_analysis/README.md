# kmer-analysis
Includes Perl and R scripts that perform various kmer analyses from FASTA files. 

### Programs required

Programs described here are running using the following versions for: 
- Perl: (v5.16.3)
- R:  3.6.0 

For the rapid kmer assesment, use the following C code from lmj: https://github.com/lmjakt/R_c_plugins

### Input data format

Use files in FASTA format. The Identifiers are not used, only the sequence content. Please make sure that the file does not contain any Ns or lowercase letters as the C code cannot operate with multisemantic letters.

### Contents

**kmer_analysis_functions.R** : Includes functions that create kmer counts and plot either density plots of kmers against two different datasets or a density plot of kmer counts and their reverse complement. 

Please refer to the Wiki for more information on the R functions.

**apply_kmer_symmetry.R** This is an R script that can be used to generate symmetry plot from a FASTA file sourcing the kmer analysis functions. The command takes four arguements: the FASTA file, the output name, the low density color and the high density color. 

Syntax example: 

`./apply_kmer_symmetry.R fasta_file.fa outputname "blue" "red"`

The command will complain if one or more of those arguments are not included. 
Pixel number (the number of rectangles plotted in each side) can only be changed by changing the source code at the  `pix=` part . 
The ouptput is a jpeg that includes a set of looged symmetry graphs with kmer sizes from 4 to 9. 

An example of the kmer symmetry graphs using the human genome: 

![image](https://user-images.githubusercontent.com/29709382/164441022-577b6810-5f20-4270-a40e-9d4eebd3ff01.png)
