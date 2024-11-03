#!/bin/bash

for i in {1..9}; do  wget "http://ftp.ensembl.org/pub/release-105/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.chromosome.$i.fa.gz";
                       done;

for i in {10..25}; do  wget "http://ftp.ensembl.org/pub/release-105/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.chromosome.$i.fa.gz";
                       done;

wget http://ftp.ensembl.org/pub/release-105/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.chromosome.MT.fa.gz
