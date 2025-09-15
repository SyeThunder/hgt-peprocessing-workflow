#!/bin/bash

#create folders
mkdir assemblies
mkdir fasta
mkdir gfa

#move gfa files
cd /home/stephen/Documents/genomes/unicycler/assemblies/$1/processed
for err in $(ls)
    do cp ${err}/assembly_graph_with_scaffolds.gfa /home/stephen/Documents/Code/pandas/gfa/
    cd /home/stephen/Documents/Code/pandas/gfa/
    mv assembly_graph_with_scaffolds.gfa ${err}.gfa
    cd /home/stephen/Documents/genomes/unicycler/assemblies/$1/processed
done

#move and linearise fastas
for err in $(ls)
	do awk '{if(NR==1) print $0; else if($0 ~ /^>/) print "\n"$0; else printf $0}' ${err}/scaffolds.fasta > /home/stephen/Documents/Code/pandas/fasta/${err}_linear.fasta
	done
	

