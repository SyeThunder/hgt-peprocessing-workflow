#!/bin/bash

#create folders
mkdir assemblies
mkdir fasta
mkdir gfa
mkdir gplas
mkdir abricate
mkdir checkm

#move gfa files
for err in $(ls /home/stephen/Documents/genomes/unicycler/assemblies/$1/processed)
    do cp /home/stephen/Documents/genomes/unicycler/assemblies/$1/processed/${err}/assembly_graph_with_scaffolds.gfa ./gfa
    mv gfa/assembly_graph_with_scaffolds.gfa gfa/${err}.gfa
done

#move and linearise fastas
for err in $(ls /home/stephen/Documents/genomes/unicycler/assemblies/$1/processed)
	do awk '{if(NR==1) print $0; else if($0 ~ /^>/) print "\n"$0; else printf $0}' /home/stephen/Documents/genomes/unicycler/assemblies/$1/processed/${err}/scaffolds.fasta > fasta/${err}_linear.fasta
	done

#move gplas files
for err in $(ls /home/stephen/Documents/genomes/gplas2/$1/results | grep results)
	do cp /home/stephen/Documents/genomes/gplas2/salmonella_typhi_ciprofloxacin_h58/results/${err} ./gplas
	done
	
#move checkm files
cp /home/stephen/Documents/genomes/checkm2/$1/report.tsv ./checkm/

#move abricate files
cp /home/stephen/Documents/genomes/abricate/$1/* ./abricate/
	

