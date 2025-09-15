HGT preprocess workflow

*collate_tables.py*
Inputs: assembly.fasta 
		    assembly.gfa
		    gplas results.tab
		    abricate concatenated output.tab
		    checkm concatenated report.tsv

Libraries: 	pandas, numpy, os, re, subprocess

Folder setup and files
      creates assembly output, fasta and gfa input folders
		  linearises fasta files in fasta folder (makes splitting assemblies easier)
		  will also need to create abricate, gplas and checkm input folders and move files in
  
Loading data and dataframes
      creates sample_names list and dataframe from original assembly files, will eventually change to reading names from a shared inputs folder in snakemake
		  loads gplas csv, finds matching nodes from assembly files and creates a new column with the node number, creates a new column with sample names. 
      Sample names need to be added to the dataframe created for each individual sample first before concatenating  
			Loads abricate files and cleans up sample column
			Loads checkm csv
			create list of unique plasmid nodes to prevent duplicating fasta sequence (only works if assumptions about gplas output are true)
      Creating table loops through samples, finds nodes from unique_nodes which match the sample name then creates a frame of rows from the abricate file which match to the sample name, 
      finds rows from that frame which match to the node. These are plasmid associated genes. 
      Separately checks for rows from the abricate file which do not match any nodes from gplas unique_nodes but do match the sample name, 
      this step contains a check if node list is empty which is then also used to fill in the “has_plasmid” column. These are chromosome associated genes
			repeats loop for different abricate databases
 			Collates relevant columns from loaded and generated dataframes into one summary.csv

Splitting assemblies	creates sample output folders within assembly/
			finds unique nodes using the same method used for abricate, feeds these as arguments to create_plasmids.sh 
      this highlights the nodes one at a time in the linearised fasta assembly and copies these to sample_plasmids.fasta,
      also adds each node to an extended regex expression which is fed into create_chromosomes.sh as an argument,
      this finds all plasmid nodes at the same time in the fasta assembly and copies all but the plasmid nodes into sample_chromosome.fasta
					
