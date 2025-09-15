#run after creating folder of sample assembly files 
import pandas as pd
import numpy as np
import os
import re
import subprocess

#setup folders and files
#subprocess.run(["bash", "setup.sh", "salmonella_typhi_ciprofloxacin_h58"])
#sample names
#will have to be adapted for snakemake samples input
sample_names = [err for err in os.listdir("/home/stephen/Documents/genomes/unicycler/assemblies/salmonella_typhi_ciprofloxacin_h58/processed")]
sample_names = sorted(sample_names)
sample_df = pd.DataFrame(sample_names, columns=["sample"])

#loading and prepping dataframes

#gplas2
#equivalent of for files in $(ls | grep results) do ...
#will have to be adapted for snakemake samples input
#could change this to a bash script in setup.sh that moves results files in
results_tabs = [file for file in os.listdir("/home/stephen/Documents/genomes/gplas2/salmonella_typhi_ciprofloxacin_h58/results")
    if "results" in file
]  
#creates a dictionary called gplas_frames
gplas_frames = {}

#adding data to that dictionary
#will have to be adapted for snakemake samples input
for tab in results_tabs:
    notab = os.path.splitext(tab)[0]
    filepath = os.path.join("/home/stephen/Documents/genomes/gplas2/salmonella_typhi_ciprofloxacin_h58/results/", tab)
    gplas_frames[notab] = pd.read_csv(filepath, sep= " ")
    
#print(gplas_frames)

for k,frame in gplas_frames.items():
    frame["Sample"] = k.split("_")[0]
    gplas_frames[k]=frame
    #print(pd.DataFrame.from_dict(frame).iloc[0]["Sample"])
    
    
#with open("gfa/salmonella_typhi_ciprofloxacin_h58/ERR023622.gfa", "r") as f:
#    lines = f.readlines()
def extract_node(number):
    match = [line for line in lines if str(number) in line]
    if not match:
        return None
    
    node_line = match[-1]
    #selecting the last line that matches the contig number
    node = re.search(r"NODE_(.*?)_", node_line)
    if node:
        return node.group(1)
    return None

#test = pd.read_csv("/home/stephen/Documents/genomes/gplas2/salmonella_typhi_ciprofloxacin_h58/results/ERR023622_results.tab", sep= " ")


#test["node"] = test["number"].apply(extract_node)
#print(test["node"])
#will have to be adapted for snakemake samples input
for k, frame in gplas_frames.items():
    gplas_frames[k] = frame
    df = pd.DataFrame.from_dict(frame).iloc[0]["Sample"]
    with open("gfa/salmonella_typhi_ciprofloxacin_h58/" + df + ".gfa", "r") as f:
        lines = f.readlines()

    def extract_node(number):
        match = [line for line in lines if str(number) in line]
        if not match:
          return None
    
        node_line = match[-1]
        #selecting the last line that matches the contig number
        node = re.search(r"NODE_(.*?)_", node_line)
        if node:
            return node.group(1)
        return None

    frame["node"] = frame["number"].apply(extract_node)

#with open("/home/stephen/Documents/genomes/unicycler/assemblies/salmonella_typhi_ciprofloxacin_h58/processed") for k, frame in gplas_frames.items():

gplas = pd.concat(gplas_frames.values(), ignore_index = True)
del gplas_frames

gplas = gplas.sort_values("Sample")
#print(gplas)
gplas.to_csv("gplas.csv", index = False)

#abricate will have to be adapted for snakemake samples input
ncbi = pd.read_csv("/home/stephen/Documents/genomes/abricate/salmonella_typhi_ciprofloxacin_h58/ncbi.tab", sep="\t")
resfinder = pd.read_csv("/home/stephen/Documents/genomes/abricate/salmonella_typhi_ciprofloxacin_h58/resfinder.tab", sep="\t")
plasmidfinder = pd.read_csv("/home/stephen/Documents/genomes/abricate/salmonella_typhi_ciprofloxacin_h58/plasmidfinder.tab", sep="\t")
ncbi["#FILE"] = ncbi["#FILE"].str.replace("/scaffolds.fasta", "")
resfinder["#FILE"] = resfinder["#FILE"].str.replace("/scaffolds.fasta", "")
plasmidfinder["#FILE"] = plasmidfinder["#FILE"].str.replace("/scaffolds.fasta", "")

#checkm will have to be adapted for snakemake samples input
checkm = pd.read_csv("/home/stephen/Documents/genomes/checkm2/salmonella_typhi_ciprofloxacin_h58/report.tsv", sep="\t")

unique_nodes = gplas.drop_duplicates("node")
#print(unique_nodes)

"""
#generating plasmid/chromosome arg columns
#these will iterate through each sample, if no match is found it returns none
#this works for no match to the args or no match to the gplas plasmid contigs
#has to work differently for chromosome contigs where if no match is found it will still check for matches to the sample name in the abricate file
#These are working as intended right now, each creates a dataframe of genes detected by an individual abricate method whether gplas detects a plasmid or not
"""
#ncbi
ncbi_plasmid = []
for s in sample_names:
    #find unique nodes from gplas which match the sample name
    sample_nodes = unique_nodes.loc[unique_nodes["Sample"] == s].copy()
    #add "NODE_" to ensure correct match with abricate file
    sample_nodes["node"] = "NODE_" + sample_nodes["node"]
    #print(sample_nodes)
    
    #find rows from abricate file where sample name match to current gplas sample
    all_arg_sample = (ncbi[ncbi["#FILE"].isin(sample_nodes["Sample"])])
    #find rows which match to node
    matches = "|".join(sample_nodes["node"].astype(str).apply(re.escape))
    arg_plasmid_sample = all_arg_sample[all_arg_sample["SEQUENCE"].str.contains(matches)]
    #creates single string containing all genes
    arg_plasmid_str = ", ".join(arg_plasmid_sample["GENE"].astype(str))
    #appends to a cell in dictionary
    ncbi_plasmid.append({"ncbi_plasmid_genes": arg_plasmid_str})

ncbi_plasmid_df = pd.DataFrame(ncbi_plasmid)

ncbi_chromosome = []
hasplasmid = []

for s in sample_names:
    sample_nodes_1 = unique_nodes.loc[unique_nodes["Sample"] == s].copy()
    sample_nodes_1["node"] = "NODE_" + sample_nodes_1["node"]
    if sample_nodes_1.empty: 
        #find rows from abricate file where sample name matches abricate[~FILE]
        noplas_arg = (ncbi[ncbi["#FILE"] == s])
        #create string of all args
        noplas_arg_str = ", ".join(noplas_arg["GENE"].astype(str))
        ncbi_chromosome.append({"ncbi_chromosome_genes": noplas_arg_str})
        hasplasmid.append({"has_plasmid": "no"})

    else:
        #find rows from abricate file where sample name match to current gplas sample
        all_arg_sample_1 = (ncbi[ncbi["#FILE"].isin(sample_nodes_1["Sample"])])
        #find rows which do not match to node
        matches_1 = "|".join(sample_nodes_1["node"].astype(str).apply(re.escape))
        #create string of arg non-matches
        arg_chromosome_sample = all_arg_sample_1[~all_arg_sample_1["SEQUENCE"].str.contains(matches_1)]
        arg_chromosome_str = ",".join(arg_chromosome_sample["GENE"].astype(str))
        ncbi_chromosome.append({"ncbi_chromosome_genes": arg_chromosome_str})
        hasplasmid.append({"has_plasmid": "yes"})

ncbi_chromosome_df = pd.DataFrame(ncbi_chromosome)
hasplasmid_df = pd.DataFrame(hasplasmid)

#Resfinder
resfinder_plasmid = []
for s in sample_names:
    #find unique nodes from gplas which match the sample name
    sample_nodes = unique_nodes.loc[unique_nodes["Sample"] == s].copy()
    #add "NODE_" to ensure correct match with abricate file
    sample_nodes["node"] = "NODE_" + sample_nodes["node"]
    
    #find rows from abricate file where sample name match to current gplas sample
    all_arg_sample = (resfinder[resfinder["#FILE"].isin(sample_nodes["Sample"])])
    #find rows which match to node
    matches = "|".join(sample_nodes["node"].astype(str).apply(re.escape))
    arg_plasmid_sample = all_arg_sample[all_arg_sample["SEQUENCE"].str.contains(matches)]
    #creates single string containing all genes
    arg_plasmid_str = ", ".join(arg_plasmid_sample["GENE"].astype(str))
    #appends to a cell in dictionary
    resfinder_plasmid.append({"resfinder_plasmid_genes": arg_plasmid_str})

resfinder_plasmid_df = pd.DataFrame(resfinder_plasmid)

resfinder_chromosome = []

for s in sample_names:
    sample_nodes_1 = unique_nodes.loc[unique_nodes["Sample"] == s].copy()
    sample_nodes_1["node"] = "NODE_" + sample_nodes_1["node"]
    if sample_nodes_1.empty: 
        #find rows from abricate file where sample name matches abricate[~FILE]
        noplas_arg = (resfinder[resfinder["#FILE"] == s])
        print(noplas_arg)
        #create string of all args
        noplas_arg_str = ", ".join(noplas_arg["GENE"].astype(str))
        resfinder_chromosome.append({"resfinder_chromosome_genes": noplas_arg_str})

    else:
        #find rows from abricate file where sample name match to current gplas sample
        all_arg_sample_1 = (resfinder[resfinder["#FILE"].isin(sample_nodes_1["Sample"])])
        #find rows which do not match to node
        matches_1 = "|".join(sample_nodes_1["node"].astype(str).apply(re.escape))
        #create string of arg non-matches
        arg_chromosome_sample = all_arg_sample_1[~all_arg_sample_1["SEQUENCE"].str.contains(matches_1)]
        arg_chromosome_str = ",".join(arg_chromosome_sample["GENE"].astype(str))
        #print(arg_chromosome_str)
        resfinder_chromosome.append({"resfinder_chromosome_genes": arg_chromosome_str})

resfinder_chromosome_df = pd.DataFrame(resfinder_chromosome)

#plasmidfinder
plasmidfinder_plasmid = []
for s in sample_names:
    #find unique nodes from gplas which match the sample name
    sample_nodes = unique_nodes.loc[unique_nodes["Sample"] == s].copy()
    #add "NODE_" to ensure correct match with abricate file
    sample_nodes["node"] = "NODE_" + sample_nodes["node"]
    
    #find rows from abricate file where sample name match to current gplas sample
    all_arg_sample = (plasmidfinder[plasmidfinder["#FILE"].isin(sample_nodes["Sample"])])
    #find rows which match to node
    matches = "|".join(sample_nodes["node"].astype(str).apply(re.escape))
    arg_plasmid_sample = all_arg_sample[all_arg_sample["SEQUENCE"].str.contains(matches)]
    #creates single string containing all genes
    arg_plasmid_str = ", ".join(arg_plasmid_sample["GENE"].astype(str))
    #appends to a cell in dictionary
    plasmidfinder_plasmid.append({"plasidfidner_plasmid_genes": arg_plasmid_str})

plasmidfider_plasmid_df = pd.DataFrame(plasmidfinder_plasmid)

plasmidfinder_chromosome = []

for s in sample_names:
    sample_nodes_1 = unique_nodes.loc[unique_nodes["Sample"] == s].copy()
    sample_nodes_1["node"] = "NODE_" + sample_nodes_1["node"]
    if sample_nodes_1.empty: 
        #find rows from abricate file where sample name matches abricate[~FILE]
        noplas_arg = (plasmidfinder[plasmidfinder["#FILE"] == s])
        print(noplas_arg)
        #create string of all args
        noplas_arg_str = ", ".join(noplas_arg["GENE"].astype(str))
        plasmidfinder_chromosome.append({"plasmidfinder_chromosome_genes": noplas_arg_str})

    else:
        #find rows from abricate file where sample name match to current gplas sample
        all_arg_sample_1 = (plasmidfinder[plasmidfinder["#FILE"].isin(sample_nodes_1["Sample"])])
        #find rows which do not match to node
        matches_1 = "|".join(sample_nodes_1["node"].astype(str).apply(re.escape))
        #create string of arg non-matches
        arg_chromosome_sample = all_arg_sample_1[~all_arg_sample_1["SEQUENCE"].str.contains(matches_1)]
        arg_chromosome_str = ",".join(arg_chromosome_sample["GENE"].astype(str))
        #print(arg_chromosome_str)
        plasmidfinder_chromosome.append({"plasmidfinder_chromosome_genes": arg_chromosome_str})

plasmidfinder_chromosome_df = pd.DataFrame(plasmidfinder_chromosome)



#creating summary dataframe

summary = pd.concat([
    sample_df["sample"],
     checkm["Genome_Size"],
     checkm["Completeness"],
     checkm["Contamination"],
     ncbi_plasmid_df,
     ncbi_chromosome_df,
     resfinder_plasmid_df,
     resfinder_chromosome_df,
     plasmidfider_plasmid_df,
     plasmidfinder_chromosome_df,
     hasplasmid_df
     #"plasmid_names":
], axis= 1, ignore_index= False)

summary.to_csv("summary.csv", index = False)


#separating chromosome and plasmid assemblies

#create sample output folders
for name in sample_names:
    os.makedirs("assemblies/" + name)
    #subprocess.run(["bash", "mkdir", "assemblies/" + s])
#create assemblies
for s in sample_names:

    #find unique nodes from gplas output which match sample name
    #ie finding nodes identified as plasmid derived
    sample_nodes_2 = unique_nodes.loc[unique_nodes["Sample"] == s].copy()
    #ensures no match with a coverage value in fasta file
    sample_nodes_2["node"] = ">NODE_" + sample_nodes_2["node"]
    plasmid_nodes = str()
    for node in sample_nodes_2["node"]:
        subprocess.run(["bash", "create_plasmids.sh", node, s])
        plasmid_nodes += node + "/|"
    #removes last "/|"
    plasmid_nodes = plasmid_nodes[:-1]
    plasmid_nodes = plasmid_nodes[:-1]
    subprocess.run(["bash", "create_chromosomes.sh", plasmid_nodes, s])
    chromosome_empty = os.path.getsize("assemblies/" + s + "/" + s + "_chromosome.fasta")
    if chromosome_empty == 0:
        subprocess.run(["bash", "fill_chromosome.sh", s])
    else:
        pass
        
#works only if there is at least 1 plasmid derived contig
#read output chromosome files into python, check if null, if null then copy ERR.fasta to ERR_chromosome.fasta

#last thing was changing output directory of assemblies to a folder with the sample aname but this seemed to break the output gen