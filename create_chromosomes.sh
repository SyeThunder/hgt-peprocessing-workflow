#!/bin/bash

grep -A1 -v "$1" fasta/$2_linear.fasta >> assemblies/$2/$2_chromosome.fasta 


