#!/bin/bash

grep -A1 -m1 "$1" fasta/$2_linear.fasta >> assemblies/$2/$2_plasmids.fasta


