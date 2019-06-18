#! /usr/bin/env bash

module load biopython
module load python/3.7.3

./scripts/chroplasmitor.py -g samples/complete_genomes/*.fna -o samples