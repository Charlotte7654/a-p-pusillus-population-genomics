#!/bin/bash
#PBS -N 05_popstats_full.R
#PBS -l select=1:ncpus=8:mem=8GB
#PBS -l walltime=23:00:00
#PBS -m be
#PBS -M 25405942@sun.ac.za

cd $PBS_O_WORKDIR
module load app/R/4.3.2

#in same 00_scripts directory as the script is being submitted 
#I believe cd $PBS_O_WORKDIR means make directory of submission working dir 
R --file=05_popstats_full.R
