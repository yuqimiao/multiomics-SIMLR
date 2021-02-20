#!/bin/sh
#$ -cwd -S /bin/bash
#$ -l mem=60G
#$ -l time=36:00:00
#$ -N R-SIMLR-3data-s4-v10 -j y 

SCRIPT=$1

R=/nfs/apps/R/3.6.0/bin/R
R_LIBS_USER=~/local/R/hpc:/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs

${R} --vanilla < ${SCRIPT}