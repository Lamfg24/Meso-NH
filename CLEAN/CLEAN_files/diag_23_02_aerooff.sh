#!/bin/bash
#MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence                                                                                     
#MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt                                                                                         
#MNH_LIC for details. version 1.                                                                                                                                 
#SBATCH -J MNH1                                                                                                                                                  
#SBATCH -N 1                                                                                                                                                     
#SBATCH -n 1                                                                                                                                                       
#SBATCH -o Examples.eo%j   #                                                                                                                                       
#SBATCH -e Examples.eo%j   #                                                                                                                                       
#SBATCH -t 2:30:00    # time limit                                                                                                                               
#SBATCH --mem=150000                                                                                                                                                                                      
ulimit -c 0                                                                                                                                                                                              
ulimit -s unlimited                                                                                                                                                                                      
# Arret du job des la premiere erreur
set -e
# Nom de la machine
hostname 
# Echo des commandes
echo $SLURM_NTASKS
. ../../conf/profile_mesonh-LXifort-R8I4-MNH-V5-4-2-ECRAD-DACCIWA-MPIINTEL-O2
export MONORUN="mpirun -np 1 "
export MPIRUN="mpirun "
export POSTRUN="time "

cd clean_23_02_aerooff/
pwd
./run_diag_xyz

