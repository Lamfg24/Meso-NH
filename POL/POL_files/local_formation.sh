#!/bin/bash
#MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
#MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
#MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
#MNH_LIC for details. version 1.
#SBATCH -J MNH1
#SBATCH -n 1                                                                                                                                                                                                                                                             
#SBATCH -o Examples.eo%j   #                                                                                                                                                                                                                                             
#SBATCH -e Examples.eo%j   #                                                                                                                                                                                                                                             
#SBATCH -t 00:15:00    # time limit                                                                                                                                                                                                                                      
#SBATCH --mem=100000                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                         
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

cd local_formation/
pwd
./run_prep_ideal_case_xyz

