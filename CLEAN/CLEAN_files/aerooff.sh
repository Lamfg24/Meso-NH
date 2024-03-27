#!/bin/bash
#MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence                                                                                     
#MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt                                                                                         
#MNH_LIC for details. version 1.              
#SBATCH -J MNH1                                                                                                                                                  
#SBATCH -N 50                                                                                                                                                    
#SBATCH -n 3600                                                                                                                                                     
#SBATCH -o Examples.eo%j   #                                                                                                                                       
#SBATCH -e Examples.eo%j   #                                                                                                                                       
#SBATCH -t 6:00:00    # time limit                                                                                                                               
#SBATCH --mem=50000         

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

#cd clean_23_02_aerooff/
#pwd
#./run_mesonh_xyz


#cd ../
#cp clean_23_02_aerooff/CLEA1.1.AEROF.004.nc clean_03_05_aerooff/.
#cp clean_23_02_aerooff/CLEA1.1.AEROF.004.des clean_03_05_aerooff/.

#cd clean_03_05_aerooff/
#pwd
#./run_mesonh_xyz



#cd ../
#cp clean_03_05_aerooff/CLEA2.1.AEROF.004.nc clean_06_08_aerooff/.
#cp clean_03_05_aerooff/CLEA2.1.AEROF.004.des clean_06_08_aerooff/.

#cd clean_06_08_aerooff/
#pwd
#./run_mesonh_xyz


#cd ../
#cp clean_06_08_aerooff/CLEA3.1.AEROF.004.nc clean_09_11_aerooff/.
#cp clean_06_08_aerooff/CLEA3.1.AEROF.004.des clean_09_11_aerooff/.

#cd clean_09_11_aerooff/
#pwd
#./run_mesonh_xyz


#cd ../
#cp clean_09_11_aerooff/CLEA4.1.AEROF.004.nc clean_12_14_aerooff/.
#cp clean_09_11_aerooff/CLEA4.1.AEROF.004.des clean_12_14_aerooff/.

#cd clean_12_14_aerooff/
#pwd
#./run_mesonh_xyz


#cd ../
#cp clean_12_14_aerooff/CLEA5.1.AEROF.004.nc clean_15_17_aerooff/.
#cp clean_12_14_aerooff/CLEA5.1.AEROF.004.des clean_15_17_aerooff/.

#cd clean_15_17_aerooff/
#pwd
#./run_mesonh_xyz

#cd ../
#cp clean_15_17_aerooff/CLEA6.1.AEROF.004.nc clean_18_20_aerooff/.
#cp clean_15_17_aerooff/CLEA6.1.AEROF.004.des clean_18_20_aerooff/.

#cd clean_18_20_aerooff/
#pwd
#./run_mesonh_xyz

#cd ../
cp clean_18_20_aerooff/CLEA7.1.AEROF.004.nc clean_21_23_aerooff/.
cp clean_18_20_aerooff/CLEA7.1.AEROF.004.des clean_21_23_aerooff/.


cd clean_21_23_aerooff/
pwd
./run_mesonh_xyz
