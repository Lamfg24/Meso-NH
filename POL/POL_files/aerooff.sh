#!/bin/bash
#MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence                                                                                     
#MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt                                                                                         
#MNH_LIC for details. version 1.              
#SBATCH -J MNH1                                                                                                                                                  
#SBATCH -N 25                                                                                                                                                    
#SBATCH -n 1800                                                                                                                                                     
#SBATCH -o Examples.eo%j   #                                                                                                                                       
#SBATCH -e Examples.eo%j   #                                                                                                                                       
#SBATCH -t 10:00:00    # time limit                                                                                                                               
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

#cd local_23_02_aerooff/
#pwd
#./run_mesonh_xyz


#cd ../
#cp local_23_02_aerooff/LOCA1.1.AEROF.004.nc local_03_05_aerooff/.
#cp local_23_02_aerooff/LOCA1.1.AEROF.004.des local_03_05_aerooff/.

#cd local_03_05_aerooff/
#pwd
#./run_mesonh_xyz


#cd ../
#cp local_03_05_aerooff/LOCA2.1.AEROF.004.nc local_06_08_aerooff/.
#cp local_03_05_aerooff/LOCA2.1.AEROF.004.des local_06_08_aerooff/.

#cd local_06_08_aerooff/
#pwd
#./run_mesonh_xyz


#cd ../
#cp local_06_08_aerooff/LOCA3.1.AEROF.004.nc local_09_11_aerooff/.
#cp local_06_08_aerooff/LOCA3.1.AEROF.004.des local_09_11_aerooff/.

#cd local_09_11_aerooff/
#pwd
#./run_mesonh_xyz

#cd ../
#cp local_09_11_aerooff/LOCA4.1.AEROF.004.nc local_12_14_aerooff/.
#cp local_09_11_aerooff/LOCA4.1.AEROF.004.des local_12_14_aerooff/.

#cd local_12_14_aerooff/
#pwd
#./run_mesonh_xyz

#cd ../
#cp local_12_14_aerooff/LOCA5.1.AEROF.004.nc local_15_17_aerooff/.
#cp local_12_14_aerooff/LOCA5.1.AEROF.004.des local_15_17_aerooff/.

#cd local_15_17_aerooff/
#pwd
#./run_mesonh_xyz


#cp local_15_17_aerooff/LOCA6.1.AEROF.004.nc local_18_20_aerooff/.
#cp local_15_17_aerooff/LOCA6.1.AEROF.004.des local_18_20_aerooff/.

#cd local_18_20_aerooff/
#pwd
#./run_mesonh_xyz


cd local_21_23_aerooff/
pwd
./run_mesonh_xyz

