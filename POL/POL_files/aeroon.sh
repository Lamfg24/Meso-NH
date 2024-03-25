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

#cd local_23_02_aeroon/
#pwd
#./run_mesonh_xyz

#cd ../
#cp local_23_02_aeroon/LOCA1.1.AERON.004.nc local_03_05_aeroon/.
#cp local_23_02_aeroon/LOCA1.1.AERON.004.des local_03_05_aeroon/.

#cd local_03_05_aeroon/
#pwd
#./run_mesonh_xyz


#cd ../
#cp local_03_05_aeroon/LOCA2.1.AERON.004.nc local_06_08_aeroon/.
#cp local_03_05_aeroon/LOCA2.1.AERON.004.des local_06_08_aeroon/.

#cd local_06_08_aeroon/
#pwd
#./run_mesonh_xyz


#cd ../
#cp local_06_08_aeroon/LOCA3.1.AERON.004.nc local_09_11_aeroon/.
#cp local_06_08_aeroon/LOCA3.1.AERON.004.des local_09_11_aeroon/.

#cd local_09_11_aeroon/
#pwd
#./run_mesonh_xyz

#cd ../
#cp local_09_11_aeroon/LOCA4.1.AERON.004.nc local_12_14_aeroon/.
#cp local_09_11_aeroon/LOCA4.1.AERON.004.des local_12_14_aeroon/.

#cd local_12_14_aeroon/
#pwd
#./run_mesonh_xyz


#cd ../
#cp local_12_14_aeroon/LOCA5.1.AERON.004.nc local_15_17_aeroon/.
#cp local_12_14_aeroon/LOCA5.1.AERON.004.des local_15_17_aeroon/.

#cd local_15_17_aeroon/
#pwd
#./run_mesonh_xyz


#cd ../
#cp local_15_17_aeroon/LOCA6.1.AERON.004.nc local_18_20_aeroon/.
#cp local_15_17_aeroon/LOCA6.1.AERON.004.des local_18_20_aeroon/.

#cd local_18_20_aeroon/
#pwd
#./run_mesonh_xyz

cd local_21_23_aeroon/
pwd
./run_mesonh_xyz










