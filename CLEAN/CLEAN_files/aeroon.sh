#!/bin/bash
#MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence                                                                                     
#MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt                                                                                         
#MNH_LIC for details. version 1.              
#SBATCH -J MNH1                                                                                                                                                  
#SBATCH -N 50                                                                                                                                                    
#SBATCH -n 3600                                                                                                                                                    
#SBATCH -o Examples.eo%j   #                                                                                                                                       
#SBATCH -e Examples.eo%j   #                                                                                                                                       
#SBATCH -t 15:10:00    # time limit                                                                                                                               
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

#cd clean_23_02_aeroon/
#pwd
#./run_mesonh_xyz


#cd ../
#cp clean_23_02_aeroon/CLEA1.1.AERON.004.nc clean_03_05_aeroon/.
#cp clean_23_02_aeroon/CLEA1.1.AERON.004.des clean_03_05_aeroon/.

#cd clean_03_05_aeroon/
#pwd
#./run_mesonh_xyz


#cd ../
#cp clean_03_05_aeroon/CLEA2.1.AERON.004.nc clean_06_08_aeroon/.
#cp clean_03_05_aeroon/CLEA2.1.AERON.004.des clean_06_08_aeroon/.

#cd clean_06_08_aeroon/
#pwd
#./run_mesonh_xyz


#cd ../
#cp clean_06_08_aeroon/CLEA3.1.AERON.004.nc clean_09_11_aeroon/.
#cp clean_06_08_aeroon/CLEA3.1.AERON.004.des clean_09_11_aeroon/.

#cd clean_09_11_aeroon/
#pwd
#./run_mesonh_xyz


#cd ../
#cp clean_09_11_aeroon/CLEA4.1.AERON.004.nc clean_12_14_aeroon/.
#cp clean_09_11_aeroon/CLEA4.1.AERON.004.des clean_12_14_aeroon/.

#cd clean_12_14_aeroon/
#pwd
#./run_mesonh_xyz


#cd ../
#cp clean_12_14_aeroon/CLEA5.1.AERON.004.nc clean_15_17_aeroon/.
#cp clean_12_14_aeroon/CLEA5.1.AERON.004.des clean_15_17_aeroon/.

#cd clean_15_17_aeroon/
#pwd
#./run_mesonh_xyz


#cd ../
cp clean_15_17_aeroon/CLEA6.1.AERON.004.nc clean_18_20_aeroon/.
cp clean_15_17_aeroon/CLEA6.1.AERON.004.des clean_18_20_aeroon/.

cd clean_18_20_aeroon/
pwd
./run_mesonh_xyz

cd ../
cp clean_18_20_aeroon/CLEA7.1.AERON.004.nc clean_21_23_aeroon/.
cp clean_18_20_aeroon/CLEA7.1.AERON.004.des clean_21_23_aeroon/.



cd clean_21_23_aeroon/
pwd
./run_mesonh_xyz
