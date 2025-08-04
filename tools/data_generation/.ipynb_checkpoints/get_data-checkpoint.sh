#!/bin/bash
#SBATCH --job-name=get_data
#SBATCH --output=get_data_%j.out
#SBATCH --error=get_data_%j.err
#SBATCH --account="project_2010692"              # Name of the project
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00

# array of months (as integers in range from 1 to 12)
list_of_files=(./CAM5-ATRAS_AP3-CTRL/renamed/aerocom3_CAM5-ATRAS_AP3-CTRL_clt_ModelLevel_2010_monthly.nc,
./INCA_AP3-CTRL/renamed/aerocom3_INCA_AP3-CTRL_clt_ModelLevel_2010_monthly.nc,
./ECHAM6.3-SALSA2.0-met2010_AP3-CTRL/renamed/aerocom3_ECHAM6.3-SALSA2.0-met2010_AP3-CTRL_clt_ModelLevel_2010_monthly.nc,
./ECHAM6.3-HAM2.3-met2010_AP3-CTRL/renamed/aerocom3_ECHAM6.3-HAM2.3-met2010_AP3-CTRL_clt_ModelLevel_2010_monthly.nc,
./OsloCTM3v1.01-met2010_AP3-CTRL/renamed/aerocom3_OsloCTM3v1.01-met2010_AP3-CTRL_clt_ModelLevel_2010_monthly.nc,
)

sftp -i /users/attelaak/.ssh/id_rsa atte.laakso@uef.fi@aerocom-users.met.no
for M in "${list_of_files[@]}"; do
    sftp> get $M /scratch/project_2010692/attelaak/
    sftp> get -r $M /scratch/project_2010692/attelaak/
done
sftp> exit
