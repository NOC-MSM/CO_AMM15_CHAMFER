#!/bin/bash
#SBATCH --job-name=WAD
#SBATCH --time=00:20:00
#SBATCH --account=n01-chamfer
#SBATCH --partition=serial
#SBATCH --qos=serial
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1

module swap craype-network-ofi craype-network-ucx
module swap cray-mpich cray-mpich-ucx
module load cray-hdf5-parallel/1.12.2.1
module load cray-netcdf-hdf5parallel/4.9.0.1
module load libfabric

export OMP_NUM_THREADS=1

./nemo.exe
