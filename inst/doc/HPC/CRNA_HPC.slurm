#!/bin/bash
#SBATCH --job-name="CRNA_HPC"
#SBATCH --time=24:00:00
#SBATCH -N 5
#SBATCH -n 100
#SBATCH --exclusive
#SBATCH --mem=120gb
#SBATCH --output=/home/jxd101/RESULTS/METHODS/CRNA/REAL/CRNA_HPC.log

  # Load modules
module load gcc/6.3.0
module load R-ext-4.0.2/4.0.2

echo "My job will use the following nodes:"
echo "-------------------------------"
echo $SLURM_JOB_NODELIST
echo "-------------------------------"
echo "My job will use the following number of tasks:"
echo "-------------------------------"
echo $SLURM_NTASKS
echo "-------------------------------"

  # Copy all files needed for the job from the submit directory ($SLURM_SUBMIT_DIR),
  # i.e. the directory where the job is submitted, 
  # to the temporary scratch directory ($PFSPDIR) associated with the job.
  # (If needed input files, workspaces, etc... are in another location, then copy from that location instead)
cp -r * $PFSDIR

  # cd to the directory where the job was submitted
cd $PFSDIR

  # Run R
Rscript /home/jxd101/CODES/R/METHODS/CRNA/REAL/CRNA_HPC.r "SOCKET"

  # Copy results from tmp to submit directory
cp -ru * $SLURM_SUBMIT_DIR

  # Delete my scratch job files
rm -Rf `ls -la "$PFSDIR"/* | grep 'jxd101' | awk ' { print $9 } '`
