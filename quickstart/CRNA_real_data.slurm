#!/bin/bash
#SBATCH --job-name="CRNA_REAL"
#SBATCH --time=138:00:00
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -p batch
#SBATCH --exclusive
#SBATCH --mem=150gb
#SBATCH --output=/home/fxe35/RESULTS/METHODS/CRNA/SIMS/CRNA_REAL%j.log

module load R
module load GMP/6.2.1-GCCcore-11.2.0

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
Rscript ./CRNA_real.r "SOCKET"

  # Copy results from tmp to submit directory
cp -ru * $SLURM_SUBMIT_DIR

  # Delete my scratch job files
rm -Rf `ls -la "$PFSDIR"/* | grep 'fxe35' | awk ' { print $9 } '`