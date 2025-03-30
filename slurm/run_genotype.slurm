#!/bin/bash
#!
#! Template SLURM Batch Script
#! <jd862@cam.ac.uk> 2021-Aug-09
#! Index bamfiles in directory specified by 1


#! ################### SLURM SETTINGS

#! Name of the job:
#SBATCH -A GREEN-SL3-CPU
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:59:59
#SBATCH --mail-type=FAIL
#SBATCH --partition=icelake
PYTHONUNBUFFERED=1 #this makes the output immediate
echo "This is run_sirt.slurm using parameter $1"
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
if [[ "$SLURM_JOB_PARTITION" == "cclake" ]]; then
        module load rhel7/default-ccl
fi
if  [[ "$SLURM_JOB_PARTITION" == "icelake" ]]; then
        module load rhel8/default-icl              # REQUIRED - loads the basic$
fi
if [[ "$SLURM_JOB_PARTITION" == "sapphire" ]]; then
	module load rhel8/default-sar
fi

#module load python/3.11.9/gcc/nptrdpll
source venv/bin/activate

if [ ! -f raw/$1.sample.dupmarked.bam.bai ]; then
	echo "Creating index for $1"
	samtools index raw/$1.sample.dupmarked.bam
fi
python geminate/src/main.py --step genotype  --bam raw/$1.sample.dupmarked.bam  --out genotypes/$1.txt.gz --threads 1
#python sirt/genotype.py insertions.genotyping.txt.gz raw/$1.sample.dupmarked.bam genotypes/$1.txt.gz
