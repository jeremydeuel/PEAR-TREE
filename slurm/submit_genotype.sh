for i in raw/M*.sample.dupmarked.bam; do
	i=$(basename $i)
	i=${i%.sample.dupmarked.bam}
	if [ -f genotypes/$i.txt.gz ]; then
		echo "$i already run, skipping"
	else
		echo "$i submitting run_genotype.slurm"
		sbatch --job-name $i -o "slurm_log/$i-genotype-%A.out" run_genotype.slurm $i
	fi
done
