for i in raw/*.sample.dupmarked.bam; do
	i=$(basename $i)
	i=${i%.sample.dupmarked.bam}
	if [ -f insertions/$i.txt.gz ]; then
		echo "$i already run, skipping"
	else
		echo "$i submitting run_extract_insertions.slurm"
		sbatch --job-name $i -o "slurm_log/$i-insertions-%A.out" run_extract_insertions.slurm $i
	fi
done
