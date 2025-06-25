# PEAR-TREE

![PEAR-TREE Logo](peartree_logo_.png)

paired ends of aberrant retrotransposons in phylogenetic trees

This repository contains the complete pipeline to find insertional mutations in whole genome sequencing data.

If you are searching for the pipeline used for targeted sequencing go to the repository PARTRIDGE

## How do I run it?

PEAR-TREE is run by calling the main.py file. See below how to run the different modes of this program.

The readme for the supplementary tools can be found [here](tools.md).

## Requirements

* python > 3.10, I recommend python3.13
* Installed packages specified in requirements.txt. See setup info
* bowtie2
* samtools
* 2bit file containing the assembly the bam files were mapped against
* indexed bowtie2 genome, ideally of a T2T reference of the organism you are studying. This assembly does not have to be the same as used to map the bamfiles against.

## Installation

* Download this repository: git clone [repo name]
* Install python3.13 (any version > 3.10 should work)
  * Get it from the [The Python homepage](https://www.python.org/downloads/)
  * Check that you have the correct version by using python -V
* Create a new venv
  * `python -m venv venv`
* Install the necessary packages
  * `source venv/bin/activate`
  * `pip install -r requirements.txt`
* Install bowtie2, I use version 2.5.4
  * [Github repo](https://github.com/BenLangmead/bowtie2)
  * use conda, homebrew or compile yourself
* Install samtools, I use version 1.21
  * [HTSLib and Samtools Homepage](https://www.htslib.org/download/)
  * use conda, homebrew or follow the instructions on the above webpage
* Download a 2bit file for the assembly your BAM-files were mapped against. See table below for sources
* Build a bowtie index for the latest assembly of the organism you are using. It is important to choose a T2T-assembly if available.

## Sources of 2bit files

Species | Assembly |  URL
--- |----------| ---
Homo sapiens | hg19     | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit  
Homo sapiens | hg38     | https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit 
Homo sapiens | hs1      | https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.2bit
Mus musculus | mm39     | https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.2bit 
Mus musculus | mm10 | https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.2bit

## Running PEAR-TREE

PEAR-TREE is run in four separate steps
1. Discovery: Run PEAR-TREE on all bam files you have, one single-threaded process is used per file. For this step no information on the genome is required and also no external files such as samtools,  bowtie2, 2bit files or indexes are required.


2. Combine Insertions: Combine the insertions detected during discovery file. For this process all requirements listed above are required.


3. Genotype: Run PEAR-TREE again on all bam files you have, one multi-threaded process is used per bam file. This sensitively re-genotypes all files based upon the consolidated list of insertions generated during step 2.


4. Combine Genotypes: This script consolidates all genotypes into one (potentially giant) csv file that can then be used for downstream analysis. The csv file contains one row per insertion and one column per bam file, specifying the genotype of that file for a given insertions. This step also performs extensive filtering defined in config.py

### Step 1: Discovery

```bash
python main.py --step discover --bam [bam] --out [outfile.txt.gz]
```

Parameter | Description
--- | ---
--step | discover for this step
--bam | Path to a single bam file
--out | Path to a single out-file, .txt.gz

This step runs on a single thread and requires at least 16 GB of ram.

*Recommended slurm settings*
```{slurm}
--mem=16G
--cpus-per-task=1
--t=11:59:59
--ntasks=1
```
### Step 2: Combine Insertions

```bash
python main.py --step combine_insertions --discovery_files [files] --out [outfile_stem] --threads 6
```

Parameter | Description
--- | ---
--step | combine_insertions for this step
--discovery_files | List of files generated in step 1, consider using bash expansions, e.g. `*.txt.gz`
--out | Stem of output files. This script actually generates a bunch of output files with this stem: {stem}.combined.txt.gz, {stem}.genotyping.txt.gz, {stem}.fa.gz, and {stem}.bam
--threads | Number of parallel threads, used for bowtie2

### Step 3: Genotyping

```bash
python main.py --step genotype --bam [bam] --out [outfile.txt.gz] --insertions [insertions.genotyping.txt.gz] --threads 6
```

 Parameter    | Description                                
--------------|--------------------------------------------
 --step       | genotype for this step                     
 --bam        | Path to a single bam file                  
 --out        | Path to the output file, is .txt.gz        
 --insertions | Path to the {stem}.genotyping.txt.gz file generated in Step 2
--threads | Number of parallel threads to use, I recommend 6, but any reasonable number (including 1) should work.


*Recommended slurm settings*
```{slurm}
--mem=16G
--cpus-per-task=16
--t=11:59:59
--ntasks=1
```
### Step 4: Combine Genotypes

```bash
python main.py --step combine_genotypes --genotypes [list of genotyping output files] --out [outfile.csv.gz] --threads 6\n',
```

Parameter | Description
--- | ---
--step | combine_genotype for this step
--genotypes | List of genotyping files generated during Step 3. Consider using bash syntax expansion, e.g. `*.txt.gz`
--out | Path to output file, .csv.gz
--threads | Number of parallel threads to use, I recommend 6-12, but any reasonable number (including 1) should work.


## Configuration

Configuration is defined in the config.py file. This table gives a detailed explanation of these settings

Parameter | Default | Description
--- | --- | ---
discovery.min_map | 60 | Minimal MAPQ value in the BAM file for a read to be even considered
discovery.min_clip_len | 12 | Minimum number of clipped bases required for a breakpoint to be discovered
discovery.min_evidence_reads_per_breakpoint | 2 | Minimum number of independent evidence reads required for each breakpoint
discovery.max_bp_window | 40 | maximum absolute number of bases delta between two breakpoint required for them to be combined to an insertion site.
discovery.max_homopolymer_len | 6 | maximum homopolymer length of mapped part of read immedietaly adjacent to the breakpoint
discovery.min_adapterlen_for_clip | 4 | minimum number of bases matching an adapter required for adapter clipping
discovery.min_good_bases | 10 | minimal number of bases of a clipped read required with a good quality (defined by the parameter below)
discovery.min_consensus_score_for_good_base | 2 | minimal score required for a good base, this corresponds to the consensus quality. 
discovery.min_breakpoints_aggregated_during_first_step | 2 | minimum number of reads with the same breakpoint aggregated needed in order to even start considering a breakpoint. At this stage both reads of a single fragment are counted as two independent reads, this is only corrected later.
discovery.max_read_count | 60 | maximum number of reads allowed to span a breakpoint region. This is designed to exclude regions of high-coverage, which are almost always artefact. This number has to be adjusted according to expected coverage and has to be set very high in case of enriched sequencing.
adapters | - | list of adapter sequences to clip. This has to be ajusted according to library prep and sequencing platform. The default uses the NebNext Adapters for Illumina
genotyping.max_bases | 12 | maximum number of bases used for genotyping.
genotyping.min_score_for_call | 10 | minimum score (+1 for match, -2 for mismatch) required to call an insertion
genotyping.reads_for_high_coverage | 60 | similar to discovery.max_read_count, but re-applied during genotyping file since not all breakpoints have been checked for all bam files during discovery phase.
combine_genotypes.min_wild-types | 20 | minimum number of bam files required to have a wild-type genotype in order for an insertion to be considered. This gets rid of insertions that don't vary between files but are not present in the reference genome.
combine_genotypes.min_insertions | 1 | minimum number of certain heterozygous or homozygous calls required for an insertion to pass filtering
combine_genotypes.max_artefact | 200 | maximum number of bam files allowed to have an "artefact" genotype (evidence of non-insertion non-reference clipped reads)
combine_genotypes.max_na | 18 | maximum number of bam files allowed to have an "NA" genotype (no coverage or high coverage).
combine_insertions.genome_2bit | - | Path to the 2bit file of the assembly used to map the bams
combine_insertions.exclude_files_with_many_insertions | 1000000 | exclude files with more than this number of insertions. Set to 1 Mio to basically switch off this feature, its probably better to do this manually, since insertional mutations only present in a single bam file will can not be detected if that bam file happens to have more than this number of insertions
combine_insertions.samtools_executable | - | Path to samtools
combine_insertions.bowtie2_executable | - | Path to bowtie2
combine_insertions.bowtie2_index | - | Path to bowtie2 index


Note: If using enriched and not whole-genome data, wild-types usually are not detected since the absence of an insertion leads to low-coverage. Thus, it might be advisable to set min_wild-types to 0 and max_na to 1000000


## License

PEAR-TREE - paired ends of aberrant retrotransposons in phylogenetic trees

Copyright (C) 2025 Jeremy Deuel <jeremy.deuel@usz.ch>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
