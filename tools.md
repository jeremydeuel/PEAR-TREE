# Tools for geminate

This readme document describes the additional tools available for geminate. They are located in the tools folder.

## Annotate

```{bash}
python tools/annotate.py
```

This tool annotates a excel output of insertions.

### Input and Configuration

This tool is configured directly in the source code by setting these constants:
* **EXCEL_IN** Excel file with at least the column "insertion" containing the names of the insertions that sould be annotated
* **EXCEL_OUT** Output path for the annotated Excel list
* **INSERTIONS_FILE** path to insertions.combined.txt.gz generated during step 2 of geminate. If necessary multiple files can be concatenated by using the `cat` command.
* **TMP_FASTA** temporary fasta file that will be used for the _hmmer_ and _bowtie2_ steps
* **DFAM_SCRIPT*** path to dfamscan.pl, downlaod from  https://dfam.org/help/tools
* **DFAM_HMM** path to the species-specific hmm file, generated with the filter_hmm.py file
* **TMP_DFAM** temporary dfam output file
* **TMP_SAM** temporary bowtie2 output file
* **BOWTIE_SCRIPT** path to bowtie2 executable
* **BOWTIE_REF** path to bowtie2 index
* **CHAINFILE** path to chainfile linking assembly of bowtie2 index to assembly used in the analysed bam files

### Output description

Remember: The insertion point is annotated as
{seqname}:{right_bp}-{left_bp}, where
* seqname: Chromosome / Contig name
* left_bp: 0-based position of the first base after the left-clipped sequence
* right_bp: 0-based position of the first base that is no longer mapped (first clipped base)
* In case of a target site duplication (TSD) left_bp is smaller than right_bp; in case of a deletion vice versa.

Column  |   Description
--- | ---
left_seq | left-clipped sequence (5' end on + strand), relative to the + strand of the insertion site
right_seq | right-clipped sequence (3' end on the + strand), relative to the + strand of the insertion site
left_dfam | best matching DFAM model for left-clipped part
right_dfam | best matching DFAM model for right-clipped part
left_align | unique (mapq>=40) end-to-end alignment for the left-clipped part, in the assembly of the bowtie2 reference
right_align | unique (mapq>=40) end-to-end alignment for the right-clipped part, in the assembly of the bowtie2 reference
hg38_ins_bp_left | last base of the uniquely mapping left-clipped sequence, in the bam file's assembly (lift-over with the provided chain-file)
hg38_ins_bp_right | first base of the uniquely mapping right-clipped sequence, in the bam file's assembly

## Check Germline Coverage

```{bash}
python tools/check_germline_coverage.py
```
Using an excel list of insertions, this tool generates an excel list of the number of bam files that insertion was identified in during step 1. This can be used to estimate sensitivity by inputing known germ-line inseritons, ideally heterozygous ones.

### Input and Configuration

This tool is configured directly in the source code by setting these constants:

* `d = pd.read_excel("all.insertions.xlsx")`: Path to the input file, should contain a column named "insertion"
* `files = list(os.listdir('PD45517'))+list(os.listdir('PD45534'))+list(os.listdir('PD48402'))`: List all the files obtained during the first step of geminate. This is only used to obtain the original bam file count.
* `with gzip.open(f'{patient}.insertions.combined.txt.gz', 'rt') as fh:` path to the insertions.combined.txt.gz files obtained during the second step of geminate
* `pd.DataFrame(output).to_excel('first.step.stats.xlsx')`: path to the desired output file


### Output description

Column | Description
--- | ---
Colonies | total number of bam files for a specific individual
n_first_step | number of bam files the inseriton has been identified in during the first step.

## filter HMM

This script filters a (huge) DFAM file with all possible HMMs for a given species. This script is configured in the command line. The output is used for the annotate tool

```{bash}
python tools/filter_hmm.py --infile Dfam.hmm --outfile Dfam.mm.hmm --species "Mus musculus"
python tools/filter_hmm.py --infile Dfam.hmm --outfile Dfam.hs.hmm --species "Homo sapiens"

```

### Input

Download the full DFAM HMM library https://dfam.org/releases/current/families/Dfam.hmm.gz. Make sure to spell the species exactly as in the DFAM HMM File, e.g. "Mus musculus" or "Homo sapiens" (Check for TaxName entries if necessary).

