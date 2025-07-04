# PEAR-TREE - paired ends of aberrant retrotransposons in phylogenetic trees
#
# Copyright (C) 2025 Jeremy Deuel <jeremy.deuel@usz.ch>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.


import sys, time, os, argparse

from config import CONFIG
def print_head(phase):
    print('\033[36;1;4mPEAR-TREE\033[0m')
    print('\033[36mpaired ends of aberrant retrotransposons in phylogenetic trees\033[0m')
    print(f"\033[90mConfiguration: {CONFIG}\033[0m")
    print(f"Current Phase: {phase}")

if __name__ == '__main__':
    clm = argparse.ArgumentParser(
        prog = 'PEAR-TREE',
        description = 'paired ends of aberrant retrotransposons in phylogenetic trees',
        usage = 'python main.py --step discover --bam [bam] --out [outfile.txt.gz]\n' +
                'python main.py --step combine_insertions --discovery_files [files] --out [outfile_stem] --threads [1]\n' +
                'python main.py --step genotype --bam [bam] --out [outfile.txt.gz] --insertions [insertions.genotyping.txt.gz] --threads [1]\n' +
                'python main.py --step combine_genotypes --genotypes [list of genotyping output files] --out [outfile.csv.gz] --threads [1]\n',

        epilog = f'Version {CONFIG["version"]}, Created by Jeremy Deuel <jeremy.deuel@usz.ch>'
    )
    clm.add_argument('--step', '-s', help="Step to execute, can be discover, genotype, combine_insertions or combine_genotypes", nargs=1)
    clm.add_argument('--bam', '-d', help='Path to input BAM file.', nargs=1)
    clm.add_argument('--out','-o', help='Path to output file, gzipped', nargs=1)
    clm.add_argument('--insertions', '-i', help='Path to summarised insertions file, gzipped', nargs=1)
    clm.add_argument('--discovery_files', '-f', help="Path to insertion files generated in the discovery step, gzipped", nargs="*")
    clm.add_argument('--genotypes', '-g', help='Path to genotype files, gzipped', nargs="*")
    clm.add_argument('--threads', '-@', help='Number of threads to use, defaults to one, ignored during discovery', nargs=1, default=[1])

    args = clm.parse_args(sys.argv[1:])
    if not args:
        clm.print_help()
        exit(1)
    if not args.step:
        clm.print_help()
        exit(1)
    if args.step[0] == 'discover':

        if not args.bam or not args.out:
            clm.print_help()
            exit(1)
        from discovery import Discovery
        print_head('Discovery')
        input_bam = args.bam[0]
        output_path = args.out[0]
        if not os.path.exists(input_bam):
            raise FileNotFoundError(f"input bam file {input_bam} does not exist!")
        print(f"input bam: {input_bam}, output file: {output_path}")
        discovery = Discovery()
        discovery.discovery(input_bam)
        discovery.output(output_path)
        exit(0)

    if args.step[0] == 'genotype':

        if not args.bam or not args.out or not args.insertions:
            clm.print_help()
            exit(1)
        input_bam = args.bam[0]
        output_path = args.out[0]
        input_insertions = args.insertions[0]
        threads = int(args.threads[0])
        if threads > os.cpu_count():
            raise ValueError(f"this machine only has {os.cpu_count()} CPUs, do not run this script with more threads, you have requested {threads}.")
        if not os.path.exists(input_bam):
            raise FileNotFoundError(f"input bam file {input_bam} does not exist!")
        if not os.path.exists(input_insertions):
            raise FileNotFoundError(f"input insertions file {input_insertions} does not exist!")
        print_head('Genotyping')
        print(f"input bam: {input_bam}, input insertions: {input_insertions} output file: {output_path}, threads: {threads}")
        from genotype import genotype
        genotype(input_insertions, input_bam, output_path, threads)
        exit(0)

    if args.step[0] == 'combine_genotypes':
        if not args.genotypes or not args.out:
            clm.print_help()
            exit(1)
        input_files = args.genotypes
        output_path = args.out[0]
        threads = int(args.threads[0])
        if threads > os.cpu_count():
            raise ValueError(f"this machine only has {os.cpu_count()} CPUs, do not run this script with more threads, you have requested {threads}.")
        for f in input_files:
            if not os.path.exists(f):
                raise FileNotFoundError(f"input genotyping file {f} does not exist!")
        print_head('Combining Genotypes')
        print(f"output file: {output_path}, threads: {threads}, input genotyping files: {input_files}")
        from combine_genotypes import collect_genotype
        collect_genotype(input_files, output_path, threads)
        exit(0)

    if args.step[0] == 'combine_insertions':
        if not args.discovery_files or not args.out:
            clm.print_help()
            exit(1)
        discovery_files = args.discovery_files
        output_stem = args.out[0]
        if output_stem[-3:] == '.gz':
            output_stem = output_stem[:-3]
        threads = int(args.threads[0])
        if threads > os.cpu_count():
            raise ValueError(
                f"this machine only has {os.cpu_count()} CPUs, do not run this script with more threads, you have requested {threads}.")
        for f in discovery_files:
            if not os.path.exists(f):
                raise FileNotFoundError(f"input insertions file {f} does not exist!")
        if not os.path.exists(CONFIG['combine_insertions']['samtools_executable']):
            raise FileNotFoundError(f"samtools not found in {CONFIG['combine_insertions']['samtools_executable']}, specify this in the config file")
        if not os.path.exists(CONFIG['combine_insertions']['bowtie2_executable']):
            raise FileNotFoundError(
                f"bowtie2 not found in {CONFIG['combine_insertions']['bowtie2_executable']}, specify this in the config file")
        if not os.path.exists(f"{CONFIG['combine_insertions']['bowtie2_index']}.1.bt2"):
            raise FileNotFoundError(
                f"bowtie2 index not found in {CONFIG['combine_insertions']['bowtie2_index']}, generate this using bowtie2 index `{CONFIG['combine_insertions']['bowtie2_executable']}-build [REFERENCE_GENOME.fasta] {CONFIG['combine_insertions']['bowtie2_index']}`")
        print_head('Combining Insertions')
        print(f"output stem: {output_stem}, threads: {threads}, input genotyping files: {discovery_files}")
        print("other output files generated:")
        print(f"* combined insertions file: {output_stem}.combined.txt.gz")
        print(f"* combined insertions files for genotyping: {output_stem}.genotyping.txt.gz")
        print(f"* fasta file used to run bowtie2: {output_stem}.fa.gz")
        print(f"* bam file generated by bowtie2: {output_stem}.bam")
        from combine_insertions import combine_insertions
        combine_insertions(discovery_files,
                           f"{output_stem}.genotyping.txt.gz",
                           f"{output_stem}.combined.txt.gz",
                           f"{output_stem}.fq.gz",
                           f"{output_stem}.bam",
                           threads)
        exit(0)

    clm.print_help()
    exit(1)