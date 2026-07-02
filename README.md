# virChap
virChap is a reference-based viral haplotype reconstruction algorithm for long reads.

## Installation

### Using Conda
```
conda install -c conda-forge -c bioconda virchap
```

### Install from YAML

Two environment YAML files are provided for reproducible installation:

- `environment.yml`: minimal recipe with core dependencies (Python 3.8, samtools, gatk4, bwa, minimap2, LoFreq and Python packages).
- `virchap_environment.yml`: fully resolved environment with exact package versions and build strings.

Create and activate the environment using `environment.yml`:
```
conda env create -f environment.yml
conda activate virchap
cd virCHap
pip install . 
```

Or reproduce the exact resolved environment:
```
conda env create -f virchap_environment.yml
conda activate virchap
```

### Installation with pip
```
conda create -n virchap python==3.8
conda activate virchap
conda install -c bioconda lofreq minimap2 samtools gatk4
cd virCHap
pip install .
```

### Run without installation
If you do not wish to install virChap, you can run it directly from the source directory after satisfying the dependencies manually.
#### Dependencies:
```
python==3.8; lofreq; minimap2; samtools; gatk4; karateclub; sortedcontainers; pysam
conda install -c bioconda lofreq minimap2 samtools gatk4 karateclub
pip install sortedcontainers pysam
```
virChap provides two main entry points depending on your input data:
#### 1. Using raw reads and a reference (pipeline mode)
#### Command:
```
python virchap_pipeline.py --reference reference.fasta --reads reads.fq --output outdir
```
#### 2. Using pre‑aligned BAM, VCF and reference (phasing mode)

#### Command:
```
python phase_pipeline.py --bam BAM/SAM --vcf VCF --reference reference.fasta --output outdir
```
We recommend using high-quality, pre-filtered VCF inputs, as virCHap does not perform filtering.
## Detailed usage 
```
Running virchap_pipeline：
usage: virchap_pipeline --reference reference.fasta --reads reads.fq/reads.fq.gz --output outdir [options]

Running virchap:
usage: virchap --bam BAM/SAM --vcf VCF --reference reference.fasta --output outdir [options]


virchap_pipeline required arguments:

  --reference REF_FILE, -r REF_FILE         Reference file (fasta).
  -i READS_FILE, --reads READS_FILE         Input reads
  --output OUTPUT_DIR, -o OUTPUT_DIR        Output folder path.


virchap required arguments:
  --bam BAM_FILE                             BAM or SAM with long reads.
  --vcf VCF_FILE                             VCF file with variants to be phased.
  --reference REF_FILE, -r REF_FILE          Reference file.
  --output OUTPUT_DIR, -o OUTPUT_DIR         Output folder path.


optional arguments:
  -h, --help            show this help message and exit

  --threads [THREADS], -t [THREADS]          Maximum number of CPU threads (default: 4).
  --sample OUT_NAME                          Name of sample to phase and the output name (default: phased).
  --minOverlapLenFirst [MIN_OVERLAP_LEN_FIRST]  Minimum length of overlap between two reads to build an overlap graph (default: 25).
  --minOvlLenFirst [MIN_OVL_LEN_FIRST]       Minimum number of common heterozygous SNPs between two reads that have fewer than MIN_OVL_LEN_FIRST heterozygous SNPs in common (default: 10).  
  --min_sim_first [MIN_SIM_FIRST]            Minimum similarity between two reads to build an overlap graph in first clustering (default: 0.95).
  --minOvl [MIN_OVL]                         Minimum percentage of overlap between two reads or two cluster consensuses that have fewer than MIN_OVL heterozygous SNPs in common (default: 0.1).
  --minOverlapLen [MIN_OVERLAP_LEN]          Minimum length of overlap between two cluster consensuses to merge clsuers (default: 25).
  --minOvlLen [MIN_OVL_LEN]                  Minimum number of common heterozygous SNPs between two cluster consensuses that have fewer than MIN_OVERLAP heterozygous SNPs in common (default: 5).
  --maxID [MAX_ID]                           Maximum identity parameter, determines how different two clusters must be to prevent them from merging (default: 0.01).
  --minSim [MIN_SIM]                         Minimum similarity between two cluster consensuses to merge clusters (default: 0.3).
  --minLen [MIN_LEN]                         Minimum number of reads each clusters must have (default: 0).
  --minAbundance [MIN_ABUNDANCE]             Minimum abundance parameter, filters haplotypes with abundance less than MIN_ABUNDANCE (default: 0.001).
  --use-supplementary [USE_SUPPLEMENTARY]    Use supplementary alignments (default: ignore supplementary alignments).
  --min_mapq [MIN_MAPQ]                      Minimum mapping quality (default: 20)
```

## Output Results
Two tab-separated value (.tsv) files and one fasta files are in OUTPUT_DIR/Phased:  
+ "*_variants.tsv": predicted name (cluster name_number of reads_abundance,Contig_name: chromosome name), chromosome, position, base.  
+ "*_clustered_read_name.tsv": predicted name, read name.  
+ "*_consensus.fasta": predicted name,  the base-level sequences.  

## Example
Pipeline mode (reads + reference)
```
virchap_pipeline --reference example/OR483991.1.fasta --reads example/reads.fastq.gz --output example/
```
Phasing mode (BAM + VCF + reference)
```
virchap --bam example/reads.sorted.bam --vcf example/reads.SNP.vcf --reference example/OR483991.1.fasta --output example/
```

## Any questions
With any questions. Please, contact: gaoyun3304@163.com
