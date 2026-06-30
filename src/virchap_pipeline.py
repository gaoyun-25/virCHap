import os
import sys
import argparse
import subprocess
try:
    from src import prepare as prepare
    from src import clustering as clustering
except ImportError:
    import prepare as prepare
    import clustering as clustering


def log_error(message):
    print_stderr(message)

# Print to STDERR
def print_stderr(message, error=True):
    if error:
        print(f"ERROR {message}", file=sys.stderr)
    else:
        print(f"{message}", file=sys.stderr)

# Utility to check if a tool is in PATH
def check_tool_in_path(tool_name):
    result = subprocess.run(f"which {tool_name}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.returncode == 0

# Check if required tools are available
def check_required_tools():
    tools = ['virchap', 'lofreq', 'minimap2', 'samtools', 'gatk']
    missing_tools = [tool for tool in tools if not check_tool_in_path(tool)]

    if missing_tools:
        log_error(f"The following tools are missing or not in PATH: {', '.join(missing_tools)}")
        sys.exit(1)

# Run a command and check if it succeeded
def run_command(command, error_message):
    p = subprocess.run(command, shell=True)
    if p.returncode != 0:
        log_error(error_message)
        sys.exit(1)

def pipeline(args,ref_len):
    # Create folder
    all_paths=[]
    out_dir = os.path.join(args.output_dir, args.out_name)
    phased_path = os.path.join(out_dir, "Phased")
    variants_long_read_path = os.path.join(out_dir, "Variants", "Reads")
    variants_path = os.path.join(out_dir, "Variants", "VCF")
    Bam_path = os.path.join(out_dir, "Mapped")

    all_paths += [out_dir, phased_path, variants_long_read_path, variants_path,Bam_path]

    for path in all_paths:
        os.makedirs(path, exist_ok=True)

    bam_file = os.path.join(Bam_path, "minimap2.sorted.bam")
    vcf_file = os.path.join(variants_path, "lofreq.vcf")

    # BAM file
    if not os.path.exists(bam_file):
        run_command(
            f"minimap2 -a {args.ref_file} {args.reads_file} -t {args.threads} | samtools sort -@ {args.threads} -o {bam_file}",
            "Alignment or BAM sorting failed"
        )

    if not os.path.exists(bam_file + ".bai"):
        run_command(f"samtools index {bam_file}", "BAM file indexing failed")

    # VCF file
    if not os.path.exists(vcf_file): 
        run_command(f"lofreq call -B -f {args.ref_file} {bam_file} -o {vcf_file} --force-overwrite", "LoFreq variant calling failed")

    # Extract heterozygous positions (vcf sam 都是1-index)
    variants_positions_file = os.path.join(variants_path, args.out_name+".hetSNPs.positions.vcf.tsv")
    prepare.get_variants_positions(vcf_file, variants_positions_file)

    # #Reduce long reads to their heterozygous SNPs
    reduced_long_reads_file = os.path.join(variants_long_read_path, args.out_name+".hetPositions.SNPxLongReads.vcf.chr.tsv")
    prepare.reduce_long_reads_to_SNPS(bam_file, variants_positions_file, reduced_long_reads_file, args.use_supplementary, args.min_mapq)

    split_number = 5
    clustering.phasing(reduced_long_reads_file, bam_file,args.threads, out_dir, args.out_name, args.max_ID, args.min_len, args.min_ovl, args.min_sim, args.min_overlap_len, args.min_ovl_len, args.min_sim_first, args.min_abundance,args.min_overlap_len_first,args.min_ovl_len_first,split_number,ref_len)

    return None


def main():
    parser = argparse.ArgumentParser(prog='virchap_pipeline',
        usage="virchap_pipeline --reference reference.fasta --reads reads.fq --output outdir [options]",
        description="Viral haplotype reconstruction tool using reads and reference as input.",
    )
    # required
    parser.add_argument('--reference', '-r', required=True, dest='ref_file', help='Reference file (fasta).')
    parser.add_argument('-i', '--reads', required=True, dest='reads_file', help='Input reads')
    parser.add_argument('--output', '-o', required=True, dest='output_dir', help='Output folder path.')

    # optional
    parser.add_argument('--threads', '-t', dest='threads', nargs="?", default=4, type=int, help='Maximum number of CPU threads (default: %(default)s).')
    parser.add_argument('--sample', dest='out_name', nargs="?", default='phased', type=str, help='Name of sample to phase and the output name (default: %(default)s).')
    parser.add_argument('--use-supplementary', dest='use_supplementary', nargs="?", default=False,help='Use also supplementary alignments (default: false, ignore supplementary alignments).')
    parser.add_argument('--min_mapq', dest='min_mapq', nargs="?", default=20, type=int, help='Minimum mapping quality (default: %(default)s)')
    parser.add_argument('--minOverlapLenFirst', dest='min_overlap_len_first', nargs="?", default=25, type=int, help='Minimum length of overlap between two reads to build an overlap graph (default: %(default)s).')
    parser.add_argument('--minOvlLenFirst', dest='min_ovl_len_first', nargs="?", default=10, type=int, help='Minimum number of common heterozygous SNPs between two reads that have fewer than MIN_OVERLAP heterozygous SNPs in common (default: %(default)s).')
    parser.add_argument('--min_sim_first', dest='min_sim_first', nargs="?", default=0.95, type=float, help='Minimum similarity between two reads to build an overlap graph in first clustering for short genomes, such as virus and baterial (default: %(default)s)')
    parser.add_argument('--minOvl', dest='min_ovl', nargs="?", default=0.1, type=float, help='Minimum percentage of overlap between two reads that have fewer than MIN_OVERLAP heterozygous SNPs in common (default: %(default)s).')
    parser.add_argument('--minOverlapLen', dest='min_overlap_len', nargs="?", default=25, type=int, help='Minimum length of overlap between two cluster consensuses to merge clsuers (default: %(default)s).')
    parser.add_argument('--minOvlLen', dest='min_ovl_len', nargs="?", default=5, type=int, help='Minimum number of common heterozygous SNPs between two cluster consensuses that have fewer than MIN_OVERLAP heterozygous SNPs in common (default: %(default)s).')
    parser.add_argument('--minSim', dest='min_sim', nargs="?", default=0.3, type=float, help='Minimum similarity between two cluster consensuses to merge clusters (default: %(default)s).')
    parser.add_argument('--maxID', dest='max_ID', nargs="?", default=0.01, type=float, help='Maximum identity parameter, determines how different two clusters must be to prevent them from merging (default: %(default)s).')
    parser.add_argument('--minAbundance', dest='min_abundance', nargs="?", default=0.001, type=float, help='Minimum abundance parameter, filters haplotypes with abundance less than MINIMUM ABUNDANCE (default: %(default)s).')
    parser.add_argument('--minLen', dest='min_len', nargs="?", default=0, type=int, help='Minimum number of reads each cluster must have (default: %(default)s).')
    parser.add_argument('--platform', '-p', dest='platform', nargs="?", default="ont", type=str, help='sequencing platform: pb or ont (default: %(default)s).')

    args = parser.parse_args()

    # 若用户未指定--min_sim_first参数，则根据sequencing platform指定--min_sim_first的值
    minSimFirst_specified = any(arg.startswith('--min_sim_first') for arg in sys.argv)
    if not minSimFirst_specified:
        if args.platform == "ont":
            args.min_sim_first = 0.95
        elif args.platform == "pb":
            args.min_sim_first = 0.98

    maxID_specified = any(arg.startswith('--maxID') for arg in sys.argv)
    if not os.path.isfile(args.ref_file+".fai"):
        p=subprocess.run(["samtools","faidx",args.ref_file],stderr=subprocess.PIPE,stdout=subprocess.PIPE, universal_newlines=True)

    ref_len = dict()
    with open(args.ref_file+".fai") as f:
        for line in f:
            line = line.strip().split("\t")
            ref_len[line[0]] = int(line[1])
    if not maxID_specified:
        if int(sum(ref_len.values())/len(ref_len)) > 120000:
            args.max_ID = 0.03
    # run nTChap
    pipeline(args,ref_len)

    print("pipeline finished")

if __name__ == "__main__":
    main()
    exit(0)

