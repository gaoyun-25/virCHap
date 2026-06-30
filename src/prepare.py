import multiprocessing
import re
from collections import defaultdict
import pysam
import os

class AlignedRead():

    def __init__(self, read_name, variant_seq, mapq, is_supplementary, is_reversed, ref_start, ref_end, ref_name):
        self.read_name = read_name
        self.variant_seq = variant_seq
        self.mapq = mapq
        self.is_supplementary = is_supplementary
        self.is_reversed = is_reversed
        self.ref_start = ref_start  # 1-index
        self.ref_end = ref_end
        self.ref_name = ref_name

    def distance(self, other):
        return max(
            other.ref_end - self.ref_start,
            other.ref_start - self.ref_end,
            0,
        )


def get_variants_positions(variants_het_snp_vcf_file, variants_positions_file):
    variants_het_snp = dict()
    with open(variants_het_snp_vcf_file) as f:
        het_positions=set()
        for line in f:
            line = line.strip("\n").split("\t")
            if "#" not in line[0]:
                chr = line[0].replace(":", "_")
                pos = line[1]
                ref = line[3]
                alts = line[4].split(",")
                variants_info = ""
                mark = 0
                if len(ref) == 1:
                    # variants_info = chr+":"+pos+"="+ref
                    variant_alt = set()
                    for alt in alts:
                        if len(alt) == 1: #SNP or del of 1 base ("*")
                            variant_alt.add(alt)
                            mark += 1
                            # het_positions.add(chr+":"+pos)
                        elif len(alt) > len(ref):  # INS
                            # het_positions.add(chr+":"+pos)
                            continue
                        elif len(alt) < len(ref):  #DEL
                            #for curPos in range(int(line[1])+1,int(line[1])+len(ref)):
                            #   het_positions.add(chr+":"+str(curPos))
                            continue
                if mark != 0:
                    chr_pos = chr+":"+pos
                    variants_het_snp.setdefault(chr_pos, set())
                    variants_het_snp[chr_pos].add(ref)
                    variants_het_snp[chr_pos].update(variant_alt)
                    # het_positions.add(variants_info)

    het_positions_write = ""
    for var in het_positions:
        var = var.split("=")
        position = var[0]
        position = position.split(":")
        het_positions_write += position[0]+"\t"+position[1]+"\t"+var[1]+"\n"

    with open(variants_positions_file, "w") as f:
        for chr_pos, bases in variants_het_snp.items():
            chr_pos = chr_pos.split(":")
            bases_str = ",".join(bases)
            f.write(chr_pos[0]+"\t"+chr_pos[1]+"\t"+bases_str+"\n")

    return "Determined positions of heterozygous SNPs"

def get_snpVCF_from_vcf(vcf,snpvcf):
    if os.path.isfile(snpvcf):
        count = 0
        with open(snpvcf) as f:
            for line in f:
                line = line.strip("\n").split("\t")
                if "#" not in line[0]:
                    count += 1
                if count > 5:
                    return None
    
    # gatk 筛选snps失败，手动筛选，只根据ref和alt是否都是单碱基筛选
    vcf_items = ""
    with open(vcf) as f:
        for line in f:
            line = line.strip("\n")
            if "#" in line[0]:
                vcf_items += line + "\n"
            else:
                line = line.split("\t")
                ref = line[3]
                alts = line[4].split(",")
                mark = 0
                if len(ref) == 1:
                    mark = 1
                    for alt in alts:
                        if len(alt) != 1:
                            mark = 2
                            break
                if mark == 1:
                    vcf_items += "\t".join(line)+"\n"

    with open(snpvcf, "w") as f:
        f.write(vcf_items)
            
    return None


def parse_cigar(cigar, sequence, ref_start, ctg_name, variant_set):

    cigarPartRE = '(\d+)([DHIMNPSX=])'

    read_pos = 0
    match_start = False
    ref_pos = ref_start
    algn_var = dict()  # {key:pos  v:base}

    for match in re.finditer(cigarPartRE, cigar):
        n = int(match.group(1))
        algnment_type = match.group(2)
        if algnment_type == "S":
            read_pos += n
        elif algnment_type in "MX=":
            match_start = True
            matches = []
            for x in range(n):
                matches.append((ref_pos, sequence[read_pos]))
                ref_pos += 1
                read_pos += 1
            for match in matches:
                pos = match[0]
                base = match[1]
                pos_name = ctg_name+":" + str(pos)
                if pos_name in variant_set:
                    if base in variant_set[pos_name]:
                        algn_var[str(pos)] = base
        elif algnment_type == "I":
            insertion_sequence = ""
            for x in range(n):
                insertion_sequence += sequence[read_pos]
                read_pos += 1
            pos = ref_pos
        elif algnment_type == "D":
            if match_start:
                for x in range(n):
                    ref_pos += 1
    ref_end = ref_pos

    return (ref_end, algn_var)


def reduce_long_reads_to_SNPS_parallel(list_var):
    algn = list_var[0]
    align_reads = list_var[1]
    variant_set = list_var[2]

    read_name = algn[0]
    ctg_name = algn[1]
    ref_start = algn[2]
    mapq = algn[3]
    cigar = algn[4]
    sequence = algn[5]
    is_supplementary = algn[6]
    is_reversed = algn[7]
    ref_end, variant_seq = parse_cigar(cigar, sequence, ref_start, ctg_name, variant_set)

    if variant_seq == {}:
        return None

    align_reads.append(AlignedRead(read_name, variant_seq, mapq, is_supplementary, is_reversed, ref_start, ref_end, ctg_name))


def create_read_from_group(group, distance_threshold):
    """
        Merge multiple AlignedReads into a single Read.

        Pick supplementary reads that have the same orientation as the primary and with
        the distance at most distance_threshold from primary, find the set of variants that fully agree
        and return these variants as a read.

        If the list does not contain primary reads, then return None.
        If the list contains more than two primary alignments, return None and report a warning.
        """
    n_primary = 0
    primary = None
    for read in group:
        if not read.is_supplementary:
            n_primary += 1
            primary = read
    if primary is None:
        return None
    if n_primary > 2:
        # logger.warning(
        #         f"Read name {group[0].read_name!r} has more than two primary alignments."
        #     )
        return None

    ref_start = primary.ref_start
    ref_name = primary.ref_name
    variants = dict()
    skip = set()
    for read in group:
        if read.is_reversed != primary.is_reversed:
            continue
        if read.ref_name != ref_name:
            continue
        if primary.distance(read) > distance_threshold:
            continue
        ref_start = min(ref_start, read.ref_start)
        variant_seq = read.variant_seq
        for var_pos, var_base in variant_seq.items():
            if var_pos in variants:
                if variants[var_pos] != var_base:
                    skip.add(var_pos)
            else:
                variants[var_pos] = var_base

    var_list = []
    for pos, base in variants.items():
        if pos not in skip:
            var_list.append(pos+"="+base)

    return var_list


def group_reads(align_reads):
    groups = defaultdict(list)
    for align_read in align_reads:
        groups[align_read.read_name].append(align_read)

    distance_threshold = 10000

    read_dict = dict()
    for read_name, group in groups.items():
        len_group = len(group)
        if len_group == 1:
            variant_seq = group[0].variant_seq
            read_dict[read_name] = []
            for var_pos, var_base in variant_seq.items():
                var = var_pos + "=" + var_base
                read_dict[read_name].append(var)
        if len_group > 1:
            variants_list = create_read_from_group(group, distance_threshold)
            if variants_list == [] or variants_list == None:
                continue
            else:
                read_dict[read_name] = variants_list

    return read_dict


def load_bam_file(bam_file,use_supplementary=False,min_mapq=20):
    alignments = []
    with pysam.AlignmentFile(bam_file,"rb") as f:
        for line in f:
            # print(line)
            # flags = line.flag
            read_name = line.query_name
            ctg_name = line.reference_name
            # pysam读出的是0-index!sam和vcf都是1-index
            ref_start = line.reference_start + 1
            mapq = line.mapping_quality
            cigar = line.cigarstring
            sequence = line.query_sequence
            is_duplicate = line.is_duplicate
            is_unmapped = line.is_unmapped
            is_secondary = line.is_secondary
            is_reversed = line.is_reverse
            is_supplementary = line.is_supplementary

            if is_duplicate:
                continue
            if is_unmapped:
                continue
            if is_secondary:
                continue
            if mapq < min_mapq:
                continue
            if sequence == "*":
                continue
            if is_supplementary and not use_supplementary: continue

            ctg_name = ctg_name.replace(":", "_")
            mp_info = [read_name,ctg_name,ref_start,mapq,cigar,sequence,is_supplementary,is_reversed]
            alignments.append(mp_info)
    print("alignments: ",len(alignments))
    return alignments

def load_sam_file(sam_file,use_supplementary=False,min_mapq=20):
    alignments = []
    with open(sam_file, "r") as f:
        for line in f:
            line = line.strip("\n").split("\t")
            if len(line) >= 11:
                flags = int(line[1])
                is_unmapped = flags & 0x4
                is_secondary = flags & 0x100
                is_duplicate = flags & 0x400
                read_name = line[0]
                ctg_name = line[2]
                ref_start = int(line[3])
                cigar = line[5]
                is_supplementary = flags & 0x800
                is_reversed = flags & 0x16
                sequence = line[9]
                mapq = int(line[4])

                if is_duplicate:
                    continue
                if is_unmapped:
                    continue
                if is_secondary:
                    continue
                if mapq < min_mapq:
                    continue
                if sequence == "*":
                    continue
                if is_supplementary and not use_supplementary: continue

                ctg_name = ctg_name.replace(":", "_")
                mp_info = [read_name,ctg_name,ref_start,mapq,cigar,sequence,is_supplementary,is_reversed]
                line[2] = line[2].replace(":", "_")
                alignments.append(mp_info)
    print("alignments: ",len(alignments))
    return alignments

def reduce_long_reads_to_SNPS(bam_file, variants_positions_file, reduced_long_reads_file, use_supplementary, min_mapq):
    # variants positions
    variant_set = dict()
    with open(variants_positions_file, "r") as f:
        for line in f:
            line = line.strip("\n").split("\t")
            contig = line[0].replace(":", "_")
            contig = contig.replace("#", "_")  # 为了方便后面给cluster取名字: read_name#i
            start = line[1]
            ref_alts = line[2].split(",")
            variant_set[contig+":"+start] = set(ref_alts)

    # bam_file_name = bam_file.split("/")[-1]
    if bam_file.endswith("bam"):
        alignments = load_bam_file(bam_file, use_supplementary, min_mapq)
    elif bam_file.endswith("sam"):
        alignments = load_sam_file(bam_file, use_supplementary, min_mapq)
    else:
        print("bam files error!")

    manager = multiprocessing.Manager()
    align_reads = manager.list()

    with multiprocessing.Pool(processes=8) as pool:
        pool.map(reduce_long_reads_to_SNPS_parallel, [[algn, align_reads, variant_set] for algn in alignments])

    reads_dict = defaultdict(dict)
    for align_read in align_reads:
        variant_seq = align_read.variant_seq
        variant_seq_set = set()
        for pos,base in variant_seq.items():
            variant_seq_set.add(str(pos)+"="+base)
        reads_dict[align_read.ref_name][align_read.read_name] = variant_seq_set

    output_read = ""
    for ref_name, read in reads_dict.items():
        output_read += ">"+ref_name+"\n"
        for read_name, values in read.items():
            if len(values) == 0:
                continue
            algn_contents = read_name+"\t"+"\t".join(values)+"\n"
            output_read += algn_contents

    with open(reduced_long_reads_file, "w") as f:
        f.write(output_read)

    return "Generate reduced reads"