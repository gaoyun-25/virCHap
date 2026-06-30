import os
import sortedcontainers
from itertools import combinations
import multiprocessing
from collections import defaultdict
import shutil
import networkx as nx
try:
    from src import community_detection
except ImportError:
    import community_detection
import pysam

def get_similarity(common_pos, consensusA, consensusB, min_ovl, min_sim, min_overlap, min_ovl_len):
    if len(consensusA) < len(consensusB):
        setA = consensusA
        setB = consensusB
    else:
        setA = consensusB
        setB = consensusA
    if len(common_pos)/len(setA) >= min_ovl and len(common_pos) >= min_ovl_len or len(common_pos) >= min_overlap:
        common_set = setA & setB
        local_similarity = len(common_set)/len(common_pos)
        if local_similarity < min_sim:
            local_similarity = 0
    else:
        local_similarity = 0

    return local_similarity


def identify_merged_cluster_id(common_pos, cluster_idA, cluster_idB):
    merged_id = {}
    for pos in common_pos:
        baseA = cluster_idA[pos]
        baseB = cluster_idB[pos]
        base_possibilities = set().union(baseA["bases"].keys(), baseB["bases"].keys())
        base_new = {"bases": {}, "total": baseA["total"] + baseB["total"]}
        for base in base_possibilities:
            if base_new["total"] > 0:
                base_new["bases"][base] = (baseA["bases"].get(base,0)*baseA["total"] + baseB["bases"].get(base,0)*baseB["total"])/base_new["total"]
            else:
                base_new["bases"][base] = 0
        merged_id[pos] = base_new
    for pos in cluster_idA.keys():
        if pos not in common_pos:
            merged_id[pos] = cluster_idA[pos]
    for pos in cluster_idB.keys():
        if pos not in common_pos:
            merged_id[pos] = cluster_idB[pos]

    return merged_id


def identity_change_bool(merged_id, clustering_id, common_pos, max_ID):
    overall_change = 0
    n_all_seq = 0
    if len(common_pos) == 0:
        return True
    for pos in common_pos:
        demographics = clustering_id[pos]
        best_base = max(demographics["bases"].values())
        n_all_seq += merged_id[pos]["total"]
        for base, proportion in demographics["bases"].items():
            if proportion == best_base:
                if merged_id[pos]["bases"][base] - proportion < 0:
                    overall_change += abs(merged_id[pos]["bases"][base] - proportion)*merged_id[pos]["total"]
            else:
                if merged_id[pos]["bases"][base] - proportion > 0:
                    overall_change += abs(proportion - merged_id[pos]["bases"][base])*merged_id[pos]["total"]
    if n_all_seq == 0:
        return True
    if overall_change/n_all_seq > max_ID:
        return True
    else:
        return False


def identity_thef_bool(clusterA, clusterB, max_ID):
    common_pos = clusterA["positions"] & clusterB["positions"]
    cluster_idA = clusterA["cluster_id"]
    cluster_idB = clusterB["cluster_id"]
    merged_id = identify_merged_cluster_id(common_pos, cluster_idA, cluster_idB)
    if identity_change_bool(merged_id, cluster_idA, common_pos, max_ID) and identity_change_bool(merged_id, cluster_idB, common_pos, max_ID):
        return True
    else:
        return False


def find_pairs(clusters_info,similarity_index, min_sim, max_ID, banned_cluster_names):
    index = 0
    indexes = []
    best_pairs = []
    for comb in reversed(similarity_index):
        index -= 1
        nameA = comb[2]
        nameB = comb[3]
        if nameA in banned_cluster_names or nameB in banned_cluster_names:
            indexes.append(index)
            continue
        similarity = comb[0]
        if similarity >= min_sim:
            clusterA = clusters_info[nameA]
            clusterB = clusters_info[nameB]
            if not identity_thef_bool(clusterA, clusterB, max_ID):
                best_pairs.append([nameA, nameB])
                return best_pairs, indexes
            else:
                indexes.append(index)
    return best_pairs, indexes

def consensus(cluster_id):
    consensus = set()
    for pos in cluster_id.keys():
        highest_pct = max(cluster_id[pos]["bases"].values())
        for base, percentage in cluster_id[pos]["bases"].items():
            if percentage == highest_pct:
                consensus.add(pos + "=" + base)
                break
    return frozenset(consensus)


def initializeG(queue, out_queue, min_ovl, min_sim, min_overlap, min_ovl_len):
    local_infos = dict()
    while True:
        combination_batch = queue.get()
        if combination_batch == "end":
            out_queue.put(local_infos)
            return None
        for comb in combination_batch:
            v_i = comb[0]
            v_j = comb[1]
            cluster_i = clusters_infos[v_i]
            cluster_j = clusters_infos[v_j]
            local_infos.setdefault(v_i, {})
            local_infos.setdefault(v_j, {})
            common_positions = cluster_i["positions"] & cluster_j["positions"]
            if len(common_positions) >= min_ovl_len:
                similarity = get_similarity(common_positions, cluster_i["consensus"], cluster_j["consensus"], min_ovl, min_sim, min_overlap, min_ovl_len)
                if similarity == 0:
                    continue
                local_infos[v_i][v_j] = similarity
                local_infos[v_j][v_i] = similarity


def start_initializeG(queue, out_queue, clusters_infos_queue, clusters_infos, min_ovl, min_sim, min_overlap, min_ovl_len, thread_number):
    all_process = list()
    for thread in range(0, thread_number-1):
        initializeG_p = multiprocessing.Process(target=initializeG, args=((queue),(out_queue),(min_ovl),(min_sim),(min_overlap),(min_ovl_len)))
        initializeG_p.daemon = True
        initializeG_p.start()
        all_process.append(initializeG_p)
    initialize_merge_p = multiprocessing.Process(target=initialize_merge, args=((out_queue),(clusters_infos_queue),(clusters_infos),))
    initialize_merge_p.daemon = True
    initialize_merge_p.start()
    all_process.append(initialize_merge_p)

    return all_process


def initialize_merge(out_queue, clusters_infos_queue, clusters_infos):
    while True:
        local_dict = out_queue.get()
        if local_dict == "end":
            clusters_infos_queue.put(clusters_infos)
            return None
        for nameA in local_dict.keys():
            for nameB, similarity in local_dict[nameA].items():
                clusters_infos[nameA]["similarities"][nameB] = similarity


def compute_reads_similarity(reduced_reads, min_ovl, min_sim, min_overlap, min_ovl_len, thread_number):
    global clusters_infos
    clusters_infos = dict()
    i = 0
    for var_seq, read_name_set in reduced_reads.items():
        reads_name = set()
        reads_name.update(read_name_set)
        num_reads = len(read_name_set)
        positions = set([var.split("=")[0] for var in var_seq])
        cluster_id = dict()
        for var in var_seq:
            var = var.split("=")
            cluster_id.setdefault(var[0], {"bases": {var[1]: 1.0}, "total": num_reads})
        clusters_infos[str(i)] = {"reads_name": reads_name, "positions": positions, "cluster_id": cluster_id, "consensus": var_seq, "similarities": dict()}
        i += 1

    batch_size = 50000

    combination_queue = multiprocessing.Queue(maxsize=3*thread_number)
    local_similarity_dict_queue = multiprocessing.Queue()
    clusters_infos_queue = multiprocessing.Queue()

    all_process = start_initializeG(combination_queue, local_similarity_dict_queue, clusters_infos_queue, clusters_infos, min_ovl, min_sim, min_overlap, min_ovl_len, thread_number)
    all_initG = all_process[0:-1]
    init_merge = all_process[-1]

    current_batch = []
    i = 0
    for comb in combinations(clusters_infos.keys(), 2):
        current_batch.append(comb)
        if len(current_batch) > batch_size:
            combination_queue.put(current_batch)
            i += 1
            current_batch = []
    if len(current_batch) > 0:
        combination_queue.put(current_batch)
        current_batch=[]

    for thread in range(0, thread_number):
        combination_queue.put("end")

    for initG_p in all_initG:
        initG_p.join()

    local_similarity_dict_queue.put("end")

    clusters_infos = clusters_infos_queue.get()
    init_merge.join()

    clusters_info = dict()
    G = nx.Graph()
    for cidA, cluster in clusters_infos.items():
        clusters_info[cidA] = cluster
        for cidB, similarity in cluster["similarities"].items():
            G.add_edge(cidA, cidB, weight=similarity)
    del clusters_infos
    print("Read count after similarity computation: ",len(clusters_info))

    return G, clusters_info


def create_reverse_dict(dictionary):
    reverse_dict = {}
    for key, val in dictionary.items():
        reverse_dict[val] = key
    return reverse_dict


def create_directory(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)

def write_final_phased_output_to_files(phased_snps_file, phased_read_name_file, all_clusters):
    phasedSNP_text = ""
    clustered_read_name_text = ""
    for ctg_name, clusters in all_clusters.items():
        for cid, cluster_info in clusters.items():
            for snp in cluster_info["consensus"]:
                pos = snp.split("=")[0]
                base = snp.split("=")[1]
                phasedSNP_text += "\t".join([cid, ctg_name, pos, base]) + "\n"
            for read_name in cluster_info["reads_name"]:
                clustered_read_name_text += "\t".join([cid, read_name]) + "\n"
    with open(phased_snps_file, "w") as f:
        f.write(phasedSNP_text)

    with open(phased_read_name_file, "w") as f:
        f.write(clustered_read_name_text)

    return None


def compute_cids_similarity(cid_vertices, clusters_info,min_ovl, min_sim, min_overlap, min_ovl_len, thread_number,min_overlap_len_first,min_ovl_len_first,split_number,reference_length):
    global clusters_infos
    interval = round(reference_length/split_number)
    clusters_infos = dict()
    used_reads = set()
    count_cids = dict()
    for i in range(split_number):
        count_cids.setdefault(i,0)
    for cid, vers in cid_vertices.items():
        reads_name_tmp = set()
        cluster_id_tmp = dict()
        positions_tmp = set()
        for ver in vers:
            c_info_tmp = clusters_info[ver]
            used_reads.add(ver)
            reads_name_tmp = reads_name_tmp | c_info_tmp["reads_name"]
            common_pos = positions_tmp & c_info_tmp["positions"]
            positions_tmp = positions_tmp | c_info_tmp["positions"]
            cluster_id_tmp = identify_merged_cluster_id(common_pos, cluster_id_tmp, c_info_tmp["cluster_id"])
        consensus_tmp = consensus(cluster_id_tmp)
        clusters_infos[cid] = {"reads_name": reads_name_tmp, "positions": positions_tmp, "cluster_id": cluster_id_tmp, "consensus": consensus_tmp, "similarities": dict(),"overlaps":set()}
        tmp_pos = set()
        for pos in positions_tmp:
            tmp_pos.add(int(pos))
        min_pos = min(tmp_pos)
        max_pos = max(tmp_pos)
        left = int(min_pos/interval)
        right = int(max_pos/interval)
        diff = right-left
        for i in range(diff):
            count_cids[left+i] += 1

    empty_regions = set()
    for i, count in count_cids.items():
        if count <= 2:
            empty_regions.add(i)

    if empty_regions != set():
        used_min_overlap_len = min_overlap
        used_min_ovl_len = min_ovl_len
        for ver, info in clusters_info.items():
            if ver in used_reads:
                continue
            tmp_positions = info["positions"]
            tmp_pos = set()
            for pos in tmp_positions:
                tmp_pos.add(int(pos))
            min_pos = min(tmp_pos)
            max_pos = max(tmp_pos)
            left = int(min_pos/interval)
            right = int(max_pos/interval)
            diff = right-left
            span_regions = {right}
            for i in range(diff):
                span_regions.add(i+left)
            if span_regions & empty_regions != set():
                info["overlaps"] = set()
                clusters_infos[ver] = info
    else:
        used_min_overlap_len = min_overlap_len_first
        used_min_ovl_len = min_ovl_len_first
    print("The number of clusters after rescuing reads:",len(clusters_infos))

    batch_size = 50000

    combination_queue = multiprocessing.Queue(maxsize=3*thread_number)
    local_similarity_dict_queue = multiprocessing.Queue()
    clusters_infos_queue = multiprocessing.Queue()

    all_process = start_initializeG(combination_queue, local_similarity_dict_queue, clusters_infos_queue, clusters_infos, min_ovl, min_sim, used_min_overlap_len, used_min_ovl_len, thread_number)
    all_initG = all_process[0:-1]
    init_merge = all_process[-1]

    current_batch = []
    i = 0
    for comb in combinations(clusters_infos.keys(), 2):
        current_batch.append(comb)
        if len(current_batch) > batch_size:
            combination_queue.put(current_batch)
            i += 1
            current_batch = []
    if len(current_batch) > 0:
        combination_queue.put(current_batch)
        current_batch=[]

    for thread in range(0, thread_number):
        combination_queue.put("end")

    for initG_p in all_initG:
        initG_p.join()

    local_similarity_dict_queue.put("end")

    clusters_infos = clusters_infos_queue.get()

    init_merge.join()

    clusters_info = dict()
    for cidA, cluster in clusters_infos.items():
        clusters_info[cidA] = cluster
        clusters_info[cidA]["overlaps"] = set(cluster["similarities"])

    del clusters_infos

    return clusters_info,used_min_overlap_len, used_min_ovl_len

def phasing_contig(reduced_reads,ctg_name,thread_number, max_ID, min_len, min_ovl, min_sim, min_overlap_len, min_ovl_len, min_sim_first,min_abundance,min_overlap_len_first,min_ovl_len_first,split_number,reference_length):
    print("Start computing read similarity...")
    G, clusters_info = compute_reads_similarity(reduced_reads, min_ovl, min_sim_first, min_overlap_len_first,min_ovl_len_first, thread_number)

    map_id = {vid:i for i,vid in enumerate(G.nodes)}
    rev_map_id = create_reverse_dict(map_id)
    G = nx.relabel_nodes(G, map_id)
    print("Start label propagation clustering.")
    partition = community_detection.find_communities(G)
    # partition node:cid
    cid_vertices = {}
    for vid, cid in partition.items():
        cid_vertices.setdefault(str(cid),set())
        cid_vertices[str(cid)].add(rev_map_id[vid])
    print("The number of clusters after LPA: ", len(cid_vertices))

    split_number = int(split_number)
    reference_length = int(reference_length)
    clusters_info_merge,used_min_overlap_len, used_min_ovl_len = compute_cids_similarity(cid_vertices, clusters_info,min_ovl, min_sim, min_overlap_len, min_ovl_len, thread_number,min_overlap_len_first,min_ovl_len_first,split_number,reference_length)
    clusters_info = ""

    print("Start merging clusters.")
    clusters_info_merge = connect_clusters(clusters_info_merge, min_ovl, min_sim, max_ID, used_min_overlap_len, used_min_ovl_len)
    len_clusters_info = len(clusters_info_merge)
    print("The number of clusters after merging: ",len_clusters_info)

    tmp_clusters = dict()
    total_reads = 0
    for cid, cluster_info in clusters_info_merge.items():
        if len(cluster_info["reads_name"]) <= min_len:
            continue
        total_reads += len(cluster_info["reads_name"])
        tmp_clusters[cid] = {"reads_name": cluster_info["reads_name"], "consensus": cluster_info["consensus"]}

    clusters_info_merge = ""

    clusters = dict()
    for cid, cid_info in tmp_clusters.items():
        read_ids = cid_info["reads_name"]
        abundance = len(read_ids)/total_reads
        if abundance < min_abundance:
            continue
        new_cid = cid + "_" +str(len(read_ids))+"_"+ str(abundance)+",Contig_name:"+ctg_name
        clusters[new_cid] = cid_info
    print("The number of clusters after filtering low abundance", len(clusters))
    cluster_info = ""
    G = ""

    return clusters


def connect_clusters(clusters_info_merge, min_ovl, min_sim, max_ID, min_overlap, min_ovl_len):
    similarity_list = list()

    for cidA, info in clusters_info_merge.items():
        for cidB, sim in info["similarities"].items():
            common_pos = info["positions"] & clusters_info_merge[cidB]["positions"]
            similarity_list.append([sim, len(common_pos), cidA, cidB])

    similarity_index = sortedcontainers.SortedList()
    similarity_index.update(similarity_list)
    banned_cluster_names = set()
    pairs, indexes = find_pairs(clusters_info_merge, similarity_index, min_sim, max_ID, banned_cluster_names)
    i = 0
    while pairs != []:
        for index in reversed(indexes):
            similarity_index.pop(index)
        for pair in pairs:
            banned_cluster_names.update(pair)
            first = pair[0]
            second = pair[1]
            new_name = first.split("#")[0] + "#" + str(i)
            new_reads_name = clusters_info_merge[first]["reads_name"] | clusters_info_merge[second]["reads_name"]
            common_positions = clusters_info_merge[first]["positions"] & clusters_info_merge[second]["positions"]
            cluster_id_1 = clusters_info_merge[first]["cluster_id"]
            cluster_id_2 = clusters_info_merge[second]["cluster_id"]
            new_cluster_id = identify_merged_cluster_id(common_positions, cluster_id_1, cluster_id_2)
            new_consensus = consensus(new_cluster_id)
            new_positions = set([var.split("=")[0] for var in new_consensus])
            new_overlaps = clusters_info_merge[first]["overlaps"] | clusters_info_merge[second]["overlaps"]
            new_overlaps.remove(first)
            new_overlaps.remove(second)
            del clusters_info_merge[first]
            del clusters_info_merge[second]
            for neighbor in new_overlaps:
                clusterB = clusters_info_merge[neighbor]
                new_common_pos = new_positions & clusterB["positions"]
                clusters_info_merge[neighbor]["overlaps"].discard(first)
                clusters_info_merge[neighbor]["overlaps"].discard(second)
                clusters_info_merge[neighbor]["overlaps"].add(new_name)
                if len(new_common_pos) >= min_ovl_len:
                    similarity = get_similarity(new_common_pos, new_consensus, clusterB["consensus"], min_ovl, min_sim, min_overlap, min_ovl_len)
                    if similarity >= min_sim:
                        similarity_index.add([similarity, len(new_common_pos), new_name, neighbor])
            clusters_info_merge[new_name] = {"reads_name": new_reads_name, "positions": new_positions, "cluster_id": new_cluster_id, "consensus": new_consensus, "similarities": dict(),"overlaps": new_overlaps}
            i += 1
        pairs, indexes = find_pairs(clusters_info_merge, similarity_index, min_sim, max_ID, banned_cluster_names)

    return clusters_info_merge

def generate_consensus_fa(bam_file, total_clusters, phased_snps_file, phased_read_name_file, phased_consensus_fa, use_supplementary=False,min_mapq=20, min_base_cov = 2,min_base_threshold=0.66):
    read_to_cids = dict()
    for ctg_name, ctg_clusters in total_clusters.items():
        for cid, cid_info in ctg_clusters.items():
            read_ids = cid_info["reads_name"]
            for read_id in read_ids:
                read_to_cids[read_id] = cid

    write_final_phased_output_to_files(phased_snps_file, phased_read_name_file, total_clusters)

    all_reads = set(read_to_cids.keys())
    cluster_records = defaultdict(list)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_duplicate or read.is_unmapped or read.is_secondary or read.query_sequence == "*" or read.mapping_quality < min_mapq:
                continue
            if read.is_supplementary and not use_supplementary: continue
            read_name = read.query_name
            if read_name in all_reads:
                cid = read_to_cids[read_name]
                cluster_records[cid].append(read)

    consensus_seqs = dict()
    for cid, reads in cluster_records.items():
        if not reads:
            continue
        ref_starts = []
        ref_ends = []
        for read in reads:
            ref_starts.append(read.reference_start)
            ref_ends.append(read.reference_end)
        
        if not ref_starts:
            continue

        ref_start = min(ref_starts)
        ref_end = max(ref_ends)
        region_length = ref_end - ref_start
        projected_info = dict()
        for read in reads:
            aligned_pairs = read.get_aligned_pairs(matches_only=False)
            for qpos, rpos in aligned_pairs:
                if rpos is None:
                    continue

                rpos = rpos - ref_start
                projected_info.setdefault(rpos, {"A":0,"C":0,"G":0,"T":0})
            
                if qpos is None:
                    continue
                else:
                    base = read.query_sequence[qpos]
                    if base in {"A", "C", "G", "T"}:
                        projected_info[rpos][base] += 1

        # consensus
        consensus_seq = ["N"] * region_length
        for pos, counts in projected_info.items():
            max_count = max(counts.values())
            for b,c in counts.items():
                if c == max_count:
                    max_base = b
            total_counts = sum(counts.values())
            if total_counts < min_base_cov:
                consensus_seq[pos] = "N"
                continue
            if (max_count / total_counts) < min_base_threshold:
                consensus_seq[pos] = "N"
                continue
            consensus_seq[pos] = max_base
        consensus_seqs[cid] = consensus_seq

    with open(phased_consensus_fa, "w") as f:
        for cid, consensus_seq in consensus_seqs.items():
            f.write(">"+cid+"\n")
            f.write("".join(consensus_seq)+"\n")


def phasing(reduced_reads_file,bam_file, thread_number, out_dir, out_name, max_ID, min_len, min_ovl, min_sim, min_overlap_len, min_ovl_len, min_sim_first,min_abundance,min_overlap_len_first,min_ovl_len_first,split_number,reference_lengths):
    outdir_name = "phased"+"_"+str(min_ovl)+"_"+str(min_sim)+"_"+str(max_ID)+"_"+str(min_len)+"_"+str(min_overlap_len)+"_"+str(min_ovl_len)+"_"+str(min_sim_first)+"_"+str(min_abundance)+"_"+str(min_overlap_len_first)+"_"+str(min_ovl_len_first)
    thread_number = int(thread_number)
    max_ID = float(max_ID)
    min_ovl = float(min_ovl)
    min_sim = float(min_sim)
    min_overlap_len = int(min_overlap_len)
    min_ovl_len = int(min_ovl_len)
    min_len = int(min_len)
    min_overlap_len_first = int(min_overlap_len_first)
    min_ovl_len_first = int(min_ovl_len_first)
    min_abundance = float(min_abundance)
    min_sim_first = float(min_sim_first)
    reads_dict = dict()
    num = 0
    filter_min_ovl_len = min(min_ovl_len, min_ovl_len_first)
    with open(reduced_reads_file, "r") as f:
        for line in f:
            if line[0] == ">":
                line=line.strip()
                ctg_name = line[1:]
                reads_dict[ctg_name] = dict()
            else:
                line = line.strip("\n").split("\t")
                reduced_read = frozenset(line[1:])
                if len(reduced_read) < filter_min_ovl_len:
                    continue
                num += 1
                reads_dict[ctg_name].setdefault(reduced_read, set())
                reads_dict[ctg_name][reduced_read].add(line[0])
    print("The number of reads: ", num)
    print("The number of unique reads: ",len(reads_dict[ctg_name]))

    phased_path = os.path.join(out_dir, "Phased")
    phased_snps_file = os.path.join(phased_path, outdir_name+"_variants.tsv")
    phased_read_name_file = os.path.join(phased_path, outdir_name+"_clustered_read_name.tsv")
    phased_consensus_fa = os.path.join(phased_path, outdir_name+"_consensus.fasta")
    if os.path.exists(phased_snps_file):
        os.remove(phased_snps_file)
    if os.path.exists(phased_read_name_file):
        os.remove(phased_read_name_file)

    total_clusters = dict()
    for ctg_name, reduced_reads in reads_dict.items():
        total_clusters.setdefault(ctg_name, dict())
        clusters = phasing_contig(reduced_reads,ctg_name,thread_number, max_ID, min_len, min_ovl, min_sim, min_overlap_len, min_ovl_len, min_sim_first,min_abundance,min_overlap_len_first,min_ovl_len_first,split_number,reference_lengths[ctg_name])
        total_clusters[ctg_name].update(clusters)

    generate_consensus_fa(bam_file, total_clusters, phased_snps_file, phased_read_name_file, phased_consensus_fa)

    return "Phaing ends"

