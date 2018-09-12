from functools import reduce
import random
from collections import Counter, defaultdict
import itertools as it

import numpy as np

from umi_tools.network import UMIClusterer


def trim_softclip(seq, cig):
    if cig[0][0] == 4:
        seq = seq[cig[0][1]:]
    if cig[-1][0] == 4:
        seq = seq[:-cig[-1][1]]
    return seq


def get_splice_juncs_and_sequence(aln, min_overhang=10):
    seq = trim_softclip(aln.seq, aln.cigar)
    leftmost = aln.reference_start
    curr_pos = 0
    splice_junctions = []
    for i, (cig_type, n_bases) in enumerate(aln.cigar):
        if cig_type in (0, 2):
            curr_pos += n_bases
        elif cig_type == 3:
            if aln.cigar[i - 1][1] >= min_overhang and \
               aln.cigar[i + 1][1] >= min_overhang:
                splice_seq = seq[curr_pos - min_overhang: curr_pos + min_overhang]
                if splice_seq:
                    splice_junctions.append((leftmost + curr_pos, leftmost + curr_pos + n_bases, splice_seq))
            curr_pos += n_bases
    return splice_junctions


def gene_kmers(seq, k):
    return set([seq[n: n + k] for n in range(0, len(seq) - k)])


def get_reference_seq(chrom, splice_junct, fasta, overhang=10):
    start, end = splice_junct
    seq = (fasta.fetch(chrom, start - overhang, start) + 
           fasta.fetch(chrom, end, end + overhang))
    return seq


def filter_unspliced_kmers(splice_junct_counts, query, fasta, overhang=10):
    contiguous_kmers = gene_kmers(fasta.fetch(*query), overhang * 2)
    filtered_splice_junct_counts = {}
    for pos, count in splice_junct_counts.items():
        seq = pos[2]
        if seq not in contiguous_kmers and 'N' not in seq:
            ref_seq = get_reference_seq(query[0], pos[:2], fasta, overhang)
            if ref_seq == seq:
                filtered_splice_junct_counts[pos] = count
    return filtered_splice_junct_counts


def get_splice_junctions_for_query(query, bam, min_overhang=10,
                                   filter_secondary=True,
                                   min_supporting_reads=2):
    '''Get all splice junctions from a bam file at a particular query position.
    '''
    splice_junctions = []
    for aln in bam.fetch(*query):
        if filter_secondary and aln.is_secondary:
            continue
        elif aln.reference_start > query[1] and aln.reference_end < query[2]:
            splice_junctions += get_splice_juncs_and_sequence(aln, min_overhang)
    splice_junct_counts = Counter(splice_junctions)
    splice_junct_counts = {
        pos: count for pos, count in splice_junct_counts.items()
        if count >= min_supporting_reads}
    return splice_junct_counts


def cluster_junctions_by_seq(splice_junct_counts, reference_name,
                             fasta, overhang=10, max_edit_distance=1):
    clustered_splice_juncts = defaultdict(set)
    clustered_splice_junct_counts = Counter()
    for (*pos, seq), count in splice_junct_counts.items():
        seq = seq.encode('utf-8')
        clustered_splice_juncts[seq].add((*pos, seq.decode()))
        clustered_splice_junct_counts[seq] += count
    clustered_seqs = UMIClusterer()(list(clustered_splice_juncts),
                                    clustered_splice_junct_counts,
                                    max_edit_distance)
    hamming_clustered_splice_junct_counts = {}
    for clust in clustered_seqs:
        pos_union = reduce(lambda x, y: x.union(y),
                           [clustered_splice_juncts[seq] for seq in clust])
        hamming_clustered_splice_junct_counts[frozenset(pos_union)] = sum(
            [splice_junct_counts[pos] for pos in pos_union])
    return hamming_clustered_splice_junct_counts


def sample_splice_junctions(junct_counts, idx,
                            min_supporting_samples=2):
    sample = [junct_counts[i] for i in idx]
    junct_sample_counts = Counter(it.chain(*sample))
    supported_splice_junctions = [
        pos for pos, count in junct_sample_counts.items()
        if count >= min_supporting_samples]
    supported_junct_counts = {pos: sum([s.get(pos, 0) for s in sample])
                              for pos in supported_splice_junctions}
    return supported_junct_counts