#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2026.01.05
Version: 1.3.0

This is a script designed to construct a high-confidence transcriptome from a merged collection of putative
transcripts (from running GffCompare together with GENCODE annotations). This script requires the following 
input files:
    1. The output file prefix used when running GffCompare
    2. A GTF file containing GENCODE annotations (the same one used when running GffCompare)
        * When running GffCompare, this file must be the first file provided as input!
    3. A FASTA file containing the reference genome sequence
    4. A BED file with refTSS v4.1 peaks (fifth column describes the number of samples with TPM > 1)
    5. A BED file with polyASite v2.0 peaks (fifth column describes the number of samples with TPM > 1)

This script will generate a GTF file containing the coordinates of transcripts in our high-confidence transcriptome.

Note: This script is intended for analyzing human transcriptomes!
'''

# =====================================================================================================================
#                                                       PACKAGES 
# =====================================================================================================================

# Load required packages
from collections import Counter
from datetime import datetime
import argparse, csv
import networkx as nx
import pandas as pd
import pysam, tabix

# =====================================================================================================================
#                                                   HELPER FUNCTIONS
# =====================================================================================================================

def outlog(myString):
    '''
    This is a function to print some output message string (myString) with the date and time
    '''
    print('[', datetime.now().strftime("%Y-%m-%d %H:%M:%S"), '] ', myString, sep = '', flush = True)

def PullFeature(infoString, feature):
    '''
    This is a function designed to pull out the value of a particular feature in a GTF info field
    '''
    return next(iter([item.split('"')[1] for item in infoString.split(';') if feature in item]), '.')

def CheckOverlap(tup1, tup2):
    '''
    This is a function that checks whether two tuples overlap
    '''
    return min(tup1[1], tup2[1]) >= max(tup1[0], tup2[0])

def PerfectSublist(list1, list2):
    '''
    This is a function that checks if list1 is a perfect sublist of list2
    '''
    if len(list1) <= len(list2):
        for i in range(len(list2) - len(list1) + 1):
            if list2[i:(i+len(list1))] == list1:
                return True
    
    return False

# =====================================================================================================================
#                                                    MAIN FUNCTIONS
# =====================================================================================================================

def FindBubbles(splice_graph, curr_path):
    bubbles, bubble_start = [], None
    for node in curr_path:
        if splice_graph.in_degree(node) == 2:
            if bubble_start is not None:
                bubble_end = node
                bubble_paths = list(nx.all_simple_paths(splice_graph, bubble_start, bubble_end))
                if PerfectSublist(bubble_paths[0], curr_path):
                    bubbles.append(bubble_paths[1])
                else:
                    bubbles.append(bubble_paths[0])
                bubble_start = None
        if splice_graph.out_degree(node) == 2:
            bubble_start = node
    return bubbles

def FormatGTF(chrom, transcript_chain, strand, geneID, txID):
    '''
    This is a function that returns GTF entries for a given transcript
    '''
    infoString = 'gene_id "' + geneID + '"; transcript_id "' + txID + '";'
    currlines = [[chrom, 'StringTie', 'transcript', transcript_chain[0], transcript_chain[-1], '.', strand, '.', infoString]]
    
    # Iterate over exons of transcript
    txExons = list(zip(transcript_chain[0::2], transcript_chain[1::2]))
    for idx, val in enumerate(txExons):
        currlines.append([chrom, 'StringTie', 'exon', val[0], val[1], '.', strand, '.', infoString + ' exon_number "' + str(idx+1) + '";'])
    
    return currlines

def FindPeak(peakfile, chrom, coord, dist, threshold, strand):
    '''
    This is a function that checks whether [coord - dist, coord + dist] harbors any peaks on the same strand
    meeting a specific peak score threshold (assumes that chrom is one of chr1-22,X,Y)
    '''
    if chrom != 'chrM':
        return len([item for item in peakfile.query(chrom, coord - dist, coord + dist) if int(item[4]) >= threshold and item[5] == strand]) > 0
    else:
        return False

def AdjustTerminalEnds(group, reftssbed, polyabed, chrom, strand):
    '''
    This is a function that adjusts the start/end coordinates of transcripts represented in a group
    This function also checks how reliable these ends are
    '''
    updated_group = group.copy()

    # Pull out first exons for transcripts in updated_group
    updated_group['first_exon'] = updated_group['transcript_chain'].apply(lambda x: tuple(x[:2]))
    updated_group = updated_group.sort_values(by = 'first_exon').reset_index(drop = True)

    # Cluster the first exons
    first_exon_clusters, curr_cluster = [], []
    for idx, exon in enumerate(updated_group['first_exon'].values):
        if len(curr_cluster) == 0:
            curr_cluster = [exon]
        else:
            if CheckOverlap(curr_cluster[0], exon):
                curr_cluster.append(exon)
            else:
                if len(first_exon_clusters) == 0:
                    first_exon_clusters += [1] * len(curr_cluster)
                else:
                    first_exon_clusters += [first_exon_clusters[-1] + 1] * len(curr_cluster)
                curr_cluster = [exon]

    if len(first_exon_clusters) > 0:
        first_exon_clusters += [first_exon_clusters[-1] + 1] * len(curr_cluster)
    else:
        first_exon_clusters = [1] * len(curr_cluster)
    updated_group['first_exon_clusters'] = first_exon_clusters

    updated_tss, trust_start = [], []
    for cluster, fegroup in updated_group.groupby('first_exon_clusters'):
        if any(fegroup['gencode_basic'].values):
            # Default to the most upstream basic TSS
            pivot = min(fegroup[fegroup['gencode_basic']]['first_exon'].apply(lambda x: x[0]).values)
            updated_tss += [pivot if pivot < fe[1] else fe[0] for fe in fegroup['first_exon'].values]
            trust_start += [True] * fegroup.shape[0]
        else:
            # Check if any of the TSS's overlap a peak
            if strand == '+':
                fegroup['near_peak'] = fegroup['first_exon'].apply(lambda x: FindPeak(reftssbed, chrom, x[0], 100, 5, strand))
            else:
                fegroup['near_peak'] = fegroup['first_exon'].apply(lambda x: FindPeak(polyabed, chrom, x[0], 100, 5, strand))
            
            if any(fegroup['near_peak']):
                # Default to the most upstream TSS classified as "near peak"
                pivot = min(fegroup[fegroup['near_peak']]['first_exon'].apply(lambda x: x[0]).values)
                updated_tss += [pivot if pivot < fe[1] else fe[0] for fe in fegroup['first_exon'].values]
                trust_start += [True] * fegroup.shape[0]
            else:
                # Assume that the TSS's are not trustworthy
                updated_tss += [-1] * fegroup.shape[0]
                trust_start += [False] * fegroup.shape[0]

    updated_group['updated_tss'], updated_group['trust_start'] = updated_tss, trust_start

    # Pull out last exons for transcripts in updated_group
    updated_group['last_exon'] = updated_group['transcript_chain'].apply(lambda x: tuple(x[-2:]))
    updated_group = updated_group.sort_values(by = 'last_exon', ascending = False).reset_index(drop = True)

    # Cluster the last exons
    last_exon_clusters, curr_cluster = [], []
    for idx, exon in enumerate(updated_group['last_exon'].values):
        if len(curr_cluster) == 0:
            curr_cluster = [exon]
        else:
            if CheckOverlap(curr_cluster[0], exon):
                curr_cluster.append(exon)
            else:
                if len(last_exon_clusters) == 0:
                    last_exon_clusters += [1] * len(curr_cluster)
                else:
                    last_exon_clusters += [last_exon_clusters[-1] + 1] * len(curr_cluster)
                curr_cluster = [exon]

    if len(last_exon_clusters) > 0:
        last_exon_clusters += [last_exon_clusters[-1] + 1] * len(curr_cluster)
    else:
        last_exon_clusters = [1] * len(curr_cluster)
    updated_group['last_exon_clusters'] = last_exon_clusters

    updated_tes, trust_end = [], []
    for cluster, legroup in updated_group.groupby('last_exon_clusters'):
        if any(legroup['gencode_basic'].values):
            # Default to the most downstream basic TES
            pivot = max(legroup[legroup['gencode_basic']]['last_exon'].apply(lambda x: x[-1]).values)
            updated_tes += [pivot if le[0] < pivot else le[1] for le in legroup['last_exon'].values]
            trust_end += [True] * legroup.shape[0]
        else:
            # Check if any of the TSS's overlap a peak
            if strand == '+':
                legroup['near_peak'] = legroup['last_exon'].apply(lambda x: FindPeak(polyabed, chrom, x[-1], 100, 5, strand))
            else:
                legroup['near_peak'] = legroup['last_exon'].apply(lambda x: FindPeak(reftssbed, chrom, x[-1], 100, 5, strand))
            
            if any(legroup['near_peak']):
                # Default to the most upstream TSS classified as "near peak"
                pivot = max(legroup[legroup['near_peak']]['last_exon'].apply(lambda x: x[-1]).values)
                updated_tes += [pivot if le[0] < pivot else le[1] for le in legroup['last_exon'].values]
                trust_end += [True] * legroup.shape[0]
            else:
                # Assume that the TES's are not trustworthy
                updated_tes += [-1] * legroup.shape[0]
                trust_end += [False] * legroup.shape[0]

    updated_group['updated_tes'], updated_group['trust_end'] = updated_tes, trust_end

    updated_group['transcript_chain'] = updated_group.apply(lambda x: ([x['updated_tss'] if x['updated_tss'] != -1 else x['transcript_chain'][0]] +
        x['transcript_chain'][1:-1] + [x['updated_tes'] if x['updated_tes'] != -1 else x['transcript_chain'][-1]]), axis = 1)
    updated_group = updated_group.drop(['first_exon', 'first_exon_clusters', 'updated_tss', 'last_exon', 'last_exon_clusters', 'updated_tes'], axis = 1)
    updated_group = updated_group.sort_values(by = ['gencode_canonical', 'gencode_basic', 'gencode_id'], ascending = False).reset_index(drop = True)

    return updated_group

def RemoveSubsetTranscripts(updated_group):
    '''
    This is a function that removes non-basic transcripts that are fully contained in at least one basic transcript
    '''
    basicChains = updated_group[updated_group['gencode_basic']].values
    nonBasicChains = updated_group[~updated_group['gencode_basic']].values
    blacklist = set()

    for nchain in nonBasicChains:
        for bchain in basicChains:
            # Check that nchain[4][1:-1] is a subset of bchain[4]
            if PerfectSublist(nchain[4][1:-1], bchain[4]):
                leftIdx = bchain[4].index(nchain[4][1])
                rightIdx = bchain[4].index(nchain[4][-2])
                if nchain[4][0] >= bchain[4][leftIdx-1] and nchain[4][-1] <= bchain[4][rightIdx+1]:
                    blacklist.add(nchain[1])
    
    return updated_group[~updated_group['transcript_id'].isin(blacklist)]

def FilterTranscripts(group, genomeFasta, reftssbed, polyabed):
    '''
    This is a function that filters group for high-confidence transcripts
    '''
    # Pull out group ID, chromosome, and strand information for group (based on first row of group)
    group_id, chrom, strand = group.iloc[0]['gene_id'], group.iloc[0][0], group.iloc[0][6]

    # Adjust TSS and TES for transcripts in group
    updated_group = AdjustTerminalEnds(group, reftssbed, polyabed, chrom, strand)

    # Remove non-basic transcripts that are mono-exonic
    # Rationale: We think mono-exonic transcripts are highly suspicious in long-read RNA-seq data unless if they're annotated
    updated_group = updated_group[~(~updated_group['gencode_basic'] & (updated_group['transcript_chain'].apply(len) == 2))]

    # Remove non-basic transcripts that are fully contained in at least one basic transcript (yes, this can still happen)
    # At this point we can assume that any non-basic transcript has at least one intron since we removed non-basic mono-exonic transcripts
    updated_group = RemoveSubsetTranscripts(updated_group)

    # Collect annotated splice sites and record the corresponding gene assignment
    anno_splice_sites = dict()

    for item in updated_group[updated_group['gencode_id'] != '.'].values:
        # Iterate over splice junctions and corresponding splice sites
        for junc in zip(item[4][1:-1][0::2], item[4][1:-1][1::2]):
            anno_splice_sites[junc[0]] =  anno_splice_sites.get(junc[0], set()) | {item[5].split('|')[0]}
            anno_splice_sites[junc[1]] =  anno_splice_sites.get(junc[1], set()) | {item[5].split('|')[0]}
    
    # Iterate over transcripts and keep track of "bad" introns that do not satisfy the following criteria:
    #   * If an intron has two novel splice sites, it cannot be retained in any other transcript in the group
    #       * This requirement protects against spurious exitrons corresponding to library preparation artifacts
    #       * UPDATE (as of 12/16/2025): We won't allow introns with two novel splice sites, ever.
    #   * If there is a novel splice site, it has a canonical dinucleotide motif (GT for donors, AG for acceptors)
    #   * If the intron features two annotated splice sites, the corresponding splice sites must be assigned to the same gene
    #       * This requirement protects against inclusion of fusion transcripts that may correspond to library preparation artifacts

    bad_introns, gene_assignment = [], []

    # Harvest all exons represented in updated_group
    unique_exons = {exon for tx_chain in updated_group['transcript_chain'].values for exon in zip(tx_chain[0::2], tx_chain[1::2])}

    for item in updated_group.values:
        if item[5] != '.':
            # We are dealing with an annotated transcript (no intron is bad)
            bad_introns.append([])
            gene_assignment.append(item[5].split('|')[0])

        else:
            # We are dealing with a novel transcript (some introns may be "bad")
            curr_bad_introns, num_anno_sites, candidate_genes = [], 0, []
            
            # We are guaranteed to be dealing with a multi-exonic transcript given our prior filtering strategy
            for intron in zip(item[4][1:-1][0::2], item[4][1:-1][1::2]):
                if intron[0] not in anno_splice_sites and intron[1] not in anno_splice_sites:
                    # Intron has two novel splice sites
                    # UPDATE (as of 12/16/2025): We will just blacklist this intron, period.
                    curr_bad_introns.append(intron)
                    # Round 1: Make sure that the intron has canonical dinucleotides
                    # dnt1 = genomeFasta.fetch(chrom, intron[0], intron[0]+2).upper()
                    # dnt2 = genomeFasta.fetch(chrom, intron[1]-3, intron[1]-1).upper()

                    # if (dnt1 + dnt2) != ('GTAG' if strand == '+' else 'CTAC'):
                    #     curr_bad_introns.append(intron)
                    # else:
                    #     # Make sure that this intron is not physically contained in any of the exons in unique_exons
                    #     if any([exon[0] <= intron[0] and intron[1] <= exon[1] for exon in unique_exons]):
                    #         curr_bad_introns.append(intron)

                elif intron[0] not in anno_splice_sites:
                    dnt = genomeFasta.fetch(chrom, intron[0], intron[0]+2).upper()
                    if dnt != ('GT' if strand == '+' else 'CT'):
                        curr_bad_introns.append(intron)
                    num_anno_sites += 1
                    candidate_genes += [gid for gid in anno_splice_sites[intron[1]]]
                elif intron[1] not in anno_splice_sites:
                    dnt = genomeFasta.fetch(chrom, intron[1]-3, intron[1]-1).upper()
                    if dnt != ('AG' if strand == '+' else 'AC'):
                        curr_bad_introns.append(intron)
                    num_anno_sites += 1
                    candidate_genes += [gid for gid in anno_splice_sites[intron[0]]]
                else:
                    num_anno_sites += 2
                    candidate_genes += [gid for site in intron for gid in anno_splice_sites[site]]
                    if len(anno_splice_sites[intron[0]] & anno_splice_sites[intron[1]]) == 0:
                        curr_bad_introns.append(intron)
            
            bad_introns.append(curr_bad_introns)

            # Determine optimal gene assignment for transcript
            gene_assignment.append(next(iter([k for k, v in dict(Counter(candidate_genes)).items() if v == num_anno_sites]), '.'))

    updated_group['bad_introns'], updated_group['gene_assignment'] = bad_introns, gene_assignment

    # Curate a set of high-confidence transcripts from updated_group. These transcripts must meet the following requirements:
    #   * Transcript has complete 5' and 3' ends and has no bad introns
    #   * Transcript must also have a defined gene assignment
    pass_group = updated_group[updated_group['trust_start'] & updated_group['trust_end'] & (updated_group['bad_introns'].apply(len) == 0) & 
        (updated_group['gene_assignment'] != '.')].reset_index(drop = True)
    
    # Pull out all exons of transcripts in pass_group
    pass_exons = {exon for tx_chain in pass_group['transcript_chain'].values for exon in zip(tx_chain[0::2], tx_chain[1::2])}

    # Rescue novel transcript isoforms that:
    #   * Have a defined gene assignment and do not have any bad introns
    #   * Transcript was flagged as having either an incomplete 5' or 3' end
    # The idea here is if their transcript ends were marked as bad but they do not overlap any exons in pass_exons, then we can still keep them around
    # Rationale: If a multi-exonic transcript is fragmented, the fragmented end would be hanging off of an exon that is trustworthy

    rescue_group = updated_group[~updated_group['transcript_id'].isin(set(pass_group['transcript_id']))].reset_index(drop = True)
    rescue_group['rescue_start'] = rescue_group.apply(lambda x: True if x['trust_start'] else not any([CheckOverlap(item, tuple(x['transcript_chain'][:2])) for item in pass_exons]), axis = 1)
    rescue_group['rescue_end'] = rescue_group.apply(lambda x: True if x['trust_end'] else not any([CheckOverlap(item, tuple(x['transcript_chain'][-2:])) for item in pass_exons]), axis = 1)
    rescue_group = rescue_group[rescue_group['rescue_start'] & rescue_group['rescue_end'] & (rescue_group['bad_introns'].apply(len) == 0) & 
        (rescue_group['gene_assignment'] != '.')].reset_index(drop = True).drop(['rescue_start', 'rescue_end'], axis = 1)
    pass_group = pd.concat([pass_group, rescue_group]).reset_index(drop = True)

    # Rescue potential full-length transcript isoforms not reported by StringTie using the following approach:
    #   For each non-basic transcript isoform (pass or no pass), compare its structure with the structure of the canonical isoform of the corresponding gene
    #   For each "bubble" in the comparison, check that the bubble is not already part of any basic transcript. If its inclusion in the context
    #   of the canonical isoform yields an isoform that is not part of pass_group, keep it. But make sure that this bubble does not have any bad introns!
    
    rescue_paths = dict()

    # Iterate over all canonical transcripts represented in pass_group
    for canonTx in pass_group[pass_group['gencode_canonical']].values:
        canon_path = ([str(canonTx[4][0]) + '-R'] + [str(canonTx[4][idx+1]) + ('-D' if idx % 2 == 0 else '-A') for idx in range(len(canonTx[4][1:-1]))] + 
            [str(canonTx[4][-1]) + '-S'])
        canon_sg = nx.DiGraph()
        canon_sg.add_edges_from([(canon_path[idx], canon_path[idx+1]) for idx in range(len(canon_path)-1)])

        # Record basic and non-basic paths for the corresponding gene
        basic_paths = pass_group[pass_group['gencode_basic'] & (pass_group['gene_assignment'] == canonTx[-1])]['transcript_chain'].values
        non_basic_paths = updated_group[~updated_group['gencode_basic'] & (updated_group['gene_assignment'] == canonTx[-1])].values
        non_basic_pass_paths = pass_group[~pass_group['gencode_basic'] & (pass_group['gene_assignment'] == canonTx[-1])]['transcript_chain'].values

        # Iterate over each path in non_basic_paths
        for non_basic_path in non_basic_paths:
            updated_path = ([str(non_basic_path[4][0]) + '-R'] + [str(non_basic_path[4][idx+1]) + ('-D' if idx % 2 == 0 else '-A') for idx in range(len(non_basic_path[4][1:-1]))] + 
                [str(non_basic_path[4][-1]) + '-S'])
            updated_sg = nx.DiGraph()
            updated_sg.add_edges_from([(updated_path[idx], updated_path[idx+1]) for idx in range(len(updated_path)-1)])
            combined_sg = nx.compose(canon_sg, updated_sg)

            # Identify bubbles in combined_sg, and make sure that none of them are part of basic transcripts
            bad_introns = non_basic_path[-2]
            bubbles = FindBubbles(combined_sg, canon_path)
            bubbles = [[int(site.split('-')[0]) for site in bubble] for bubble in bubbles]
            bubbles = [bubble for bubble in bubbles if not any([PerfectSublist(bubble, basic_path) for basic_path in basic_paths])]
            bubbles = [bubble for bubble in bubbles if not any([PerfectSublist(list(bad_intron), bubble) for bad_intron in bad_introns])]
            for bubble in bubbles:
                bubble_start, bubble_end = canonTx[4].index(bubble[0]), canonTx[4].index(bubble[-1])
                rescue_path = canonTx[4][:bubble_start] + bubble + canonTx[4][(bubble_end+1):]
                # Make sure rescue_path is not already represented in a non-basic transcript that already passed
                if not any([PerfectSublist(rescue_path[1:-1], check_path) for check_path in non_basic_pass_paths]):
                    rescue_paths[tuple(rescue_path)] = canonTx[-1]

    # Format each set of transcript coordinates in GTF format
    outlines, novel_ctr = [], 0
    for item in pass_group.values:
        if item[5] != '.':
            # This is an annotated transcript
            outlines += FormatGTF(chrom, item[4], strand, item[-1], item[5].split('|')[1])
        else:
            # This is a novel transcript
            novel_ctr += 1
            outlines += FormatGTF(chrom, item[4], strand, item[-1], group_id.replace('XLOC_', 'NovelTx-') + '-' + str(novel_ctr).zfill(6))
    
    for k, v in rescue_paths.items():
        # Each rescue transcript is a novel transcript
        novel_ctr += 1
        outlines += FormatGTF(chrom, list(k), strand, v, group_id.replace('XLOC_', 'NovelTx-') + '-' + str(novel_ctr).zfill(6))

    return outlines

def main():
    message = 'Constructs a high-confidence transcriptome from a merged collection of putative transcripts' 
    parser = argparse.ArgumentParser(description = message)

    # Add arguments
    parser.add_argument('-i', metavar = '/path/to/gffcompare/outprefix', required = True,
        help = 'path to output file prefix used when running GffCompare')
    parser.add_argument('-g', metavar = '/path/to/reference/GTF/file', required = True,
        help = 'path to GTF file with reference transcript annotations from GENCODE')
    parser.add_argument('-f', metavar = '/path/to/reference/genome/FASTA', required = True,
        help = 'path to FASTA file with reference genome sequence')
    parser.add_argument('-x', metavar = '/path/to/refTSS/BED/file', required = True,
        help = 'path to BED file with refTSS peaks')
    parser.add_argument('-y', metavar = '/path/to/polyASite/BED/file', required = True,
        help = 'path to BED file with polyASite peaks')
    parser.add_argument('-o', metavar = '/path/to/output/GTF/file', required = True,
        help = 'path to output GTF file')
    
    # Parse command-line arguments
    args = parser.parse_args()
    inprefix, gencode, genome, reftss, polyasite, outfile = args.i, args.g, args.f, args.x, args.y, args.o

    # Read in gencode and determine which transcripts are basic
    gencodeDF = pd.read_csv(gencode, sep = '\t', header = None, comment = '#')
    gencodeDF = gencodeDF[gencodeDF[2] == 'transcript']
    gencodeDF['gene_id'] = gencodeDF[8].apply(lambda x: PullFeature(x, 'gene_id'))
    gencodeDF['transcript_id'] = gencodeDF[8].apply(lambda x: PullFeature(x, 'transcript_id'))
    gencodeDF['canonical'] = gencodeDF[8].str.contains('tag "Ensembl_canonical"')
    gencodeDF['basic'] = gencodeDF[8].str.contains('tag "basic"')
    basicTx = dict(zip(gencodeDF['gene_id'] + '|' + gencodeDF['transcript_id'], gencodeDF['basic']))
    canonicalTx = dict(zip(gencodeDF['gene_id'] + '|' + gencodeDF['transcript_id'], gencodeDF['canonical']))

    # Read in inprefix.tracking file and determine which consensus transcripts are derived from GENCODE
    trackDF = pd.read_csv(inprefix + '.tracking', sep = '\t', header = None)
    trackDF = trackDF[trackDF[4] != '-'][[0, 1, 4]]
    trackDF[4] = trackDF[4].apply(lambda x: '|'.join(x.split(':')[1].split('|')[:2]))
    trackDict = dict(zip(trackDF[0], trackDF[4]))

    # Read in infile and filter for putative transcripts residing on chr1-22,X,Y,M
    inDF = pd.read_csv(inprefix + '.combined.gtf', sep = '\t', header = None, comment = '#')
    inDF = inDF[inDF[0].isin({'chr' + str(i) for i in range(1, 23)} | {'chrX', 'chrY', 'chrM'}) & (inDF[2] == 'exon')]
    inDF['gene_id'] = inDF[8].apply(lambda x: PullFeature(x, 'gene_id'))
    inDF['transcript_id'] = inDF[8].apply(lambda x: PullFeature(x, 'transcript_id'))

    # Construct exon chains for each transcript and annotate each transcript with their respective annotation status
    inDF = inDF.groupby(['gene_id', 'transcript_id', 0, 6])[[3, 4]].agg(lambda x: sorted(list(x))).reset_index()
    inDF['transcript_chain'] = inDF.agg(lambda x: [site for item in zip(x[3], x[4]) for site in item], axis = 1) 
    inDF = inDF.drop([3, 4], axis = 1)
    inDF['gencode_id'] = inDF['transcript_id'].map(trackDict).fillna('.')
    inDF['gencode_canonical'] = inDF['gencode_id'].map(canonicalTx).fillna(False)
    inDF['gencode_basic'] = inDF['gencode_id'].map(basicTx).fillna(False)

    # Read in genome as a pysam FastaFile object
    genomeFasta = pysam.FastaFile(genome)

    # Read in BED file of refTSS peaks and polyASite peaks
    reftssbed, polyabed = tabix.open(reftss), tabix.open(polyasite)

    # Iterate over rows of inDF
    groupDF, ctr, output = inDF.groupby('gene_id'), 0, []
    for gene_id, group in groupDF:
        ctr += 1

        # Confirm that there exists at least one basic transcript in the group
        if any(group['gencode_basic'].values):
            # Filter group for high-confidence transcripts
            output += FilterTranscripts(group, genomeFasta, reftssbed, polyabed)
        
        if ctr % 10 == 0:
            outlog('Done processing ' + str(ctr) + ' transcript bundles...')
    
    # Save output to outfile
    pd.DataFrame(output).to_csv(outfile, sep = '\t', index = False, header = False, quoting = csv.QUOTE_NONE)
    
if __name__ == '__main__':
    main()

