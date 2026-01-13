#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2026.01.05
Version: 1.4.0

This is a script designed to directly quantify annotated transcripts from alignments of long RNA-seq reads to 
the reference genome. This script requires the following two input files:
    1. An indexed BAM file containing alignments of long RNA-seq reads to the reference genome
    2. A tabix-compressed GTF file containing reference transcript annotations

This script will generate an output file transcript.counts.tsv, which has the following three columns:
    1. Gene ID
    2. Transcript ID
    3. Estimated read counts
'''

# =====================================================================================================================
#                                                       PACKAGES 
# =====================================================================================================================

# Load required packages
from datetime import datetime
import cvxpy as cp
import argparse, os, re, bisect
import numpy as np
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

def GetExons(read):
    '''
    This is a function that retrieves exon-level coordinates for a given read alignment
    '''
    blocks = re.findall('(\d+)(\D+)', read.cigarstring)
    cumlen = np.cumsum([0] + [int(item[0]) if item[1] in {'M', 'D', 'N'} else 0 for item in blocks][:-1])
    introns = sorted([(cumlen[i] + int(read.reference_start) + 1, cumlen[i+1] + int(read.reference_start)) 
        for i in range(len(blocks)) if blocks[i][1] == 'N'])
    
    return [read.reference_start] + [site for intron in introns for site in [intron[0] - 1, intron[1] + 1]] + [read.reference_end]

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

def QuantifyTranscripts(currChains, readAssign):
    '''
    This is a function designed to estimate read counts for transcripts
    '''
    transcripts = [gene + '|' + transcript for gene in currChains.keys() for transcript in currChains[gene][1]]

    if len(readAssign) > 0:
        # Set up read-to-isoform compatibility matrix
        Y = np.array([[int(tuple(transcript.split('|')) in readAssign[read]) for transcript in transcripts] for read in readAssign])

        # Remove rows of Y that are all zeros (non-informative reads)
        Y = Y[~np.all(Y == 0, axis = 1)]

        if Y.shape[0] > 0:
            # Collapse Y into equivalence groups
            P, w = np.unique(Y, axis = 0, return_counts = True)

            # Use splitting cone solver to estimate transcript isoform proportions
            N, K = Y.shape
            theta = cp.Variable(K, nonneg = True)
            nll = cp.sum(cp.multiply(w / N, -cp.log((P @ theta))))
            cons = [cp.sum(theta) == 1]

            prob = cp.Problem(cp.Minimize(nll), cons)
            _ = prob.solve(solver = cp.SCS, max_iters = 10000, eps = 0.0001)

            return [transcripts[idx].split('|') + [theta.value[idx] * N] for idx, value in enumerate(theta.value)]

    return [item.split('|') + [0] for item in transcripts]

def GetCandidateTranscripts(readExons, retainedIntrons, txChains):
    '''
    This is a function designed to retrieve transcripts in txChains that are most compatible with readExons
    '''
    candidateTx = []

    if len(readExons) == 2:
        # Read is mono-exonic
        for tx in txChains:
            if len(txChains[tx]) == 2:
                # Allow a wiggle-room of 20 bp around the ends of the mono-exonic transcript
                # when checking if the read is fully contained
                if txChains[tx][0] - 20 <= readExons[0] and readExons[1] <= txChains[tx][1] + 20:
                    # Check if the mono-exonic transcript contains any retained introns
                    containedIntrons = [intron for intron in retainedIntrons if txChains[tx][0] <= intron[0] and
                        intron[1] <= txChains[tx][1]]
                    
                    if len(containedIntrons) == 0:
                        candidateTx.append((1, tx))
                    else:
                        # Check that all introns in containedIntrons are reliable
                        if all([retainedIntrons[intron] for intron in containedIntrons]):
                            # Finally, make sure that your read actually overlaps all introns in containedIntrons
                            if all([CheckOverlap(intron, readExons) for intron in containedIntrons]):
                                candidateTx.append((1, tx))
                            elif any([CheckOverlap(intron, readExons) for intron in containedIntrons]):
                                candidateTx.append((0, tx))

            else:
                # Permit wiggle-room of 20 bp around ends of the multi-exonic transcript
                updatedChain = [txChains[tx][0] - 20] + txChains[tx][1:-1] + [txChains[tx][-1] + 20]
                if updatedChain[0] <= readExons[0] and readExons[1] <= updatedChain[-1]:
                    bisectIdx = bisect.bisect_right(updatedChain, readExons[0])
                    if bisectIdx % 2 == 1:
                        # Make sure that the current exon fully contains the mono-exonic read
                        if updatedChain[bisectIdx-1] <= readExons[0] and readExons[1] <= updatedChain[bisectIdx]:
                            # Check if the current exon (without adjusted coordinates) contains any retained introns of the gene
                            containedIntrons = [intron for intron in retainedIntrons if txChains[tx][bisectIdx-1] <= intron[0] and 
                                intron[1] <= txChains[tx][bisectIdx]]
                            
                            if len(containedIntrons) > 0:
                                # Check that all introns in containedIntrons are reliable:
                                if all([retainedIntrons[intron] for intron in containedIntrons]):
                                    # Make sure that your read actually overlaps all introns in containedIntrons
                                    if all([CheckOverlap(intron, readExons) for intron in containedIntrons]):
                                        candidateTx.append((0, tx))

    else:
        # Read is multi-exonic
        for tx in txChains:
            # Check how well the terminal exons of the read support the given transcript
            startScore, endScore = -1, -1

            # First check if the internal structure of the read matches the internal structure of the given transcript
            if PerfectSublist(readExons[1:-1], txChains[tx][1:-1]):
                # Locate index of readExons[1] in txChains[tx] and make sure it marks the start of an intron in transcript
                startAnchor = txChains[tx].index(readExons[1])
                if startAnchor % 2 == 1:
                    if startAnchor == 1:
                        # First exon of read matches first exon of transcript
                        # Check if the first exon of transcript contains any introns in retainedIntrons
                        containedIntrons = [intron for intron in retainedIntrons if txChains[tx][startAnchor-1] <= intron[0] and 
                            intron[1] <= txChains[tx][startAnchor]]
                        
                        if len(containedIntrons) == 0:
                            startScore = 1
                        else:
                            # Check if all introns in containedIntrons are reliable:
                            if all([retainedIntrons[intron] for intron in containedIntrons]):
                                # Make sure that the first exon of your read actually overlaps all introns in containedIntrons
                                if all([CheckOverlap(intron, (readExons[0], readExons[1])) for intron in containedIntrons]):
                                    startScore = 1
                                else:
                                    startScore = 0

                    else:
                        # First exon of read matches internal exon of transcript
                        # Start of the read can only be either within the overlapping internal exon or the intron upstream of it
                        if txChains[tx][startAnchor-2] <= readExons[0]:
                            if txChains[tx][startAnchor-1] <= readExons[0]:
                                # Read starts within the overlapping internal exon
                                # Check if this internal exon contains any introns in retainedIntrons
                                containedIntrons = [intron for intron in retainedIntrons if txChains[tx][startAnchor-1] <= intron[0] and 
                                    intron[1] <= txChains[tx][startAnchor]]
                                
                                if len(containedIntrons) == 0:
                                    startScore = 0
                                else:
                                    # Check if all introns in containedIntrons are reliable:
                                    if all([retainedIntrons[intron] for intron in containedIntrons]):
                                        # Make sure that the first exon of your read actually overlaps all introns in containedIntrons
                                        if all([CheckOverlap(intron, (readExons[0], readExons[1])) for intron in containedIntrons]):
                                            startScore = 0
                                        else:
                                            startScore = 0

                            else:
                                # Read starts within the intron upstream of the overlapping internal exon
                                # We can still use this read to support the transcript as long as it does not 
                                # bleed into the upstream intron by more than 6 bp
                                if txChains[tx][startAnchor-1] - readExons[0] <= 6:
                                    startScore = 0
                    
                # Locate index of readExons[-2] in txChains[tx] and make sure it marks the start of an exon in transcript
                endAnchor = txChains[tx].index(readExons[-2])
                if endAnchor % 2 == 0:
                    if endAnchor == len(txChains[tx]) - 2:
                        # Last exon of read matches last exon of transcript
                        # Check if the last exon of transcript contains any introns in retainedIntrons
                        containedIntrons = [intron for intron in retainedIntrons if txChains[tx][endAnchor] <= intron[0] and 
                            intron[1] <= txChains[tx][endAnchor+1]]
                        
                        if len(containedIntrons) == 0:
                            endScore = 1
                        else:
                            # Check if all introns in containedIntrons are reliable:
                            if all([retainedIntrons[intron] for intron in containedIntrons]):
                                # Make sure that the last exon of your read actually overlaps all introns in containedIntrons
                                if all([CheckOverlap(intron, (readExons[-2], readExons[-1])) for intron in containedIntrons]):
                                    endScore = 1
                                else:
                                    endScore = 0
                    
                    else:
                        # Last exon of read matches internal exon of transcript
                        # End of the read can only be either within the overlapping internal exon or the intron downstream of it
                        if txChains[tx][endAnchor+2] >= readExons[-1]:
                            if txChains[tx][endAnchor+1] >= readExons[-1]:
                                # Read ends within the overlapping internal exon
                                # Check if this internal exon contains any introns in retainedIntrons
                                containedIntrons = [intron for intron in retainedIntrons if txChains[tx][endAnchor] <= intron[0] and 
                                    intron[1] <= txChains[tx][endAnchor+1]]
                                
                                if len(containedIntrons) == 0:
                                    endScore = 0
                                else:
                                    # Check if all introns in containedIntrons are reliable:
                                    if all([retainedIntrons[intron] for intron in containedIntrons]):
                                        # Make sure that the last exon of your read actually overlaps all introns in containedIntrons
                                        if all([CheckOverlap(intron, (readExons[-2], readExons[-1])) for intron in containedIntrons]):
                                            endScore = 0
                                        else:
                                            endScore = 0
                            
                            else:
                                # Read ends within the intron downstream of the overlapping internal exon
                                # We can still use this read to support the transcript as long as it does not
                                # bleed into the downstream intron by more than 6 bp
                                if readExons[-1] - txChains[tx][endAnchor+1] <= 6:
                                    endScore = 0

                if min(startScore, endScore) >= 0:
                    candidateTx.append((min(startScore, endScore), tx))

    if len(candidateTx) > 0:
        # Check if we have any transcripts with a minimum score of 1
        if any([item[0] == 1 for item in candidateTx]):
            return [item[1] for item in candidateTx if item[0] == 1]
        else:
            return [item[1] for item in candidateTx]
    
    else:
        return []

def CheckReads(samfile, region, currChains, blacklist):
    '''
    This is a function designed to pull out reads mapping to a specific genomic region and check their
    compatibility with transcripts represented in currChains
    '''
    readAssign = dict()

    # Iterate over primary/supplementary alignments of any mapping quality in region
    for read in samfile.fetch(region[0], region[1], region[2]):
        if not read.is_secondary:
            # Get exon-level coordinates for read alignment
            readExons = GetExons(read)
            compatibleTx = set()

            for gene in currChains:
                # Pull out retained introns for gene as well as corresponding transcript chains
                retainedIntrons, txChains = currChains[gene]

                # Retrieve candidate transcripts for the read and update compatibleTx
                candidateTx = GetCandidateTranscripts(readExons, retainedIntrons, txChains)
                compatibleTx = compatibleTx | {(gene, tx) for tx in candidateTx if tx not in blacklist}

            # Update readAssign with compatibleTx
            if read.query_name not in readAssign:
                readAssign[read.query_name] = compatibleTx
            else:
                readAssign[read.query_name] = readAssign[read.query_name] & compatibleTx

    return readAssign

def CheckRetainment(transcript, intron):
    '''
    This is a function designed to check whether a given intron is retained in a transcript
    '''
    return any([exon[0] < intron[0] and intron[1] < exon[1] for exon in list(zip(transcript[0::2], transcript[1::2]))])

def BuildChains(tbfile, region):
    '''
    This is a function designed to construct genomic chains for transcripts of each gene contained within a specific region
    This function also reports introns that could potentially be retained
    '''
    chainDF = pd.DataFrame([record for record in tbfile.query(region[0], region[1], region[2]) if record[2] == 'exon'])
    chainDF[3], chainDF[4] = chainDF[3].astype(int), chainDF[4].astype(int)
    chainDF['gene_id'] = chainDF[8].apply(lambda x: PullFeature(x, 'gene_id'))
    chainDF['transcript_id'] = chainDF[8].apply(lambda x: PullFeature(x, 'transcript_id'))
    chainDF = chainDF.groupby(['gene_id', 'transcript_id'])[[3, 4]].agg(lambda x: sorted(list(x))).reset_index()
    chainDF['transcript_chain'] = chainDF.apply(lambda x: [coord for item in zip(x[3], x[4]) for coord in item], axis = 1)

    chainDict = dict()
    genes = chainDF['gene_id'].unique()
    for gene in genes:
        geneChains = chainDF[chainDF['gene_id'] == gene].copy()

        # Identify unique introns in geneChains
        geneIntrons = sorted(list(set(geneChains['transcript_chain'].apply(lambda x: list(zip(x[1:-1][0::2], x[1:-1][1::2]))).sum())))

        # Filter geneIntrons for introns that could potentially be retained based on reported transcript structures
        # We will also indicate whether the intron is retained in at least one annotated transcript or not
        retainedIntrons = [intron for intron in geneIntrons if any([CheckRetainment(transcript, intron) for transcript in geneChains['transcript_chain']])]
        annoRetainedIntrons = [intron for intron in geneIntrons if any([CheckRetainment(transcript, intron) for transcript in geneChains[~geneChains['transcript_id'].str.contains('NovelTx')]['transcript_chain']])]
        retainedIntrons = [(intron, intron in annoRetainedIntrons) for intron in retainedIntrons]

        chainDict[gene] = (retainedIntrons, dict(zip(geneChains['transcript_id'], geneChains['transcript_chain'])))

    return chainDict

def AssessRetainedIntrons(samfile, region, currChains):
    '''
    This is a function designed to evaluate potential retained introns based on whether they actually show strong evidence of retention
    in raw sequencing data. Specifically, if an intron is known to be retained (relative to an annotated transcript), then by default
    we will believe that it is retained (regardless of what the sequencing data may suggest). Alternatively, if an intron is not known
    to be retained, then we will do the following:
        1. Compute read coverage over the intron and partition this vector into bins of 20
        2. Compute the ratio between the bin with the smallest coverage and the larger of the two flanking intronic bins
        3. If this ratio is less than 75%, then we will assume that this intron is not reliably retained
    '''
    outDict = dict()

    for gene in currChains:
        results, retainedIntrons = dict(), currChains[gene][0]
        
        # Pull out array of read coverage for each retained intron
        for ri in retainedIntrons:
            if ri[1]:
                results[ri[0]] = True
            else:
                # Make sure coverage over entire intron is uniform
                coverageInfo = np.array(samfile.count_coverage(region[0], ri[0][0], ri[0][1]-1, quality_threshold = 0)).sum(axis = 0)
                coverageBins = np.array([np.mean(bin) for bin in np.array_split(coverageInfo, 20)]) if len(coverageInfo) > 20 else coverageInfo 
                results[ri[0]] = np.min(coverageBins)/(max(coverageBins[0], coverageBins[-1]) + 1) >= 0.75

        outDict[gene] = (results, currChains[gene][1])
    
    return outDict

def CheckNovelInternalExons(samfile, region, currChains):
    '''
    This is a function designed to evaluate the reliability of novel internal exons. Specifically, we want to check if these internal exons
    have adequate read coverage (to eliminate spurious novel exons that are super long and sponge up read support during transcript
    quantification). Specifically, we will do the following:
        1. Compute read coverage over the exon and partition this vector into bins of 20
        2. Compute the ratio between the bin with the smallest coverage and the smaller of the two flanking exonic bins
        3. If this ratio is less than 50%, then we will assume that this exon is not reliably used
    
    Furthermore, we will require this exon to show up in at least two reads
    '''
    blacklist, internalExonCounter = set(), dict()

    # Keep track of all internal exons for reads mapping to region
    for read in samfile.fetch(region[0], region[1], region[2]):
        if not read.is_secondary:
            # Get exon-level coordinates for read alignment
            readExons = GetExons(read)
            if len(readExons) > 4:
                for exon in list(zip(readExons[2:-2][0::2], readExons[2:-2][1::2])):
                    internalExonCounter[exon] = internalExonCounter.get(exon, 0) + 1

    for gene in currChains:
        badNovelInternalExons = set()

        # Gather all novel internal exons
        annoInternalExons = list({exon for k, v in currChains[gene][1].items() if 'NovelTx' not in k and len(v) > 4 for exon in zip(v[2:-2][0::2], v[2:-2][1::2])})
        novelInternalExons = list({exon for k, v in currChains[gene][1].items() if 'NovelTx' in k and len(v) > 4 for exon in 
            zip(v[2:-2][0::2], v[2:-2][1::2]) if exon not in annoInternalExons})

        # Pull out array of read coverage for each novel internal exon
        for exon in novelInternalExons:
            coverageInfo = np.array(samfile.count_coverage(region[0], exon[0]-1, exon[1], quality_threshold = 0)).sum(axis = 0)
            coverageBins = np.array([np.mean(bin) for bin in np.array_split(coverageInfo, 20)]) if len(coverageInfo) > 20 else coverageInfo 
            if np.min(coverageBins)/(min(coverageBins[0], coverageBins[-1])+1) < 0.5:
                badNovelInternalExons.add(exon)
            if internalExonCounter.get(exon, 0) < 2:
                badNovelInternalExons.add(exon)

        # Figure out which novel transcripts harbor at least one exon in badNovelInternalExons
        blacklist = blacklist | {k for k, v in currChains[gene][1].items() if any([exon in badNovelInternalExons for exon in zip(v[0::2], v[1::2])])}
    
    return blacklist
            
def SplitAnnotations(gtffile):
    '''
    This is a function designed to split reference transcript annotations into non-overlapping regions
    '''
    # Read in gtffile as a pandas dataframe and filter for transcript-level annotations
    transcripts = pd.read_csv(gtffile, sep = '\t', header = None, compression = 'gzip')
    transcripts = transcripts[transcripts[2] == 'transcript'].copy()
    transcripts['gene_id'] = transcripts[8].apply(lambda x: PullFeature(x, 'gene_id'))
    transcripts = transcripts.groupby(['gene_id', 0])[[3, 4]].agg(lambda x: sorted(list(x))).reset_index()
    transcripts[3], transcripts[4] = transcripts[3].str[0], transcripts[4].str[-1]
    transcripts = transcripts.drop('gene_id', axis = 1).sort_values(by = [0, 3, 4]).reset_index(drop = True)

    # Collapse genes into non-overlapping genomic regions
    regions, cluster, genes = [], None, transcripts.values

    for idx in range(len(genes)):
        if cluster is None:
            cluster = list(genes[idx])
        else:
            if genes[idx][0] != cluster[0]:
                regions.append(cluster)
                cluster = list(genes[idx])
            else:
                if CheckOverlap((genes[idx][1], genes[idx][2]), (cluster[1], cluster[2])):
                    cluster = [cluster[0], min(cluster[1], genes[idx][1]), max(cluster[2], genes[idx][2])]
                else:
                    regions.append(cluster)
                    cluster = list(genes[idx])
    regions.append(cluster)

    return regions
        
def main():
    message = 'Quantifies annotated transcripts from long RNA-seq read alignments to the genome'
    parser = argparse.ArgumentParser(description = message)

    # Add arguments
    parser.add_argument('-i', metavar = '/path/to/input/BAM/file', required = True,
        help = 'path to BAM file with long RNA-seq read alignments to genome')
    parser.add_argument('-g', metavar = '/path/to/reference/GTF/file', required = True,
        help = 'path to GTF file with reference transcript annotations')
    parser.add_argument('-o', metavar = '/path/to/output/file', required = True,
        help = 'path to output file')
    
    # Parse command-line arguments
    args = parser.parse_args()
    infile, gtffile, outfile = args.i, args.g, args.o

    # Split transcript annotations into non-overlapping genomic regions
    outlog('Splitting reference transcript annotations into non-overlapping regions...')
    regions = SplitAnnotations(gtffile)

    # Iterate over regions and quantify transcripts in each region
    tbfile = tabix.open(gtffile)
    samfile = pysam.AlignmentFile(infile, 'rb')

    outlog('Starting transcript quantification...')
    output, errorlog = [], []
    for idx, region in enumerate(regions):
        # Construct genomic "chains" for all transcripts of each gene within region
        currChains = BuildChains(tbfile, region)

        # For each gene in currChains, pull out the set of candidate retained introns and determine
        # which ones look reliable for downstream quantification
        currChains = AssessRetainedIntrons(samfile, region, currChains)

        # Blacklist transcripts with novel internal exons that have very poor read support
        # We will assign their corresponding transcripts a read count of 0 (because they should not exist)
        blacklist = CheckNovelInternalExons(samfile, region, currChains)

        # For reads mapping to region, determine their compatibility with transcripts in currChains
        readAssign = CheckReads(samfile, region, currChains, blacklist)
                
        # Estimate read counts for each transcript in currChains
        output += QuantifyTranscripts(currChains, readAssign)

        if idx % 10 == 0 and idx > 0:
            outlog('Done processing ' + str(idx) + ' regions...')
    
    # Save output to outfile
    outDF = pd.DataFrame(output, columns = ['gene_id', 'transcript_id', 'count'])
    outDF = outDF.sort_values(by = ['gene_id', 'transcript_id'])
    outDF.to_csv(outfile, sep = '\t', index = False, float_format = '%.2f')

if __name__ == '__main__':
    main()

