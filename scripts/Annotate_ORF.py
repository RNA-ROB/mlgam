#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2026.01.05
Version: 2.3.0

This is a script designed to update a GTF file produced by a transcript assembly tool to include open reading frame
annotations (e.g., CDS, UTR, start/stop codon coordinates) for each transcript where applicable. This script requires
the following three files as input:
    * A GTF file describing genomic features of discovered transcripts
    * A GTF file of reference gene annotations from GENCODE
    * A FASTA file containing the reference genome sequence

We will use the following approach to identify open reading frames (at least 30 amino acids) for each transcript:
* If there exists at least one annotated start codon that overlaps the given transcript sequence:
  * Go with the longest open reading frame starting from an annotated start codon
* If no annotated start codon overlaps the transcript sequence:
  * Go with the longest open reading frame from any ATG in the transcript sequence
  * The % identity between the putative ORF and the canonical ORF must be at least 60%
'''

# =====================================================================================================================
#                                                       PACKAGES 
# =====================================================================================================================

# Load required packages
import argparse, re, csv
from Bio import Align
from datetime import datetime
import numpy as np
import pandas as pd
import pysam

# =====================================================================================================================
#                                                   GLOBAL VARIABLES
# =====================================================================================================================

# Establish the following global variables:
#   * revcompDict: dictionary mapping a nucleotide with its complement
#   * translateDict: dictionary mapping a codon triplet with its amino acid

revcompDict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
translateDict = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L',
    'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G',
    'GGA': 'G', 'GGG': 'G'}

# Initialize a PairwiseAligner (following what BLASTP uses)
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = Align.substitution_matrices.load('BLOSUM62')
aligner.open_gap_score = -11
aligner.extend_gap_score = -1
aligner.mode = 'local'

# =====================================================================================================================
#                                                   HELPER FUNCTIONS
# =====================================================================================================================

def outlog(myString):
    '''
    This is a function to print some output message string (myString) with the date and time
    '''
    print('[', datetime.now().strftime("%Y-%m-%d %H:%M:%S"), '] ', myString, sep = '', flush = True)

def revcomp(seq):
    '''
    This is a function that will generate the reverse complement of an input nucleotide sequence
    '''
    return ''.join([revcompDict.get(nt, 'N') for nt in seq[::-1]])

def clipseq(seq, n):
    '''
    This is a function that will clip a sequence to have a length that is a multiple of n
    '''
    return seq[0:n*int(len(seq)/n)]

def genome2TxCoord(inputCoord, exonCoord, strand):
    '''
    This is a function that will convert input genome-level coordinates to transcript-level coordinates
    relative to exonCoord (returns -1 if inputCoord does not overlap transcript)
    '''
    # Locate which exon harbors the input coordinate
    inputIdx = next((i for i, v in enumerate(exonCoord) if v[0] <= inputCoord and v[1] >= inputCoord), None)

    if inputIdx is not None:
        if strand == '+':
            # Compute intron lengths (with respect to transcript strand)
            intronLenCum = np.cumsum([0] + [exonCoord[i+1][0] - exonCoord[i][1] - 1 for i in range(len(exonCoord)-1)])
            offset = exonCoord[0][0]

            # Compute transcript-level coordinates for input coordinate
            txCoord = inputCoord - offset - intronLenCum[inputIdx]

        else:
            # Compute intron lengths (with respect to transcript strand)
            intronLenCum = np.cumsum([0] + [exonCoord[i+1][0] - exonCoord[i][1] - 1 for i in range(len(exonCoord)-1)][::-1])[::-1]
            offset = exonCoord[-1][1]

            # Compute transcript-level coordinates for input coordinate
            txCoord = offset - intronLenCum[inputIdx] - inputCoord

        return int(txCoord)
        
    else:
        # Input coordinate does not overlap transcript
        return -1

def tx2GenomeCoord(inputCoord, exonCoord, strand):
    '''
    This is a function that will convert input transcript-level coordinates to genome-level coordinates
    relative to exonCoord
    '''
    if strand == '+':
        # Compute intron lengths (with respect to transcript strand)
        intronLenCum = np.cumsum([0] + [exonCoord[i+1][0] - exonCoord[i][1] - 1 for i in range(len(exonCoord)-1)])
        offset = exonCoord[0][0]

        # Compute transcript-level coordinates for exons 
        adjExonCoord = [(exonCoord[i][0] - offset - intronLenCum[i], exonCoord[i][1] - offset - intronLenCum[i]) for i in range(len(exonCoord))]

        # Locate which exon harbors the input transcript-level coordinate
        inputIdx = next(i for i, v in enumerate(adjExonCoord) if v[0] <= inputCoord and v[1] >= inputCoord)

        # Compute genome-level coordinates for input coordinate
        genCoord = inputCoord + offset + intronLenCum[inputIdx]

    else:
        # Compute intron lengths (with respect to transcript strand)
        intronLenCum = np.cumsum([0] + [exonCoord[i+1][0] - exonCoord[i][1] - 1 for i in range(len(exonCoord)-1)][::-1])[::-1]
        offset = exonCoord[-1][1]

        # Compute transcript-level coordinates for exons 
        adjExonCoord = [(offset - exonCoord[i][0] - intronLenCum[i], offset - exonCoord[i][1] - intronLenCum[i]) for i in range(len(exonCoord))]

        # Locate which exon harbors the input transcript-level coordinate
        inputIdx = next(i for i, v in enumerate(adjExonCoord) if v[1] <= inputCoord and v[0] >= inputCoord)

        # Compute genome-level coordinates for input coordinate
        genCoord = offset - intronLenCum[inputIdx] - inputCoord
    
    return int(genCoord)

# =====================================================================================================================
#                                                    MAIN FUNCTIONS
# =====================================================================================================================

def getAnnoStartCodon(annoDF):
    '''
    This is a function to extract genomic start coordinates of annotated start codons in annoDF
    '''
    # Filter annoDF for entries corresponding to start codons
    annoStartCodon = annoDF[annoDF[2] == 'start_codon'][[0,3,4,6,8]]
    annoStartCodon.columns = ['chrom', 'startCodonStart', 'startCodonEnd', 'strand', 'txInfo']

    # Extract geneID and txID from txInfo column, then group start codons by geneID, txID, chrom, and strand
    # Also check if the start codon is derived from a basic transcript or not
    annoStartCodon['geneID'] = annoStartCodon['txInfo'].apply(lambda x: [item.split('"')[1] for item in x.split(';') if 'gene_id' in item][0])
    annoStartCodon['txID'] = annoStartCodon['txInfo'].apply(lambda x: [item.split('"')[1] for item in x.split(';') if 'transcript_id' in item][0])
    annoStartCodon['basic'] = annoStartCodon['txInfo'].str.contains('tag "basic";')
    annoStartCodon = annoStartCodon.groupby(['geneID', 'txID', 'basic', 'chrom', 'strand'])[['startCodonStart', 'startCodonEnd']].agg(lambda x: sorted(list(x))).reset_index()

    # Make sure that start codons spanning splice junctions are treated as individual entries in annoStartCodon
    annoStartCodon['startCodonStart'], annoStartCodon['startCodonEnd'] = annoStartCodon['startCodonStart'].apply(min), annoStartCodon['startCodonEnd'].apply(max)

    # Identify the genomic start coordinate of each start codon with respect to transcript strand
    annoStartCodon['startCoord'] = annoStartCodon.apply(lambda x: x['startCodonStart'] if x['strand'] == '+' else x['startCodonEnd'], axis = 1)
    annoStartCodon = annoStartCodon.groupby(['geneID', 'chrom', 'startCoord'])['basic'].agg(any).reset_index()
   
    return dict(zip(annoStartCodon['geneID'] + '_' + annoStartCodon['chrom'] + '_' + annoStartCodon['startCoord'].astype(str), annoStartCodon['basic']))

def getCanonProteinSeq(annoDF, genome):
    '''
    This is a function to extract canonical protein sequences for all protein-coding genes in annoDF
    '''
    # Filter annoDF for CDS entries derived from canonical protein-coding transcripts
    cdsCanon = annoDF[(annoDF[2] == 'CDS') & annoDF[8].str.contains('Ensembl_canonical') & annoDF[8].str.contains('transcript_type "protein_coding"')].copy()
    cdsCanon['geneID'] = cdsCanon[8].apply(lambda x: [item.split('"')[1] for item in x.split(';') if 'gene_id' in item][0])
    cdsCanon = cdsCanon.groupby(['geneID', 0, 6])[[3, 4]].agg(lambda x: sorted(list(x))).reset_index()

    # Retrieve protein sequences for genes in cdsCanon
    def getProteinSeq(x, genome):
        genomeSeq = [genome.fetch(x[0], x[3][i] - 1, x[4][i]).upper() for i in range(len(x[3]))]
        txSequence = revcomp(''.join(genomeSeq)) if x[6] == '-' else ''.join(genomeSeq)
        return ''.join([translateDict.get(txSequence[i:i+3], 'X') for i in range(0, len(txSequence), 3)])

    cdsCanon['protein_seq'] = cdsCanon.apply(lambda x: getProteinSeq(x, genome), axis = 1)
    return dict(zip(cdsCanon['geneID'], cdsCanon['protein_seq']))

def getLongestORF(transcript, annoStartCodon, canonProteinSeq, genome):
    '''
    This is a function that will retrieve the longest open reading frame for a given transcript
    starting from any ATG in the transcript sequence (with priority placed on annotated start codons)

    *Update (as of 1/5/2026): Higher priority will be placed on start codons from GENCODE basic transcripts!
    '''
    # Retrieve genome and transcript sequences
    genomeSeq = [genome.fetch(transcript['chrom'], transcript['exonCoord'][i][0] - 1, transcript['exonCoord'][i][1]).upper() for i in range(transcript['numExon'])]
    txSequence = revcomp(''.join(genomeSeq)) if transcript['strand'] == '-' else ''.join(genomeSeq)
    
    # Retrieve putative ORFs from all three possible open reading frames
    candidates = []
    for shift in range(3):
        frameSeq = clipseq(txSequence[shift:], 3)
        frameSeq = ''.join([translateDict.get(frameSeq[i:i+3], 'X') for i in range(0, len(frameSeq), 3)])

        # Locate all M's and *'s in frameSeq
        myIdx = [(i, l) for i, l in enumerate(frameSeq) if l in {'M', '*'}]

        # Flag which M's in myIdx are annotated or not
        myIdx = [(tup[0], tup[1], transcript['geneID'] + '_' + transcript['chrom'] + '_' + str(tx2GenomeCoord(tup[0]*3 + shift, transcript['exonCoord'], transcript['strand'])) 
            in annoStartCodon) if tup[1] == 'M' else (tup[0], tup[1], False) for tup in myIdx]
        myIdxAnno = [tup for tup in myIdx if (tup[2] == True) or (tup[1] == '*')]
        myIdxOther = [tup for tup in myIdx if tup[2] == False]

        # Pass 1. Attempt translation from overlapping annotated start codons
        findM, currSeg, segments = True, [], []
        for tup in myIdxAnno:
            if findM and tup[1] == 'M':
                # Found an M to define the beginning of a putative ORF (i.e., "segment")
                currSeg.append(tup[0])
                findM = False
            elif not findM and tup[1] == '*':
                # Segment is now completed
                currSeg.append(tup[0])

                # Transform items in currSeg into genome coordinates
                transSeg = [tx2GenomeCoord(item*3 + shift, transcript['exonCoord'], transcript['strand']) for item in currSeg]
                segments.append(transSeg + [currSeg[1] - currSeg[0], True, frameSeq[currSeg[0]:currSeg[1]]])
                currSeg, findM = [], True

        # If currSeg is non-empty, then our transcript does not harbor a putative stop codon
        if len(currSeg) == 1:
            transSeg = [tx2GenomeCoord(currSeg[0]*3 + shift, transcript['exonCoord'], transcript['strand']), 'NA']
            segments.append(transSeg + [len(frameSeq)-currSeg[0]+1, True, frameSeq[currSeg[0]:]])

        # Add segments to candidates
        candidates += segments

        # Pass 2. Attempt translation from any other putative start codon
        findM, currSeg, segments = True, [], []
        for tup in myIdxOther:
            if findM and tup[1] == 'M':
                # Found an M to define the beginning of a putative ORF (i.e., "segment")
                currSeg.append(tup[0])
                findM = False
            elif not findM and tup[1] == '*':
                # Segment is now completed
                currSeg.append(tup[0])
                # Transform items in currSeg into genome coordinates
                transSeg = [tx2GenomeCoord(item*3 + shift, transcript['exonCoord'], transcript['strand']) for item in currSeg]
                segments.append(transSeg + [currSeg[1] - currSeg[0], False, frameSeq[currSeg[0]:currSeg[1]]])
                currSeg, findM = [], True

        # If currSeg is non-empty, then our transcript does not harbor a putative stop codon
        if len(currSeg) == 1:
            transSeg = [tx2GenomeCoord(currSeg[0]*3 + shift, transcript['exonCoord'], transcript['strand']), 'NA']
            segments.append(transSeg + [len(frameSeq)-currSeg[0]+1, False, frameSeq[currSeg[0]:]])

        # Add segments to candidates
        candidates += segments
        
    # Transform candidates into a dataframe and filter for ORFs that are at least 30 amino acids long
    currDF = pd.DataFrame(candidates, columns=[0, 1, 2, 3, 4])
    currDF = currDF[currDF[2] >= 30]

    if currDF.shape[0] > 0:
        # Compute percentage identity between ORFs in currDF and the canonical ORF and flag ORFs with at least 60% identity
        if transcript['geneID'] in canonProteinSeq:
            currDF[5] = currDF.apply(lambda x: format(aligner.align(canonProteinSeq[transcript['geneID']], x[4])[0]).split('\n')[1].count('|')/len(x[4]) >= 0.6, axis = 1)
        else:
            currDF[5] = False
        
        if any(currDF[3]):
            # Add sixth column indicating whether the start codon is derived from a basic transcript
            currDF = currDF[currDF[3]].copy()
            currDF[6] = [annoStartCodon[item] for item in transcript['geneID'] + '_' + transcript['chrom'] + '_' + currDF[0].astype(str)]
            return tuple(currDF.sort_values(by = [6, 2], ascending = False).iloc[0, 0:2])
        else:
            if any(currDF[5]):
                return tuple(currDF[currDF[5]].sort_values(by = 2, ascending = False).iloc[0, 0:2])
            else:
                return None
    else:
        return None

def injectCDS(transcript, breakpoints):
    '''
    This is a function to update transcript annotations to include CDS features
    '''
    # Determine exon numbering from transcript strand
    exonNumArr = list(range(1, transcript['numExon']+1)) if transcript['strand'] == '+' else list(range(transcript['numExon'], 0, -1))

    # Transform transcript into a dataframe and add in exon number
    myDF = transcript.to_frame().T.explode('exonCoord').drop('numExon', axis = 1)
    myDF['exonNum'] = exonNumArr

    # Figure out where breakpoints should be inserted in myDF
    cdsStartIdx = next(i for i, v in enumerate(myDF['exonCoord']) if breakpoints[0] >= v[0] and breakpoints[0] <= v[1])
    cdsEndIdx = next(i for i, v in enumerate(myDF['exonCoord']) if breakpoints[1] >= v[0] and breakpoints[1] <= v[1]) if breakpoints[1] != 'NA' else None

    # Assemble 5' UTR annotations
    utr5DF = myDF.iloc[:cdsStartIdx].copy() if transcript['strand'] == '+' else myDF.iloc[(cdsStartIdx+1):].copy()
    if transcript['strand'] == '+' and myDF.iloc[cdsStartIdx]['exonCoord'][0] < breakpoints[0]:
        utr5Append = myDF.iloc[cdsStartIdx].copy()
        utr5Append['exonCoord'] = (utr5Append['exonCoord'][0], breakpoints[0] - 1)
        utr5DF = pd.concat([utr5DF, utr5Append.to_frame().T])
    elif transcript['strand'] == '-' and myDF.iloc[cdsStartIdx]['exonCoord'][1] > breakpoints[0]:
        utr5Append = myDF.iloc[cdsStartIdx].copy()
        utr5Append['exonCoord'] = (breakpoints[0] + 1, utr5Append['exonCoord'][1])
        utr5DF = pd.concat([utr5Append.to_frame().T, utr5DF])
    utr5DF['biotype'], utr5DF['phase'] = 'UTR', '.'

    # Assemble start codon annotations
    startCodonSeries = myDF.iloc[cdsStartIdx].copy()
    startCodonEndCoord = tx2GenomeCoord(genome2TxCoord(breakpoints[0], transcript['exonCoord'], transcript['strand']) + 2, transcript['exonCoord'], transcript['strand'])
    if abs(startCodonEndCoord - breakpoints[0]) == 2:
        startCodonSeries['exonCoord'] = (breakpoints[0], startCodonEndCoord) if transcript['strand'] == '+' else (startCodonEndCoord, breakpoints[0])
        startCodonDF = startCodonSeries.to_frame().T
        startCodonDF['biotype'], startCodonDF['phase'] = 'start_codon', 0
    else:
        # Start codon is separated by a splice junction
        startCodonSeries1, startCodonSeries2 = startCodonSeries.copy(), startCodonSeries.copy()
        startCodonMidCoord = tx2GenomeCoord(genome2TxCoord(breakpoints[0], transcript['exonCoord'], transcript['strand']) + 1, transcript['exonCoord'], transcript['strand'])
        if abs(startCodonMidCoord - breakpoints[0]) == 1:
            startCodonSeries1['exonCoord'] = (breakpoints[0], startCodonMidCoord) if transcript['strand'] == '+' else (startCodonMidCoord, breakpoints[0])
            startCodonSeries2['exonCoord'] = (startCodonEndCoord, startCodonEndCoord)
            startCodonDF1, startCodonDF2 = startCodonSeries1.to_frame().T, startCodonSeries2.to_frame().T
            startCodonDF1['biotype'], startCodonDF1['phase'] = 'start_codon', 0
            startCodonDF2['biotype'], startCodonDF2['phase'] = 'start_codon', 1
        else:
            startCodonSeries1['exonCoord'] = (breakpoints[0], breakpoints[0])
            startCodonSeries2['exonCoord'] = (startCodonMidCoord, startCodonEndCoord) if transcript['strand'] == '+' else (startCodonEndCoord, startCodonMidCoord)
            startCodonDF1, startCodonDF2 = startCodonSeries1.to_frame().T, startCodonSeries2.to_frame().T
            startCodonDF1['biotype'], startCodonDF1['phase'] = 'start_codon', 0
            startCodonDF2['biotype'], startCodonDF2['phase'] = 'start_codon', 2
        startCodonDF = pd.concat([startCodonDF1, startCodonDF2])
    
    # Assemble CDS annotations
    if cdsEndIdx is not None:
        if cdsStartIdx != cdsEndIdx:
            cdsDF = myDF.iloc[(cdsStartIdx+1):cdsEndIdx].copy() if transcript['strand'] == '+' else myDF.iloc[(cdsEndIdx+1):cdsStartIdx].copy()
            cdsStartAppend = myDF.iloc[cdsStartIdx].copy()
            cdsStartAppend['exonCoord'] = (breakpoints[0], cdsStartAppend['exonCoord'][1]) if transcript['strand'] == '+' else (cdsStartAppend['exonCoord'][0], breakpoints[0])
            cdsDF = pd.concat([cdsStartAppend.to_frame().T, cdsDF]) if transcript['strand'] == '+' else pd.concat([cdsDF, cdsStartAppend.to_frame().T])
            if transcript['strand'] == '+' and myDF.iloc[cdsEndIdx]['exonCoord'][0] < breakpoints[1]:
                cdsEndAppend = myDF.iloc[cdsEndIdx].copy()
                cdsEndAppend['exonCoord'] = (cdsEndAppend['exonCoord'][0], breakpoints[1] - 1)
                cdsDF = pd.concat([cdsDF, cdsEndAppend.to_frame().T])
            elif transcript['strand'] == '-' and myDF.iloc[cdsEndIdx]['exonCoord'][1] > breakpoints[1]:
                cdsEndAppend = myDF.iloc[cdsEndIdx].copy()
                cdsEndAppend['exonCoord'] = (breakpoints[1] + 1, cdsEndAppend['exonCoord'][1])
                cdsDF = pd.concat([cdsEndAppend.to_frame().T, cdsDF])
        else:
            # ORF is found entirely within a single exon
            cdsSeries = myDF.iloc[cdsStartIdx].copy()
            cdsSeries['exonCoord'] = (breakpoints[0], breakpoints[1] - 1) if transcript['strand'] == '+' else (breakpoints[1] + 1, breakpoints[0])
            cdsDF = cdsSeries.to_frame().T

    else:
        cdsDF = myDF.iloc[(cdsStartIdx+1):].copy() if transcript['strand'] == '+' else myDF.iloc[:cdsStartIdx].copy()
        cdsStartAppend = myDF.iloc[cdsStartIdx].copy()
        cdsStartAppend['exonCoord'] = (breakpoints[0], cdsStartAppend['exonCoord'][1]) if transcript['strand'] == '+' else (cdsStartAppend['exonCoord'][0], breakpoints[0])
        cdsDF = pd.concat([cdsStartAppend.to_frame().T, cdsDF]) if transcript['strand'] == '+' else pd.concat([cdsDF, cdsStartAppend.to_frame().T])
    cdsDF['biotype'], cdsDF['featureLen'] = 'CDS', cdsDF['exonCoord'].apply(lambda x: x[1] - x[0] + 1)
    cdsDF = cdsDF.sort_values(by = 'exonCoord', ascending = True) # Make sure rows of cdsDF are sorted by coordinate
    cdsDF['featureCumLen'] = [0] + np.cumsum(cdsDF['featureLen']).to_list()[:-1] if transcript['strand'] == '+' else ([0] + np.cumsum(cdsDF['featureLen'][::-1]).to_list()[:-1])[::-1]
    cdsDF['phase'] = cdsDF['featureCumLen'].apply(lambda x: 0 if x % 3 == 0 else 3 - (x % 3))
    cdsDF = cdsDF.drop(['featureLen', 'featureCumLen'], axis=1)

    # Assemble stop codon annotations
    if cdsEndIdx is not None:
        stopCodonSeries = myDF.iloc[cdsEndIdx].copy()
        stopCodonEndCoord = tx2GenomeCoord(genome2TxCoord(breakpoints[1], transcript['exonCoord'], transcript['strand']) + 2, transcript['exonCoord'], transcript['strand'])
        if abs(stopCodonEndCoord - breakpoints[1]) == 2:
            stopCodonSeries['exonCoord'] = (breakpoints[1], stopCodonEndCoord) if transcript['strand'] == '+' else (stopCodonEndCoord, breakpoints[1])
            stopCodonDF = stopCodonSeries.to_frame().T
            stopCodonDF['biotype'], stopCodonDF['phase'] = 'stop_codon', 0
        else:
            # Stop codon is separated by a splice junction
            stopCodonSeries1, stopCodonSeries2 = stopCodonSeries.copy(), stopCodonSeries.copy()
            stopCodonMidCoord = tx2GenomeCoord(genome2TxCoord(breakpoints[1], transcript['exonCoord'], transcript['strand']) + 1, transcript['exonCoord'], transcript['strand'])
            if abs(stopCodonMidCoord - breakpoints[1]) == 1:
                stopCodonSeries1['exonCoord'] = (breakpoints[1], stopCodonMidCoord) if transcript['strand'] == '+' else (stopCodonMidCoord, breakpoints[1])
                stopCodonSeries2['exonCoord'] = (stopCodonEndCoord, stopCodonEndCoord)
                stopCodonDF1, stopCodonDF2 = stopCodonSeries1.to_frame().T, stopCodonSeries2.to_frame().T
                stopCodonDF1['biotype'], stopCodonDF1['phase'] = 'stop_codon', 0
                stopCodonDF2['biotype'], stopCodonDF2['phase'] = 'stop_codon', 1
            else:
                stopCodonSeries1['exonCoord'] = (breakpoints[1], breakpoints[1])
                stopCodonSeries2['exonCoord'] = (stopCodonMidCoord, stopCodonEndCoord) if transcript['strand'] == '+' else (stopCodonEndCoord, stopCodonMidCoord)
                stopCodonDF1, stopCodonDF2 = stopCodonSeries1.to_frame().T, stopCodonSeries2.to_frame().T
                stopCodonDF1['biotype'], stopCodonDF1['phase'] = 'stop_codon', 0
                stopCodonDF2['biotype'], stopCodonDF2['phase'] = 'stop_codon', 2
            stopCodonDF = pd.concat([stopCodonDF1, stopCodonDF2])
    else:
        stopCodonDF = pd.DataFrame(columns=['txID', 'annoStatus', 'geneID', 'chrom', 'strand', 'exonCoord', 'exonNum', 'biotype', 'phase'])

    # Assemble 3' UTR annotations
    if cdsEndIdx is not None:
        utr3DF = myDF.iloc[(cdsEndIdx+1):].copy() if transcript['strand'] == '+' else myDF.iloc[:cdsEndIdx].copy()
        if transcript['strand'] == '+' and myDF.iloc[cdsEndIdx]['exonCoord'][1] > breakpoints[1]:
            utr3Append = myDF.iloc[cdsEndIdx].copy()
            utr3Append['exonCoord'] = (breakpoints[1], utr3Append['exonCoord'][1])
            utr3DF = pd.concat([utr3Append.to_frame().T, utr3DF])
        elif transcript['strand'] == '-' and myDF.iloc[cdsEndIdx]['exonCoord'][0] < breakpoints[1]:
            utr3Append = myDF.iloc[cdsEndIdx].copy()
            utr3Append['exonCoord'] = (utr3Append['exonCoord'][0], breakpoints[1])
            utr3DF = pd.concat([utr3DF, utr3Append.to_frame().T])
        utr3DF['biotype'], utr3DF['phase'] = 'UTR', '.'
    else:
        utr3DF = pd.DataFrame(columns=['txID', 'annoStatus', 'geneID', 'chrom', 'strand', 'exonCoord', 'exonNum', 'biotype', 'phase'])
    
    # Combine all five pieces together
    outDF = pd.concat([utr5DF, startCodonDF, cdsDF, stopCodonDF, utr3DF]).reset_index(drop=True)

    # Re-structure outDF as a GTF
    outDF[['exonStart', 'exonEnd']] = pd.DataFrame(outDF['exonCoord'].tolist())
    outDF['score'], outDF['info'] = '.', 'gene_id \"' + outDF['geneID'] + '\"; transcript_id \"' + outDF['txID'] + '\"; exon_number ' + outDF['exonNum'].astype(str) + ';'
    outDF = outDF[['chrom', 'annoStatus', 'biotype', 'exonStart', 'exonEnd', 'score', 'strand', 'phase', 'info']]
    outDF.columns = list(range(outDF.shape[1]))

    return outDF

def predictORF(inDF, annoDF, annoStartCodon, canonProteinSeq, genome):
    '''
    This is a function to update inDF with ORF annotations for each input transcript (where applicable)
    '''
    # Generate a copy of inDF and extract gene/transcript IDs from the INFO field (last column)
    outDF = inDF.copy()
    outDF['txID'] = outDF[8].apply(lambda x: [item.split('"')[1] for item in x.split(';') if 'transcript_id' in item][0])
    outDF['geneID'] = outDF[8].apply(lambda x: [item.split('"')[1] for item in x.split(';') if 'gene_id' in item][0])

    # Blacklist transcripts where there exists at least one feature whose end coordinate is smaller than the start coordinate
    blacklist = set(outDF[outDF[4] < outDF[3]]['txID'])
    outDF = outDF[~outDF['txID'].isin(blacklist)]

    # Reconstruct the INFO field of inDF to include gene_id followed by transcript_id
    outDF[8] = 'gene_id \"' + outDF['geneID'] + '\"; transcript_id \"' + outDF['txID'] + '\";'

    # Only retain exon-level entries for outDF; then group exons by transcript ID
    txDF = outDF[outDF[2] == 'exon'].groupby([0, 1, 2, 5, 6, 7, 8, 'txID', 'geneID'])[[3, 4]].agg(lambda x: sorted(list(x))).reset_index()
    
    # Determine exon numbering by transcript strand
    txDF['exonNumber'] = txDF.apply(lambda x: list(range(1, 1+len(x[3]))) if x[6] == '+' else list(range(len(x[3]), 0, -1)), axis = 1)
    exonDF = txDF.explode([3, 4, 'exonNumber']).copy()
    exonDF[8] = exonDF[8] + ' exon_number ' + exonDF['exonNumber'].astype(str) + ';'
    outDF = pd.concat([outDF[outDF[2] == 'transcript'], exonDF[outDF.columns]])

    # Combine exon start coordinates with exon end coordinates
    txDF['exonCoord'] = txDF.apply(lambda x: [(x[3][i], x[4][i]) for i in range(len(x[3]))], axis = 1)
    txDF = txDF.drop([2, 5, 7, 8, 3, 4, 'exonNumber'], axis = 1)
    txDF.columns = ['chrom', 'annoStatus', 'strand', 'txID', 'geneID', 'exonCoord']

    outlog('Found a total of ' + str(txDF.shape[0]) + ' transcripts in input GTF file...')

    # Determine number of exons and re-order columns of txDF
    txDF['numExon'] = txDF['exonCoord'].apply(len)
    txDF = txDF[['txID', 'numExon', 'annoStatus', 'geneID', 'chrom', 'strand', 'exonCoord']]

    # Restructure outDF to only contain columns needed in the final GTF
    outDF = outDF.drop(['txID', 'geneID'], axis = 1)

    # Iterate over transcripts in txDF
    appendDF = pd.DataFrame(columns=list(range(9)))

    appendItems = []
    for idx, transcript in txDF.iterrows():
        # Try translating from any ATG in the transcript sequence but prioritize those that are annotated
        anyRes = getLongestORF(transcript, annoStartCodon, canonProteinSeq, genome)
        if anyRes is not None:
            appendItems += injectCDS(transcript, anyRes).values.tolist()

        # Send update message for every 10 entries processed
        if (idx + 1) % 10 == 0:
            outlog('Finished detecting ORFs for ' + str(idx + 1) + ' transcripts...')
    
    # Merge outDF with appendItems; then sort rows of outDF by genomic coordinates
    outDF = pd.concat([outDF, pd.DataFrame(appendItems)]).sort_values(by = [0, 3, 4]).reset_index(drop = True)

    return outDF

def main():
    message = 'Updates a GTF file to include ORF annotations for transcripts'
    parser = argparse.ArgumentParser(description = message)

    # Add arguments
    parser.add_argument('-i', metavar = '/path/to/transcript/assembly/GTF', required = True,
        help = 'path to GTF file describing genomic features of discovered transcripts')
    parser.add_argument('-a', metavar = '/path/to/GENCODE/reference/GTF', required = True,
        help = 'path to GTF file containing reference gene annotations from GENCODE')
    parser.add_argument('-f', metavar = '/path/to/reference/genome/FASTA', required = True,
        help = 'path to FASTA file containing reference genome sequence')
    parser.add_argument('-o', metavar = '/path/to/output/GTF', required = True,
        help = 'path to output GTF')

    # Parse command-line arguments
    args = parser.parse_args()
    infile, anno, reference, outfile = args.i, args.a, args.f, args.o

    # Read in infile and anno as Pandas dataframes
    inDF = pd.read_csv(infile, sep = '\t', header = None, comment = '#')
    inDF = inDF[inDF[2] != 'gene']
    annoDF = pd.read_csv(anno, sep = '\t', header = None, comment = '#')

    # Read in reference as a pysam FastaFile object
    genome = pysam.FastaFile(reference)

    # Extract genomic start coordinates of annotated start codons in annoDF
    annoStartCodon = getAnnoStartCodon(annoDF)

    # Extract canonical protein sequences for all protein-coding genes in annoDF
    canonProteinSeq = getCanonProteinSeq(annoDF, genome)

    # Predict open reading frames for each transcript in inDF
    outlog('Annotating ORFs for input transcripts...')
    orfDF = predictORF(inDF, annoDF, annoStartCodon, canonProteinSeq, genome)

    # Write orfDF to outfile
    outlog('Writing updated transcript annotations to output file...')
    orfDF.to_csv(outfile, sep = '\t', index = False, header = False, quoting = csv.QUOTE_NONE)

    # Close genome FastaFile object
    genome.close()

if __name__ == '__main__':
    main()

