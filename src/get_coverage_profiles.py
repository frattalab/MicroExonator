#!/usr/bin/env python3

'''
Create coverage files (BigWig/Bedgraph) from BAM files produced in round2
Microexonator's custom tags prevent using SAMs/BAMs on top of standard human genome for visualisation
Use pysam pileup to determine read coverage across sequence tags and reformat these to 'standardised' format e.g.
chr1 10001 10002 4
chr1 10002 10003 3
etc (for all microexons)


Notes:
https://pysam.readthedocs.io/en/latest/faq.html#pysam-coordinates-are-wrong
Pysam reports all coordinates & intervals following Python's 0-based, half-open notation

https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile.pileup
first base returned by pileup function is the first base of the first read
NOT guaranteed to be the first base or query region

'''

import sys
import pysam
import os
from collections import OrderedDict

def ref_name_add_brackets(string):
    '''
    returns string of reference sequence name with brackets added around integer after 'MULTI'
    Removed '(' or ')' from the reference sequence names in header of BAM file so reference sequence names match SAM specification regex (see Round2.skm sam_to_bam)
    But corresponding reference sequence name in alignment lines (column 3) still contains brackets
    This helper function adds the brackets back into str
    Example: ref_name_add_brackets('chr6:61732764-61894634|MULTI2|100_CATTGTGATAACGAAGTACCACAAACTGG_100')
    returns chr6:61732764-61894634|MULTI(2)|100_CATTGTGATAACGAAGTACCACAAACTGG_100
    '''
    split_str = string.split('|')
    #2nd element of split_str - |MULTI<number> | is offending party need to add back in
    fixed_multi = split_str[1][:5] + '(' + split_str[1][6:] + ')'

    return '|'.join(split_str[0], fixed_multi, split_str[2])

def initial_pileup(bam):
    '''
    Calculate per nucleotide coverage for all sequence tags using pysam pileup
    Return nested dict of {Reference name: {position: coverage}}
    (first pass - returns values for bases with coverage only)
    '''
    pileup_dict = {}

    #bam.references returns large list of reference sequence names
    #Try bundling reference_tags into generator expression, so not all sequence names are stored in memory
    for reference_tag in (ref for ref in bam.references):

        if 'MULTI' in reference_tag:
            # This means I removed '(' or ')' from the header of BAM file (see Round2.skm sam_to_bam)
            # bam.references fetches references sequence names from @SQ SN:<name> lines in header (by looks of it...)
            # I did not edit the Reference sequence names in column 3 of SAM/BAM alignment lines
            # pileup would return empty values if don't add brackets back into sequence name for query

            if '(' not in reference_tag and ')' not in reference_tag:
                #add brackets back in
                fixed_ref_tag = ref_name_add_brackets(reference_tag)
                reference_cov = {pileupcol.pos: pileupcol.n for pileupcol in bam.pileup(fixed_ref_tag)}
                pileup_dict[reference_tag] = reference_cov

            else:
                #proceed as normal
                reference_cov = {pileupcol.pos: pileupcol.n for pileupcol in bam.pileup(reference_tag)}
                pileup_dict[reference_tag] = reference_cov

        else:
            reference_cov = {pileupcol.pos: pileupcol.n for pileupcol in bam.pileup(reference_tag)}
            pileup_dict[reference_tag] = reference_cov

    return pileup_dict

def fill_pileup_dict(pileup_dict):
    '''
    Returns complete nested dict of {reference_dict: OrderedDict{position: coverage}} with position & coverage key:value pair for every position in reference tag
    pysam.pileup returns 0-based coordinates - this dictionary preservess this convention
    OrderedDict stores keys in ascending order (i.e. position 1 is first in dict and represented by key 0 etc.)
    '''

    def seq_range_from_tag(string):
        '''
        Return list of all positions for reference sequence from reference sequence tag
        e.g. sequence tag chrY:2966939-2981837|ENST00000417305.1|100_96
        has a reference tag of 196 nucleotides/bases
        return list of positions of [1..196]
        e.g. sequence tag chrY:5501255+5737271|ENST00000400457.3|100_CTTTCATACCTGGACTAAAGAAAG_100
        has a reference tag of 224 nucleotides/bases
        returns list of positions [1..224]
        '''

        length_tag = string.split("|")[2]
        length_tag = length_tag.split("_")

        if len(length_tag) == 2:
            tag_length = int(length_tag[0]) + int(length_tag[1])

        elif len(length_tag) == 3:
            #2nd element is string representing microexon insertion
            tag_length = int(length_tag[0]) + len(length_tag[1]) + int(length_tag[2])

        #positions in this list are 0 based i.e. 0 corresponds to first base in reference tag
        positions_list = [pos for pos in range(tag_length)]

        return positions_list

    #nested dict to return at end - nested dictionary for each reference tag is OrderedDict of with position 1 first key
    complete_pileup_dict = {}

    for reference_tag, initial_cov_dict in pileup_dict.items():
        # all nucleotide positions want a coverage value for in OrderedDict for reference tag
        complete_cov_dict = OrderedDict()
        full_reference_range = seq_range_from_tag(reference_tag)

        # check if position has coverage value reported in initial_cov_dict (from initial_pileup)
        # if it does, extract value and report in OrderedDict for reference_tag
        # if not then no reads cover that position
        for position in full_reference_range:
            if position in initial_cov_dict.keys():
                complete_cov_dict[position] = initial_cov_dict.get(position)
            else:
                complete_cov_dict[position] = 0

        # add reference sequence name and completed coverage dict to complete dictionary to return at end of function
        complete_pileup_dict[reference_tag] = complete_cov_dict

    return complete_pileup_dict



if __name__ == '__main__':

    #bamfile = "/home/sam/cluster/sbs_projects/MicroExonator_fork/MicroExonator/Round2/Ward_CTL_1.bam"
    bamfile = sys.argv[1]
    bamfile = pysam.AlignmentFile(bamfile,"rb")

    print_list = ["chr16:28152785-28156073|MULTI2|100_AATATGGGGCTCCTGGTGAGGAACAGAAAG_100","chrY:9545006-9545124|ENST00000421178.2|100_90","chr17:42981074+42987255|ENST00000361677.5|100_GGCTTG_100"]

    #generate dictionary of {reference_sequence_name: {position: n_reads}} using pysam.pileup
    bam_pileup_dict = initial_pileup(bamfile)

    for ref_name, pileup_dict in bam_pileup_dict.items():
        if ref_name in print_list:
            print(ref_name,pileup_dict)

    bamfile.close()


    # Initial pileup dict only returns positions where coverage > 0
    # Need to generate coverage values for all positions in sequence tag
    bam_pileup_dict = fill_pileup_dict(bam_pileup_dict)

    #Now have nested dict of {ref_tag_name: OrderedDict(position: coverage)} where OrderedDict stores positions in ascending order!
