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
import csv
import sys
import pysam
import os
from collections import OrderedDict, Counter


######-------------------------
# Helper functions for main functions contained within this block
######--------------------------

def get_data_dict(row_list,key_idx_dict={'transcript_id': 1, 'splice_junction': 3,'mex_sequence': 6}):
    '''
    returns dictionary of data of choice to return from row of microexon results table (row data, row_list, is stored in list by csv.reader).
    Key names and index of column should be provided in key_idx_dict
    '''
    data_dict = {key: row_list[val] for key, val in key_idx_dict.items()}
    return data_dict

def mexonator_csv_to_dict(path):
    '''
    return nested dictionary of {microexon_coords: {'splice_junction': coords, 'mex_sequence': string, transcript = 'transcript'}} from microexonator results table
    Table columns to extract are 1st (microexon_coords), 2nd (transcript id), 4th (splice_junction_coords) & 7th (mex_sequence)
    '''

    with open(path) as infile:
        table_dict = {row[0]: get_data_dict(row) for row in csv.reader(infile, delimiter = "\t") if row[0] != 'ME'}
        #!= ME my hacky way to skip header line
    return table_dict


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
    fixed_multi = split_str[1][:5] + '(' + split_str[1][5:] + ')'

    return '|'.join([split_str[0], fixed_multi, split_str[2]])


def filter_me_reads_list(round2_filter_path):
    '''
    Returns list of valid microexon reads produced by Round2_filter rule - only these reads will contribute to coverage in bedgraph
    '''
    #read ids are in 1st column of SAM file format/output from Round2_filter
    with open(round2_filter_path) as infile:
        reads = [line.split('\t')[0] for line in infile]

    sys.stderr.write("{0}\n".format('\n'.join(reads)))
    return reads

#############----------------------
# Main functions - loosely chronological order (check main for final order)
#############----------------------


def all_microexons_dict(filepaths_list, expected_files_list=['out.ambiguous.txt','out_filtered_ME.txt','out_low_scored_ME.txt','out_shorter_than_3_ME.txt','out.high_quality.txt']):
    '''
    Returns nested dictionary of all microexons as {microexon_coords: {'splice_junction': coords, 'mex_sequence': string, transcript = 'transcript'}}
    Collated from all output files under Report/
    All files listed in expected_files_list should have the same 'column indexes' for data want to extract
    '''
    full_dict = {}
    for path in filepaths_list:
        if os.path.basename(path) not in expected_files_list:
            raise Exception("{0} could not be found in expected_files_list - {1} - which are known to have same column indexes. Was the correct file passed?".format(path,','.join(expected_files_list)))
        else:
            #dict of {microexon_coords: {'splice_junction': coords, 'mex_sequence': string, transcript = 'transcript'}} for given filepath
            path_dict = mexonator_csv_to_dict(path)

            for mex, mex_data in path_dict.items():
                if mex not in full_dict.keys():
                    full_dict[mex] = mex_data
                else:
                    #microexon has already been added to full dict - check if mex_data is identical to dict associated with mex in full_dict
                    #if same skip to next microexon, otherwise stop execution
                    if mex_data.get('splice_junction') == full_dict.get(mex).get('splice_junction'):
                        if mex_data.get('mex_sequence') == full_dict.get(mex).get('mex_sequence'):
                            continue
                        else:
                            raise Exception("microexon {0} has differently reported microexon sequences in different output files - please clarify".format(mex))

                    elif mex_data.get('mex_sequence') == full_dict.get(mex).get('mex_sequence'):
                        if mex_data.get('splice_junction') == full_dict.get(mex).get('splice_junction'):
                            continue
                        else:
                            raise Exception("microexon {0} has differently reported splice junctions in different output files - please clarify".format(mex))
                    else:
                        raise Exception("microexon {0} has differently reported splice junctions and sequences throughout different output files. please check".format(mex))

    return full_dict


def get_mex_reference_names(microexon_dict,pysam_bam):
    '''
    Return dict of {(sj_coords, me_seq): 'ref_seq_name'} mapping identified microexons to reference sequence names
    Do this to reduce iterations of pysam.pileup - only really interested in coverages profiles across microexons
    Reference names have general format:
    chr17:42981074+42987255|ENST00000361677.5|100_GGCTTG_100
    splice_junction | tr_id | 5' upstream extension _ microexon_seq _ 3' upstream extension

    splice_junction ('microexon_coords') values inside nested dict of microexon_dict would be represented as
    chr17:42981074+42987255 - can match to ref_name.split('|')[0]

    microexon sequence ('mex_sequence') values inside nested dict would be represented as
    GGCTTG - can match to ref_name.split('|')[2].split('_')[1]

    if both satisfied then return reference sequence name in dict
    '''

    # 1. Generate list of tuples of (splice_junction, mex_sequence) from nested dicts in microexon_dict
    #make this nested dict of {microexon: (splice_junction,mex_sequence)}
    #sj_seq_tuples = [tuple([nested_dict.get('splice_junction'), nested_dict.get('mex_sequence')]) for nested_dict in microexon_dict.values()]
    sj_seq_tuples = {tuple([nested_dict.get('splice_junction'), nested_dict.get('mex_sequence')]): microexon for microexon, nested_dict in microexon_info_dict.items()}

    # 2. Iterate over tuple of reference sequence names, report if splice junctions and microexon sequences match any tuples in sj_seq_tuples
    # also initiate dict where will store {(splice_junction,me_seq): microexon_coords}
    ref_seq_dict = {}
    sj_seq_mex_dict = {}

    for ref_name in pysam_bam.references:
        #this should then iterate over .items of dict
        #if all conditions matched, update ref_seq_dict and tuple_microexon_dict (tuple: microexon(key))
        for mex_tup, microexon in sj_seq_tuples.items():
            #sj_coords can contain multiple entries, separated by comma - check each one
            if ',' in mex_tup[0]:
                sj_coords = mex_tup[0].split(',')
                for coords in sj_coords:
                    #do (1) sj coords and (2) mex_sequences match?
                    if coords == ref_name.split('|')[0] and mex_tup[1] == ref_name.split('|')[2].split('_')[1]:
                        ref_seq_dict[tuple([coords, mex_tup[1]])] = ref_name
                        sj_seq_mex_dict[tuple([coords, mex_tup[1]])] = microexon
            else:
                #do (1) sj coords and (2) mex_sequences match?
                if mex_tup[0] == ref_name.split('|')[0] and mex_tup[1] == ref_name.split('|')[2].split('_')[1]:
                    ref_seq_dict[mex_tup] = ref_name
                    sj_seq_mex_dict[mex_tup] = microexon
                else:
                    continue

    return ref_seq_dict, sj_seq_mex_dict


def initial_pileup(bam, seq_name_list, valid_reads_path):
    '''
    Calculate per nucleotide coverage for all reference sequence names in seq_name_list using pysam pileup
    Return nested dict of {Reference name: {position: coverage}}
    (first pass - returns values for bases with coverage only)
    '''
    #valid read ids that do not have a primary alignment to the whole genome
    valid_reads_list = filter_me_reads_list(valid_reads_path)

    #bam.references returns large list of reference sequence names
    #Try bundling reference_tags into generator expression, so not all sequence names are stored in memory
    pileup_dict = {}
    for reference_tag in (ref for ref in seq_name_list):
        # only want reads to count towards coverage if also found in valid_reads_list
        reference_cov = {pileupcolumn.pos: (sum(1 for read in pileupcolumn.get_query_names() if read in valid_reads_list)) for pileupcolumn in bam.pileup(reference_tag)}

        #reference_cov = {pileupcol.pos: pileupcol.n for pileupcol in bam.pileup(reference_tag)}
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

def pileup_to_genome_coordinates(pileup_dict,tuples_reference_dict,tuples_microexon_dict):
    '''
    return dict of {reference_tag: {genome_coord: coverage}} where coordinate integer corresponds to actual coordinate (NOT an index to coordinate)
    For each reference tag, convert index based coverage data to corresponding genome coordinates
    reference tag chrY:5501255+5737271|ENST00000400457.3|100_CTTTCATACCTGGACTAAAGAAAG_100
    sj_seq_tuple and microexon_coords ('chrY:5501255+5737271', 'CTTTCATACCTGGACTAAAGAAAG'): 'chrY_+_5581774_5581798'
    Assumptions: coordinates follow BED conventions
    '''
    genome_cov_dict = {}

    for sj_seq_tuple, ref_name in tuples_reference_dict.items():
        # 1. Parse out 5' upstream & 3' downstream of sj distances from end of ref_name
        upstream_ext_length = int(ref_name.split('|')[2].split('_')[0])
        downstream_ext_length = int(ref_name.split('|')[2].split('_')[-1])

        #2. extract 5' splice junction coord and 3' splice junction coord from start of ref name
        #defining as 5' & 3' based on genomic coordinate (not relative to gene/strand)
        if '+' in sj_seq_tuple[0]:
            sj_coord5 = int(sj_seq_tuple[0].split('+')[0].split(':')[1])
            sj_coord3 = int(sj_seq_tuple[0].split('+')[1])
        elif '-' in sj_seq_tuple[0]:
            sj_coord5 = int(sj_seq_tuple[0].split('-')[0].split(':')[1])
            sj_coord3 = int(sj_seq_tuple[0].split('-')[1])

        #3. Get range of genomic coordinates for regions upstream of 5' splice coord & downstream of 3' splice coord
        #range(1,5) includes 1 but not 5, but BED conventions state genomic range of 1-5 encompasses nucleotides 2-5 (0-based start, 1-based end)
        #To go from pythonic to (BED) genomic - need to add one to both 'start' & 'end' (i.e. for grange 1-5, pass range 2-6 to return 2,3,4,5)
        genome_coords5 = [coord for coord in range((sj_coord5 - upstream_ext_length) + 1, sj_coord5 + 1)]
        genome_coords3 = [coord for coord in range((sj_coord3 + 1), (sj_coord3 + downstream_ext_length) + 1)]

        # 4. Get genome coordinates of microexon given splice junction & sequence tuples
        mex_genome_coords_str = tuples_microexon_dict.get(sj_seq_tuple)

        #Again assume provided microexon coordinates follow BED12 convention - start is 0 based and end is 1-based
        mex_coord_start = int(mex_genome_coords_str.split('_')[-2]) + 1
        mex_coord_end = int(mex_genome_coords_str.split('_')[-1])

        #for range() function to include mex_coord_end value, need to add 1
        mex_genome_coords = [coord for coord in range(mex_coord_start, mex_coord_end + 1)]

        # 5. Join lists together in most 5' - most 3' order - indexes in coverage dict start from leftmost genome coordinate in reference_tag
        genome_coords5.extend(mex_genome_coords)
        genome_coords5.extend(genome_coords3)

        #6. Map genome coordinates to coverage values in nested OrderedDict in pileup_dict for given reference name
        #Expect number of genome coordinates to match number of keys (coordinates) in OrderedDict
        coverage_dict = pileup_dict.get(ref_name)

        if len(genome_coords5) != len(coverage_dict.values()):
            raise Exception("Number of generated genome coordinates does not match the number of coordinates computed for reference region. Double check how generating genome coordinates")

        else:
            coords_cov_dict = OrderedDict()
            #use enumerate over reference's ordered dict - use idx to extract genome coordinate from genome_coords5 and key to get coverage value
            for idx, key in enumerate(coverage_dict):
                #coordinate: coverage
                coords_cov_dict[genome_coords5[idx]] = coverage_dict.get(key)

            #ref_name: coords_cov_dict
            genome_cov_dict[ref_name] = coords_cov_dict

    return genome_cov_dict

'''
def resolve_coverage_conflicts(chrs_coord_ref_dict, conflicting_tags_dict, threshold=0.5):
    #
    #chrs_coord_ref_dict = {chr :{coord: reference_name}}
    #conflicting_tags_dict = {(tag1, tag2): {coord: (tag1_cov, tag2_cov)}} where tag already present in chrs_coord_tag_dict
    #Some reference tags share common genomic coordinates but differ in read coverage at these positions
    #Function selects reference tag in which majority (threshold) of positions have higher coverage than overlapping tag
    #Returns dictionary of {<current_tag_name>: <tag to replace>}
    #

    # 1. Reference tag (1st position in tuple key) could have conflicts with multiple sequence tags (i.e. multiple tuples have same tuple[0] tag)
    # Generate list of tuples of [(tag1,tag2),(tag1,tag3), (tag4,tag6)]} from keys in conflicting_tags_dict
    #tag_conflict_pairs_list = [conflict_tup]

    #for conflict_tuple in conflicting_tags_dict.keys():
    #    if conflict_tuple[0] not in tag_conflict_pairs_dict:
    #         tag_conflict_pairs_dict[conflict_tuple[0]] = [conflict_tuple]
    #    else:
    #        tag_conflict_pairs_dict[conflict_tuple[0]].append(conflict_tuple)

    #2. Iterate over each combo of conflicting tags in tag_conflict_pairs_dict to pick tag with consistently highest coverage over conflicting coordinates
    # In this way select the tag/microexon and host transcripts that is best supported by reads
    # generate dict of {<current_tag>: <tag_to_replace>} where current_tag is name of tag present in conflicting_tags_dict
    tags_to_replace_dict = {}

    for conflict_tuple, coverage_dict in conflicting_tags_dict.items():
        if conflict_tuple[0] not in tags_to_replace_dict:
            # as it stands tuple[0] (existing tag) does 'not need replacing' or is not worse supported than another conflicting tag
            #how many positions is coverage for coordinates in existing reference tag >= coverage in conflicting reference tag?
            n_pos_higher = sum(1 for coverage_tuple in coverage_dict.values() if coverage_tuple[0] >= coverage_tuple[1])

            if n_pos_higher / len(coverage_dict.keys()) >= threshold:
                sys.stderr.write("{0} selected over {1} to represent conflicting coordinates as it has greater coverage over at least {2} conflicting positions\n".format(conflict_tuple[0], conflict_tuple[1], threshold))
                continue
            else:
                #conflicting tag not present in existing coverage dict has consistently higher coverage than existing tag - want to replace
                tags_to_replace_dict[conflict_tuple[0]] = conflict_tuple[1]
        else:
            #Need to compare to highest coverage tag for position NOT conflict_tuple[0] (already proved to be more poorly supported than another tag)
            current_tag = tags_to_replace_dict.get(conflict_tuple[0])
            current_tag_covs_dict = conflicting_tags_dict.get(tuple([conflict_tuple[0], current_tag]))

            #how many coordinates are shared between the comparison reference tags?
            compare_tag_coords = [tup[1] for tup in coverage_dict.values()]
            shared_coords = [coord for coord in current_tag_covs_dict.keys() if coord in compare_tag_coords]
'''


def coord_coverage_to_chrs(coord_pileup_dict):
    '''
    return dictionary of {chr :{coord: coverage}} from coord_pileup dict of {reference_name: OrderedDict((coord: coverage))}
    return dictionary of {chr :{coord: reference_name}} - track which reference name is contributing coverage value at that coordinate
    returned as tuple of coverage dictionaries (coverage then reference name)
    e.g. reference tag chrY:5501255+5737271|ENST00000400457.3|100_CTTTCATACCTGGACTAAAGAAAG_100
    '''

    chrs_cov_dict = {} # {chr: {coord: coverage}}
    chrs_coord_tag_dict = {} # {chr: {coord: [reference_name]}}
    #conflicting_tags_dict = {} # {(tag1, tag2): {coord: (tag1_cov, tag2_cov)}} where reference_name already present in chrs_coord_tag_dict

    for ref_name, coverage_dict in coord_pileup_dict.items():
        chr = ref_name.split(':')[0]

        if chr not in chrs_cov_dict:
            chrs_cov_dict[chr] = dict(coverage_dict) #coverage_dict is OrderedDict - want a normal one
            chrs_coord_tag_dict[chr] = {coord: [ref_name] for coord in coverage_dict.keys()}
            continue
        else:
            for coord, coverage in coverage_dict.items():
                if coord not in chrs_cov_dict.get(chr):
                    chrs_cov_dict[chr][coord] = coverage
                    chrs_coord_tag_dict[chr][coord] = [ref_name]

                else:
                    #tags cover same genomic coordinates but have different read coverages
                    #As valid reads correspond to primary alignments to reference tag (and don't align to genome) - know they are distinct reads
                    #Can sum coverage at that position
                    new_coverage = chrs_cov_dict[chr][coord] + coverage
                    chrs_cov_dict[chr][coord] = new_coverage
                    chrs_coord_tag_dict[chr][coord].append(ref_name)


                    #if tuple([chrs_coord_tag_dict[chr][coord], ref_name]) not in conflicting_tags_dict:
                    #    conflicting_tags_dict[tuple([chrs_coord_tag_dict[chr][coord], ref_name])] = {coord: tuple([chrs_cov_dict[chr][coord], coverage])}
                    #else:
                    #    conflicting_tags_dict[tuple([chrs_coord_tag_dict[chr][coord], ref_name])][coord] = tuple([chrs_cov_dict[chr][coord], coverage])



    #conflicting_tags_dict - how many conflicting regions?
    #sys.stderr.write("{0} regions covered by multiple reference sequence tags have conflicting coverages\n".format(len(set([key[0] for key in conflicting_tags_dict.keys()]))))
    #For all regions, how many regions have conflicts between multiple pairs of reference tags?
    #1. how many times/conflicts does each tag1 appear in conflicting
    #2. how many tags have n conflicts
    #conflict_counts = Counter(Counter([key[0] for key in conflicting_tags_dict.keys()]).values())
    #for n_conflicts, n_regions in conflict_counts.items():
    #    sys.stderr.write("{0} regions/sequence tags have {1} conflicting tag pairs\n".format(n_regions, n_conflicts))



                    #sys.stderr.write("conflicting_coverages\t{1}\t{0}\t{3}\t{2}\t{5}\t{4}\n".format(coord, chr, chrs_cov_dict[chr][coord], chrs_coord_ref_dict[chr][coord], coverage, ref_name))
                    #raise Exception("coordinate {0} on {1} has different read coverages across reference tags {2} & {3}".format(coord, chr, chrs_coord_ref_dict[chr][coord], ref_name))

    return chrs_cov_dict, chrs_coord_tag_dict




def chrs_coverage_to_bed(chrs_cov_dict,header="type=bedGraph"):
    '''
    Print coverage values in tab-delimited BEDgraph half-open coordinate format (start 0 based, end 1 based)
    <chr> \t <start> \t <end> \t <coverage value>
    Chromosomes and coordinates sorted in ascending order (i.e. smallest first)
    chr1..chr<last> then named chrs in alphabetical order
    Takes dict of {chr: {coordinate: coverage} }, where coordinate key represents value of corresponding genome coordinate
    '''

    numbered_chrs = []
    lettered_chrs = []

    for chr in chrs_cov_dict.keys():
        try:
            numbered_chrs.append(int(chr.split('chr')[-1]))
        except ValueError:
            lettered_chrs.append(chr.split('chr')[-1])

    # Before printing values, need a header line (currently minimum is "type=bedGraph")
    print(header)

    #first print ascending order sorted numbered_chrs
    for chr in sorted(numbered_chrs):
        coverage_dict = chrs_cov_dict.get('chr'+str(chr))
        # items tuple approach coverage dict will sort keys (coordinates, integer) in ascending order
        for coord, cov in sorted(coverage_dict.items()):
            #BED are half open - so for coverage to represent actual coordinate need coord to be end (and start is coord -1)
            print('\t'.join(['chr'+str(chr), str(coord - 1), str(coord), str(cov)]))

    #then add on lettered chromosomes (in alphabetical order)
    for chr in sorted(lettered_chrs):
        coverage_dict = chrs_cov_dict.get('chr'+str(chr))
        # items tuple approach coverage dict will sort keys (coordinates, integer) in ascending order
        for coord, cov in sorted(coverage_dict.items()):
            #BED are half open - so for coverage to represent actual coordinate need coord to be end (and start is coord -1)
            print('\t'.join(['chr'+str(chr), str(coord - 1), str(coord), str(cov)]))



    #for chr, coverage_dict in sorted(chrs_cov_dict.items()):
        #chromosomes sorted in ascending order - now for each chromosome need coordinates to be sorted in ascending order
        #tuple approach again with coverage dict will sort keys (coordinates) in ascending order
    #    for coord, cov in sorted(coverage_dict.items()):
            #BED are half open - so for coverage to represent actual coordinate need coord to be end (and start is coord -1)
    #        print('\t'.join([chr, str(coord - 1), str(coord), str(cov)]))


if __name__ == '__main__':

    #bamfile = "/home/sam/cluster/sbs_projects/MicroExonator_fork/MicroExonator/Round2/Ward_CTL_1.bam"
    bamfile = sys.argv[1]
    results_csvs = sys.argv[2]
    processed_reads_path = sys.argv[3]

    #results csvs are all microexon results tables (under Report/), passed to script as comma-separated string
    results_csvs = results_csvs.split(',')

    #load results tables into dictionary of {microexon_coords: {'splice_junction': coords, 'mex_sequence': string, transcript = 'transcript'}}
    microexon_info_dict = all_microexons_dict(results_csvs)
    #print(microexon_info_dict)
    #print(len(microexon_info_dict.keys()))

    #dict of {(sj_coords, me_seq): 'ref_seq_name'} mapping identified microexons to reference sequence names - values passed to pysam.pileup() to generate intiial coverages
    bamfile = pysam.AlignmentFile(bamfile,"rb")

    #(1) = {(sj_coords, me_seq): 'ref_seq_name'}, (2) = {(sj_coords, me_seq): 'microexon'}
    microexon_ref_seq_names, sj_seq_microexon_dict  = get_mex_reference_names(microexon_info_dict,bamfile)
    #print(microexon_ref_seq_names)
    #print(len(microexon_ref_seq_names.values()))
    #print(sj_seq_microexon_dict)

    #print_list = ["chrY:9545006-9545124|ENST00000421178.2|100_90","chr17:42981074+42987255|ENST00000361677.5|100_GGCTTG_100","chr2:98546694+98552785|ENST00000409851.8|100_CAGTTTTGAGGAGTGTTG_100", "chr7:128393028-128394277|ENST00000469328.5|100_CTATACCTTTCTGCCGT_100"]
    #print_tuple_list = [key for key, val in microexon_ref_seq_names.items() if val in print_list]
    #print_dict = {key: val for key, val in microexon_ref_seq_names.items() if val in print_list}
    #sj_print_dict = {key: val for key, val in sj_seq_microexon_dict.items() if key in print_tuple_list}


    bam_pileup_dict = initial_pileup(bamfile, microexon_ref_seq_names.values(), processed_reads_path)
    bamfile.close()

    #for ref in print_list:
    #    print(ref, bam_pileup_dict.get(ref))

    # Initial pileup dict only returns positions where coverage > 0
    # Need to generate coverage values for all positions in sequence tag
    bam_pileup_dict = fill_pileup_dict(bam_pileup_dict)

    #for ref_name in print_list:
    #    print(ref_name,bam_pileup_dict.get(ref_name))

    # convert coverages positions for reference tags (currently stored as index) to genome coordinates for tag
    #{ref_name: OrderedDict((coord, coverage))}
    #bam_pileup_dict = pileup_to_genome_coordinates(bam_pileup_dict, print_dict, sj_print_dict)
    bam_pileup_dict = pileup_to_genome_coordinates(bam_pileup_dict,microexon_ref_seq_names, sj_seq_microexon_dict)

    #for ref in print_list:
    #    print(ref,bam_pileup_dict.get(ref))

    #get dictionary of {chr: {coord: coverage}} & {chr: {coord: ref_name}} from reference sequence names
    chrs_pileup_dict, chrs_pileup_references_dict = coord_coverage_to_chrs(bam_pileup_dict)
    #print(chrs_pileup_dict)

    #prints to STDOUT
    chrs_coverage_to_bed(chrs_pileup_dict)

    '''





    #generate dictionary of {reference_sequence_name: {position: n_reads}} using pysam.pileup
    bam_pileup_dict = initial_pileup(bamfile)

    for ref_name in print_list:
        print(ref_name,bam_pileup_dict.get(ref_name))

    bamfile.close()




    #Now have nested dict of {ref_tag_name: OrderedDict(position: coverage)} where OrderedDict stores positions in ascending order!

    #sj_seq_tuple to microexon coords dict
    {('chrX:81118759-81121534', 'GCTGCAGGTCAAGGTGATATGAGGCAGGAG'): 'chrX_-_81119787_81119817', ('chrX:81118780-81121534', 'GCTGCAGGTCAAGGTGATATGAGGCAGGAG'): 'chrX_-_81119787_81119817', ('chrX:92201455+92387734', 'CGGAAATCTGAAGGGAAAGTGGCAGGAAAG'): 'chrX_+_92263113_92263143', ('chrX:92387933+92618263', 'CTTTCATACCTGGACTAAAGAAAG'): 'chrX_+_92468298_92468322', ('chrY:13869128+13871632', 'ACTTAAGGAGCATGAG'): 'chrY_+_13869661_13869677', ('chrY:5501255+5737271', 'CTTTCATACCTGGACTAAAGAAAG'): 'chrY_+_5581774_5581798', ('chrY:57210792+57211760', 'TGAGAGCCACGAGCCAAG'): 'chrY_+_57211551_57211569'}



    '''
