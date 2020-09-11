
import sys
import csv


#cryptics run - get error _csv.Error: field larger than field limit (131072) - this sets limit to maximum integer value variable can take
csv.field_size_limit(sys.maxsize)

def main(ME_centric):

    header = ["ME", "U2_score", "Vertebrate_conservation", "ME_len", "ME_max_U2"]

    print("\t".join(header))

    for row in csv.reader(open(ME_centric), delimiter = '\t'):

        ME, transcript, sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, U2_scores, mean_conservations, P_MEs, total_ME = row

        ME_strand, ME_start, ME_end = ME.split("_")[-3:]
        ME_chrom =  "_".join(ME.split("_")[:-3])

        #print(total_ME.split("|"))

        for ME_match in total_ME.split(","):

            print("\t".join(ME_match.split("|") + [len_micro_exon_seq_found, U2_scores]))

if __name__ == '__main__':
        main (sys.argv[1])
