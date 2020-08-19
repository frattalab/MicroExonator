import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pybedtools
from pybedtools import BedTool
import pyBigWig
from collections import defaultdict



Genome = {}


def percent (c, total):
	try:
		return (100*float(c))/float(total)
	except ZeroDivisionError:
		return 0

def Genomictabulator(fasta):

	#print >> sys.stderr, "Cargando genoma en la memoria RAM ...",

	f = open(fasta)

	for chrfa in SeqIO.parse(f, "fasta"):
		Genome[chrfa.id] = chrfa.seq

	#print >> sys.stderr, "OK"

	f.close()



def PWM_to_dict(file):
	reader = csv.reader(open(file), delimiter = '\t')
	header = reader.next()
	header_dict = {}
	col = 0

	matrix = {}

	for name in header:
		header_dict[name] = col
		col += 1

	A_frec = []
	C_frec = []
	G_frec = []
	T_frec = []
	N_freq = []

	for row in reader:
		A = row[header_dict["A"]]
		C = row[header_dict["C"]]
		G = row[header_dict["G"]]
		T = row[header_dict["T"]]

		A_frec.append(float(A))
		C_frec.append(float(C))
		G_frec.append(float(G))
		T_frec.append(float(T))
		N_freq.append(0)

	matrix["A"] = A_frec
	matrix["C"] = C_frec
	matrix["G"] = G_frec
	matrix["T"] = T_frec
	matrix["N"] = N_freq

	return matrix


def main(ME_centric, bed12, U2_GTAG_5_file, U2_GTAG_3_file, phylop, ME_len, ME_DB=False, out_canonical, out_db):




	n = 100

	min_intron_lenght = 80


	if phylop!="NA":
		phylop_bw = pyBigWig.open(phylop)


	U2_GTAG_5 = PWM_to_dict(U2_GTAG_5_file)
	U2_GTAG_3 = PWM_to_dict(U2_GTAG_3_file)

	U2_GTAG_5_max_score = 0
	U2_GTAG_3_max_score = 0

	for index in range(13):
		U2_GTAG_5_max_score += max(U2_GTAG_5['A'][index], U2_GTAG_5['C'][index], U2_GTAG_5['T'][index], U2_GTAG_5['G'][index])

	for index in range(17):
		U2_GTAG_3_max_score += max(U2_GTAG_3['A'][index], U2_GTAG_3['C'][index], U2_GTAG_3['T'][index], U2_GTAG_3['G'][index])

	TOTAL_U2_max_score = U2_GTAG_5_max_score + U2_GTAG_3_max_score



	found_ME = set([])
	ME_chroms = set([])

	for row in csv.reader(open(ME_centric), delimiter = '\t'):

		ME, transcript, sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, U2_scores, mean_conservations, P_MEs, total_ME = row

		ME_strand, ME_start, ME_end = ME.split("_")[-3:]
		ME_chrom =  "_".join(ME.split("_")[:-3])

		found_ME.add(ME)
		ME_chroms.add(ME_chrom)



	introns = set([])

	non_detected_ME = defaultdict(list) # a microexon can be derived from more than one transcript. The idea is to collapese the transcript



	SJ_start_seqs = {}
	SJ_end_seqs = {}


	for row in csv.reader(open(bed12), delimiter = '\t'):




		blocksizes = list(map(int, row[10].strip(",").split(",")))
		qstarts = list(map(int, row[11].strip(",").split(",")))


		start = int(row[1])
		end = int(row[2])
		strand = row[5]
		bn = int(row[9])
		chrom = row[0]
		transcript = row[3]

		f_seq = ""
		r_seq= ""

		if chrom in Genome:



			for q1, q2, b1 in zip(qstarts, qstarts[1:], blocksizes):

				istart = start + q1 + b1
				iend = start + q2

				SJ_ID =  transcript + str(istart)


				intron = " ".join([chrom, str(istart), str(iend), "SJ", "0", strand])

				#if chrom in ME_chroms:

				introns.add(intron)


				# Indexing tag library


				estart = start + q1
				eend = start + q1 + b1

				f_seq += str(Genome[chrom][estart:eend])


				if (chrom, eend) in SJ_start_seqs:

					if f_seq[-100:] > len(SJ_start_seqs[(chrom, eend )]):

						SJ_start_seqs[(chrom, eend )] = f_seq[-100:]

				else:

					SJ_start_seqs[(chrom, eend )] = f_seq[-100:]




			for q1, b1 in zip(qstarts[::-1],  blocksizes[::-1]):

				estart = start + q1
				eend =  start + q1 + b1

				r_seq = str(Genome[chrom][estart:eend]) + r_seq


				if (chrom, estart) in SJ_end_seqs:

					if r_seq[:100] > len(SJ_end_seqs[(chrom, estart )]):

						SJ_end_seqs[(chrom, estart )] = r_seq[:100]

				else:

					SJ_end_seqs[(chrom, estart )] = r_seq[:100]






			for q1, q2, q3, b1, b2, b3 in zip(qstarts, qstarts[1:] , qstarts[2:], blocksizes, blocksizes[1:], blocksizes[2:]):

				estart = start + q2
				eend = start + q2 + b2
				elength = eend - estart
				exon = "_".join([chrom, strand, str(estart),  str(eend)])

				SJ_start = start + q1 + b1
				SJ_end = start + q3
				ME_intron = " ".join([chrom, str(SJ_start), str(SJ_end), "SJ", "0", strand])



				dn = Genome[chrom][(estart-2):estart] + Genome[chrom][eend:(eend+2)]

				if strand=="-":
					dn = dn.reverse_complement()

				dn = str(dn).upper()



				if elength <= ME_len and dn=="AGGT" and exon not in found_ME:

					#if chrom in ME_chroms:

					introns.add(ME_intron)

					non_detected_ME[(chrom, estart, eend, strand, elength)].append(transcript)




	##### Microexon database ######


	if ME_DB!=False:


		for row in csv.reader(open(ME_DB), delimiter = '\t'):


			if len(row)==12:


				blocksizes = list(map(int, row[10].strip(",").split(",")))
				qstarts = list(map(int, row[11].strip(",").split(",")))


				start = int(row[1])
				end = int(row[2])
				strand = row[5]
				bn = int(row[9])
				chrom = row[0]

				if chrom in Genome:


					for q1, q2, q3, b1, b2, b3 in zip(qstarts, qstarts[1:] , qstarts[2:], blocksizes, blocksizes[1:], blocksizes[2:]):


						estart = start + q2
						eend = start + q2 + b2
						elength = eend - estart
						exon = "_".join([chrom, strand, str(estart), str(eend)])
						transcript = row[3]

						SJ_start = start + q1 + b1
						SJ_end = start + q3
						ME_intron = " ".join([chrom, str(SJ_start), str(SJ_end), "SJ", "0", strand])



						dn = Genome[chrom][(estart-2):estart] + Genome[chrom][eend:(eend+2)]

						if strand=="-":
							dn = dn.reverse_complement()

						dn = str(dn).upper()





						if elength <= ME_len and dn=="AGGT" and exon not in found_ME:

							#introns.add(ME_intron)

							non_detected_ME[(chrom, estart, eend, strand, elength)].append(transcript)



	introns_str =  "\n".join(list(introns))


	intron_bed = BedTool(introns_str , from_string=True)
	intron_bed = intron_bed.sort()


	TOTAL_SJ_starts = set([])
	TOTAL_SJ_ends = set([])


	with open(out_canonical, 'w') as out_tags, open(out_db, 'w') as out_ME_centric :



		for i in non_detected_ME.items():

			ME_info, transcripts = i

			chrom, estart, eend, strand, elength = ME_info

			transcript = transcripts[0]


			#ME = "_".join([chrom, str(estart), strand, str(eend)])
			ME = "_".join([chrom, strand, str(estart),  str(eend)])


			if elength <= ME_len and dn=="AGGT" and exon not in found_ME:


				if phylop=="NA":

					mean_conservation=0

				else:

					try:
						mean_conservation= phylop_bw.stats(chrom, estart-2, eend+2, type="mean")[0]
					except RuntimeError:
						mean_conservation=0

					if mean_conservation==None:
						mean_conservation=0


				ME5 = str(Genome[chrom][estart-14:estart+3]).upper()
				ME3 = str(Genome[chrom][eend-3:eend+10]).upper()

				micro_exon_seq_found = str(Genome[chrom][estart:eend]).upper()


				if strand == "-":

					ME5 = str(Genome[chrom][eend-3:eend+14].reverse_complement()).upper()
					ME3 = str(Genome[chrom][estart-10:estart+3].reverse_complement()).upper()

					micro_exon_seq_found = str(Genome[chrom][estart:eend].reverse_complement()).upper()



				U2_score = 0

				i = 0

				for N in ME5:
					U2_score += U2_GTAG_3[N][i]
					i += 1

				i = 0

				for N in ME3:
					U2_score += U2_GTAG_5[N][i]
					i += 1

				U2_score = percent(U2_score, TOTAL_U2_max_score)



				ME_bed = BedTool(" ".join([chrom, str(estart), str(eend - 1), "ME", "0", strand]) , from_string=True)

				SJs_bed = intron_bed.intersect(ME_bed, wa=True, s=True, F=1, nonamecheck=True)


				SJs = set([])

				SJ_starts = []
				SJ_ends = []



				if len(SJs_bed)!=0:


					for sj in SJs_bed:

						SJ_chrom, SJ_start, SJ_end, ID, score, SJ_strand = str(sj).strip("\n").split("\t")
						SJ = SJ_chrom + ":" + SJ_start + SJ_strand + SJ_end

						SJ_starts.append(int(SJ_start))
						SJ_ends.append(int(SJ_end))

						SJs.add(SJ)

						TOTAL_SJ_starts.add((chrom, SJ_start))
						TOTAL_SJ_ends.add((chrom, SJ_end))


						### TAG creation


						UP_TAG =  SJ_start_seqs[(SJ_chrom, int(SJ_start) )]
						DOWN_TAG =  SJ_end_seqs[(SJ_chrom, int(SJ_end) )]


						ME_TAG = UP_TAG +  Genome[chrom][estart:eend] + DOWN_TAG

						tag_pos = "_".join(map(str, [len(UP_TAG), micro_exon_seq_found, len(DOWN_TAG)]))



						if strand == "-":

							 ME_TAG = ME_TAG.reverse_complement()

							 tag_pos = "_".join(map(str, [len(UP_TAG), micro_exon_seq_found, len(DOWN_TAG)][::-1]))

						ME_TAG = str(ME_TAG).upper()



						ME_TAG_ID = chrom + ":" + "".join([ str(estart), strand, str(eend)])



						out_tags.write(">" + "|".join([ SJ, transcript ,  tag_pos]) + "\n" )
						out_tags.write(ME_TAG + "\n")

						# print ">" + "|".join([ ME_TAG_ID, transcript,  tag_pos ])
						# print ME_TAG





					total_SJs = ",".join(SJs)



					min_intron_seq = str(Genome[chrom][max(SJ_starts):min(SJ_ends)]).upper()

					if strand == "-":

						min_intron_seq = str(Genome[chrom][max(SJ_starts):min(SJ_ends)].reverse_complement()).upper()



					total_number_of_micro_exons_matches = min_intron_seq.count("AG" + micro_exon_seq_found + "GT")


					P_ME = 	1 - ( 1 - (float(1)/float(4**len(micro_exon_seq_found)+4)))**( len(min_intron_seq) - (len(micro_exon_seq_found)+4))


					info = ME, transcript, 0, total_SJs, 0, elength, micro_exon_seq_found, total_number_of_micro_exons_matches, U2_score, mean_conservation, P_ME, "|".join(map(str,  [ME, U2_score, mean_conservation]))


					out_ME_centric.write("\t".join(map(str, info)) + "\n")












#ME, transcript, sum_total_coverage, total_SJs, total_coverages, len_micro_exon_seq_found, micro_exon_seq_found, total_number_of_micro_exons_matches, U2_scores, mean_conservations, P_MEs, total_ME = row


if __name__ == '__main__':
	Genomictabulator(sys.argv[1])
	main(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], int(sys.argv[7]), sys.argv[8], sys.argv[9], sys.argv[10])


#python2 ~/my_src/Micro-Exonator/Get_annotated_microexons.py ../../../../../Genome/mm10/mm10.fa Round1/TOTAL/TOTAL.sam.row_ME.filter1.ME_centric ../../../../../Genome/mm10/Tracks/Gene_annotation/gencode.vM11.annotation.bed12 ../../../../../Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_5.good.matrix ../../../../../Genome/mm10/Tracks/SpliceRack/mm10_GT_AG_U2_3.good.matrix ../../../../../Genome/mm10/Tracks/Phylop/mm10.60way.phyloP60way.bw 30
