# /usr/bin/python
from Bio import SeqIO
import os
import sys

def fun(protein_fasta, threshold, ad_predict_1):
	# motif pattern for the prediction of H2-Ed ligands in five positions
	mp_H2_Ed_4 = [['F',10],['I',6],['L',6],['M',2],['V',6],['W',10],['Y',10]]
	mp_H2_Ed_7 = [['I',4],['K',6],['R',6]]
	mp_H2_Ed_9 = [['G',2],['I',2],['L',2],['V',2]]
	mp_H2_Ed_11 = [['K',8],['R',8]]
	mp_H2_Ed_12 = [['K',8],['R',8]]

	# define value for non-map 
	default = 0
	# the length of predicted peptide is 15
	fifteen = 15

	# give each peptide a id number
	num = 1

	# building map
	map_H2_Ed_4 = dict(mp_H2_Ed_4)
	map_H2_Ed_7 = dict(mp_H2_Ed_7)
	map_H2_Ed_9 = dict(mp_H2_Ed_9)
	map_H2_Ed_11 = dict(mp_H2_Ed_11)
	map_H2_Ed_12 = dict(mp_H2_Ed_12)

	#creat the output file
	fout = open(ad_predict_1,"w")

	# open the protein file
	fin = open(protein_fasta,"r")
	score = 0
	for seq_record in SeqIO.parse(fin,"fasta"):
		#print seq_record.id
		sequences = seq_record.seq 
		for i in range(0, len(seq_record)-14):
			# peptide, 15 amino acids
			aa15 = sequences[i:i+fifteen]
			# specific position scores
			pos4 = map_H2_Ed_4.get(aa15[3],default)
			pos7 = map_H2_Ed_7.get(aa15[6],default)
			pos9 = map_H2_Ed_9.get(aa15[8],default)
			pos11 = map_H2_Ed_11.get(aa15[10],default)
			pos12 = map_H2_Ed_12.get(aa15[11],default)
			# total scores
			score = pos4+pos7+pos9+pos11+pos12	
			# selecte expected score
			#print score
			if score >= threshold:
				title_line = ">%s_%d %s\n" % (seq_record.description,num, score)
				fout.write(title_line)
				aa_line = "%s\n" % aa15
				fout.write(aa_line)
				num+=1
			score = 0

	#fout.close()
	#fin.close()
