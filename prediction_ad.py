# /usr/bin/python
from Bio import SeqIO
import os
import sys
"""
if len(sys.argv)<2:
	sys.exit("[Usage:\n\t python <scriptsname> <protein_fasta> <score> <output>")

if os.path.exists(sys.argv[1]):
	protein_fasta = sys.argv[1]
else:
	sys.exit("%s doesn't exist" % sys.argv[1])

threshold = int(sys.argv[2])

if os.path.exists(sys.argv[3]):
	sys.exit("%s exists already, please rename it" % sys.argv[3])
else:
	output = sys.argv[3]
"""

def fun(protein_fasta, threshold, ad_predict_1):
	# motif pattern for the prediction of H2-Ad ligands in four positions 
	mp_H2_Ad_4 = [['A',1],['C',1],['D',1],['E',3],['F',1],['G',-2],\
	['H',1],['I',1],['K',-2],['L',2],['M',2],['N',-2],['P',-2],['Q',0],\
	['R',-2],['S',4],['T',4],['V',1],['W',1],['X',0],['Y',4]]

	mp_H2_Ad_7 = [['A',8],['G',-3],['H',-2],['I',8],['L',8],['V',10]]

	mp_H2_Ad_9 = [['A',10],['G',-3],['K',-3],['L',6],['R',-3],['V',8]]

	mp_H2_Ad_12 = [['A',10],['G',6],['I',6],['K',-3],['L',6],['S',10],['V',6]]

	# define value for non-map 
	default = 0
	# the length of predicted peptide is 15
	fifteen = 15

	# building map
	map_H2_Ad_4 = dict(mp_H2_Ad_4)
	map_H2_Ad_7 = dict(mp_H2_Ad_7)
	map_H2_Ad_9 = dict(mp_H2_Ad_9)
	map_H2_Ad_12 = dict(mp_H2_Ad_12)

	#creat the output file
	fout = open(ad_predict_1,'w')

	# open the protein file
	fin = open(protein_fasta,'r')

	# give each peptide a id number
	num = 1
	score = 0
	for seq_record in SeqIO.parse(fin,"fasta"):
		#print seq_record.id
		sequences = seq_record.seq 
		for i in range(0, len(seq_record)-14):
			# peptide, 15 amino acids
			aa15 = sequences[i:i+fifteen]
			# specific position scores
			pos4 = map_H2_Ad_4.get(aa15[3],default)
			pos7 = map_H2_Ad_7.get(aa15[6],default)
			pos9 = map_H2_Ad_9.get(aa15[8],default)
			pos12 = map_H2_Ad_12.get(aa15[11],default)
			# total scores
			score = pos4+pos7+pos9+pos12	
			# selecte expected score
			if score >= threshold:
				title_line = ">%s_%d %s\n" % (seq_record.description,num,score)
				fout.write(title_line)
				aa_line = "%s\n" % aa15
				fout.write(aa_line)
				num+=1
			score = 0

"""
ok = 1 
protein = ""
for line in open(protein_fasta,'r'):
	if line[0]=='>':
		if ok==0:
			list_aa15 = []
			for i in range (0,len(protein)-15):
			        aa15 = protein[i:i+15]
        			pos4 = map_H2_Ad_4.get(aa15[3],default)
			        pos7 = map_H2_Ad_7.get(aa15[6],default)
			        pos9 = map_H2_Ad_9.get(aa15[8],default)
			        pos12 = map_H2_Ad_12.get(aa15[11],default)
			        score = pos4+pos7+pos9+pos12
				if score >= threshold:
			        	list_aa15.append([score,aa15,i+1])
			if not list_aa15:
				list_aa15.sort(reverse=True)
				for peptide in list_aa15:
					out_line = "%s\t%d\t%d\n" % (peptide[1],peptide[0],peptide[2])
					fout.write(out_line)
				fout.write(line)
				protein = ""
		else:
			fout.write(line)
			ok = 0
	else:
		protein = protein + line.strip()
		

list_aa15 = []
for i in range (0,len(protein)-15):
	aa15 = protein[i:i+15]
	pos4 = map_H2_Ad_4.get(aa15[3],default)
	pos7 = map_H2_Ad_7.get(aa15[6],default)
	pos9 = map_H2_Ad_9.get(aa15[8],default)
	pos12 = map_H2_Ad_12.get(aa15[11],default)
	score = pos4+pos7+pos9+pos12
	list_aa15.append([score,aa15,i+1])

list_aa15.sort(reverse=True)
for peptide in list_aa15:
	out_line = "%s\t%d\t%d\n" % (peptide[1],peptide[0],peptide[2])
	fout.write(out_line)
"""




