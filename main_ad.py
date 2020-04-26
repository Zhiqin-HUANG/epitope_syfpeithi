#! /usr/bin/python 
import prediction_ad
import os
import sys
import time
from Bio.Blast.Applications import BlastallCommandline
from Bio.Blast import NCBIXML
from Bio.Blast import Record
from Bio import SeqIO

start = time.time()

# check arguments
if len(sys.argv) < 3:
	sys.stderr.write('Usage: python [scriptsname] [inputfasta.file] [output directory] [threshold score]\n')
	sys.exit(1)

if os.path.exists(sys.argv[1]):
	input_fasta = sys.argv[1]
else:
	sys.stderr.write('ERROR: Fasta file %s was not found!\n' % sys.argv[1])
	sys.exit(1)

if os.path.exists(sys.argv[2]):
	outputdir = sys.argv[2]
	outputdir = os.path.dirname(outputdir)
else:
	sys.stderr.write('ERROR: Output directory %s was not found!\n' % sys.argv[2])
	sys.exit(1)

expected_value = int(sys.argv[3])

#####################
# input your parameters
# expected_value = 20
# lito_genome database
litoDB = "/home/zhiqin/vaccine/db/Litomosoides_sigmodontis_genomeDB"
# mouse database
mouseDB = "/home/zhiqin/vaccine/db/est_mouse"
#####################

# make the result file 
resultPath = outputdir+"/results/"
if not os.path.isdir(resultPath):
	os.mkdir(resultPath)

# make the ad_predicted fasta for the original proteins
ad_1 = outputdir + "/ad_predict_1"

# SYFPEITHY epitope prediction for protein candidates
prediction_ad.fun(input_fasta,expected_value,ad_1)

# blast fasta to lito_database
blast_xml = ad_1+"_lito_blast.xml"
command = "blastall -i %s -d %s -p tblastn -m 7 -o %s" % (ad_1,litoDB,blast_xml)
os.system(command)

print "After blasting against lito(s): %d " % (time.time()-start)

# using biopython package to parse blast_xml file
result_handle = open(blast_xml,'r')
blast_records = NCBIXML.parse(result_handle)

# create a new fasta file for store all the 15aa without predicted score
fpep = outputdir + "/peptide15aa.fasta"
fpeptide = open(fpep,'w')

# parse the blast.xml, output a new fasta
for blast_record in blast_records:
	#print blast_record.query
	#print blast_record.num_hits
	#print dir(blast_record)
	#print dir(blast_record.alignments)
	for alignment in blast_record.alignments:
		for hsp in alignment.hsps:
			if len(hsp.sbjct)==15: # only the 15aa are selected
				if hsp.identities!=15:
					fpeptide.write(">"+blast_record.query+" " + hsp.query+ "\n")
				else:
					fpeptide.write(">"+blast_record.query+"\n")
				fpeptide.write(hsp.sbjct+"\n")
fpeptide.close()
	
# make the ad_predicted fasta for the original proteins
ad_2 = outputdir + "/ad_predict_2"

# SYFPEITHY epitope prediction for seletect 15aa candidates
prediction_ad.fun(fpep,expected_value,ad_2)

# blast fasta to mouse genome database
blast_xml_2 = ad_2+"_mouse_blast.xml"
#command = "blastall -i %s -d %s -p tblastn -m 7 -o %s" % (ad_2,mouseDB,blast_xml_2) #bug
#os.system(command)
cline = BlastallCommandline(cmd="blastall", program="tblastn", infile=ad_2, database=mouseDB, expectation=0.001, outfile=blast_xml_2, align_view='7')
#print cline
command = str(cline)
os.system(command)

# using biopython package to parse blast_xml file
result_handle = open(blast_xml_2,'r')
blast_records2 = NCBIXML.parse(result_handle)

# parse the blast.xml, output list of peptide-titles, which does present in mouse genome.
peptide_in_mouse = []
for blast_record in blast_records2:
	if blast_record.alignments:
		#for hsp in alignment.hsps:
		peptide_in_mouse.append(blast_record.query)

# delete repeative elements
#print len(peptide_in_mouse)
#peptide_in_mouse = list(set(peptide_in_mouse))
#print len(peptide_in_mouse)

# create the final ad_prediction file
ad_3 = outputdir + "/ad_predicted_unsorted.fasta"
f_final = open(ad_3,'w')


a = 0
b = 0
c = 0
for seq_record in SeqIO.parse(open(ad_2,'r'),"fasta"):
	a+=1
	if seq_record.description not in peptide_in_mouse:
		#print seq_record.description
		b+=1
		buff = ">"+ seq_record.description + "\n" + str(seq_record.seq) + "\n"
		f_final.write(buff)
	else:
		c+=1
f_final.close()
#print a, b , c

# sorting the final sorted ad_prediction file
fin = open(ad_3,'r')
# create a empty list to store peptide information
lisPep = []
for seq_record in SeqIO.parse(fin,"fasta"):
    title = seq_record.description
    score = int(title.split()[-1])
    sequences = seq_record.seq
    lisPep.append([score,title,sequences])

# create the final sorted file
ad_4 = resultPath + "/ad_predicted_final.fasta"
f_ad_4 = open(ad_4,'w')
lisPep.sort(reverse=True)
for element in lisPep:
    f_ad_4.write(">"+str(element[1])+"\n"+str(element[2])+"\n")

print "Total runtime(s): %d" % (time.time()-start)

"""
# parse lito_blast files 
filelist = os.listdir(idPath)
for filename in filelist:
	if fnmatch.fnmatch(filename,"*blast"):
		parse_blast.parse(idPath+filename)

# deal with the not identity peptide compared to lito_genome
subMain_ForNotId.subMain(idPath,resultPath,expected_value,outputdir)

# join the fasta file in resultPath according motif Ad and Ed, respectively
command = "cat "+resultPath+"*Ad*"+" > "+ resultPath+"ad_fasta_total"
os.system(command)
command = "cat "+resultPath+"*Ed*"+" > "+ resultPath+"ed_fasta_total"
os.system(command)

#blast mouse genome using the total Ad and Ed fasta 
filelist = os.listdir(resultPath)
for filename in filelist:
	if "_fasta_total" in filename:
		command="blastall -i "+resultPath+filename+" -d "+mouseDB+" -p tblastn -o "+resultPath+filename+"_mouse_blast"
		print "now blast to mouse...Be patient..."
		os.system(command)
		print "Time: " + str(time.time()-start)

# parse the mouse blast file and extrac the title of peptide not found in mouse 
# join with peptide in total fasta file
integrate_title.exclude(resultPath)

# order the peptide file according to score
sorted_score.sort_score(resultPath)


# complete the query information
add_details_notId.fun(resultPath,idPath)

print "Total running time:"+ str(time.time()-start)
"""







