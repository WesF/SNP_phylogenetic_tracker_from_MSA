#! /usr/bin/env python

	# Version 1.75 - Return an ordered list of sequences
	# Version 1.7  - Return multiple matches in tie amongh sequence lengths
	# Version 1.61 - Take user input in addition to file name
	# Version 1.5  - Split on X as well as * for size tests
	# Version 1.41 - Minor bug fix
	# Version 1.4  - Added checkfile option for faster connection
	# Version 1.3  - Added min length of ORF to retain. 
	#                Default is 1. Use 0 to print all
	# Version 1.2  - Reports more of the ambiguity info
	# Version 1.1  - Requires the mybio module for some functions
	# 			   www.mbari.org/staff/haddock/scripts/


usage="""
TRANSLATEDNA.PY - version 1.6

This program takes a fasta file of nucleotide sequences,
does a six-frame translation, and tries to determine the 
longest open reading frame from each sequence. 

If it finds one, it prints that out in fasta format. If 
it finds two long ORFs, the it prints both, appending the 
number of the ORF that was selected. 

The criteria used to select the best ORFs are (1) length of
the longest ORF and (2) the number of individual ORFs. 

PRINT BEST ORF 0/[1] 
To print all the ORFs and not just the best, add an argument 
of 0 after the file name. (Adding 1 will print best match.)

CHECK FILE [0]/1
For large (clean) sequence files you can run w/o checking 
dupes or file names by adding 3rd parameter of zero [no check]
(this is the default) or one for checkfile. It goes about 150 
times faster w/o checking.

The names are appended with the ORF#, or in the case of a tie
that is followed by _x and the number of additional equally
valid translations

To save to a file, put this after the command name:
 	> outputname.txt 

You can also input a sequence (w/o spaces or returns) instead
of a file name and it will translate that to the screen

Bug reports  - beroe [at] mac {dot} com

Requires: mybio.py mini-library at www.mbari.org/~haddock/scripts/

Usage: 
	translatedna.py filename.fasta [bestORFonly] [checkfile] 
	

"""

import re
import sys
import os
# needed for maketrans -- not sure if this is obsolete
import string
import mybio as mb
standardWithAmbig = {
	'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAR':'K', 'AAT':'N', 'AAY':'N', 'ACA':'T', 'ACB':'T', 
	'ACC':'T', 'ACD':'T', 'ACG':'T', 'ACH':'T', 'ACK':'T', 'ACM':'T', 'ACN':'T', 'ACR':'T', 
	'ACS':'T', 'ACT':'T', 'ACV':'T', 'ACW':'T', 'ACY':'T', 'AGA':'R', 'AGC':'S', 'AGG':'R', 
	'AGR':'R', 'AGT':'S', 'AGY':'S', 'ATA':'I', 'ATC':'I', 'ATG':'M', 'ATH':'I', 'ATM':'I', 
	'ATT':'I', 'ATW':'I', 'ATY':'I', 'CAA':'Q', 'CAC':'H', 'CAG':'Q', 'CAR':'Q', 'CAT':'H', 
	'CAY':'H', 'CCA':'P', 'CCB':'P', 'CCC':'P', 'CCD':'P', 'CCG':'P', 'CCH':'P', 'CCK':'P', 
	'CCM':'P', 'CCN':'P', 'CCR':'P', 'CCS':'P', 'CCT':'P', 'CCV':'P', 'CCW':'P', 'CCY':'P', 
	'CGA':'R', 'CGB':'R', 'CGC':'R', 'CGD':'R', 'CGG':'R', 'CGH':'R', 'CGK':'R', 'CGM':'R', 
	'CGN':'R', 'CGR':'R', 'CGS':'R', 'CGT':'R', 'CGV':'R', 'CGW':'R', 'CGY':'R', 'CTA':'L', 
	'CTB':'L', 'CTC':'L', 'CTD':'L', 'CTG':'L', 'CTH':'L', 'CTK':'L', 'CTM':'L', 'CTN':'L', 
	'CTR':'L', 'CTS':'L', 'CTT':'L', 'CTV':'L', 'CTW':'L', 'CTY':'L', 'GAA':'E', 'GAC':'D', 
	'GAG':'E', 'GAR':'E', 'GAT':'D', 'GAY':'D', 'GCA':'A', 'GCB':'A', 'GCC':'A', 'GCD':'A', 
	'GCG':'A', 'GCH':'A', 'GCK':'A', 'GCM':'A', 'GCN':'A', 'GCR':'A', 'GCS':'A', 'GCT':'A', 
	'GCV':'A', 'GCW':'A', 'GCY':'A', 'GGA':'G', 'GGB':'G', 'GGC':'G', 'GGD':'G', 'GGG':'G', 
	'GGH':'G', 'GGK':'G', 'GGM':'G', 'GGN':'G', 'GGR':'G', 'GGS':'G', 'GGT':'G', 'GGV':'G', 
	'GGW':'G', 'GGY':'G', 'GTA':'V', 'GTB':'V', 'GTC':'V', 'GTD':'V', 'GTG':'V', 'GTH':'V', 
	'GTK':'V', 'GTM':'V', 'GTN':'V', 'GTR':'V', 'GTS':'V', 'GTT':'V', 'GTV':'V', 'GTW':'V', 
	'GTY':'V', 'MGA':'R', 'MGG':'R', 'MGR':'R', 'TAA':'*', 'TAC':'Y', 'TAG':'*', 'TAR':'*', 
	'TAT':'Y', 'TAY':'Y', 'TCA':'S', 'TCB':'S', 'TCC':'S', 'TCD':'S', 'TCG':'S', 'TCH':'S', 
	'TCK':'S', 'TCM':'S', 'TCN':'S', 'TCR':'S', 'TCS':'S', 'TCT':'S', 'TCV':'S', 'TCW':'S', 
	'TCY':'S', 'TGA':'*', 'TGC':'C', 'TGG':'W', 'TGT':'C', 'TGY':'C', 'TRA':'*', 'TTA':'L', 
	'TTC':'F', 'TTG':'L', 'TTR':'L', 'TTT':'F', 'TTY':'F', 'YTA':'L', 'YTG':'L', 'YTR':'L',
	'---': '-', '...': '-', '~~~': '-'
}

def revcomp(dna):
	#reverse complement of a DNA sequence
	comp = dna.translate(string.maketrans("ATGCatgcRYMKrymkHBVDhbvd", "TACGTACGYRKMYRKMDVBHDVBH"))

	lcomp = list(comp)
	lcomp.reverse()
	return "".join(lcomp)

def dna_translate(cdna, code = standardWithAmbig):
	# translate a cDNA sequence to a protein
	allprot=[]
	allprot.append("".join([ code.get(cdna[i:i+3], "X") for i in xrange(0,len(cdna),3) ]) )
	allprot.append("".join([ code.get(cdna[i:i+3], "X") for i in xrange(1,len(cdna),3) ]) )
	allprot.append("".join([ code.get(cdna[i:i+3], "X") for i in xrange(2,len(cdna),3) ]) )
	revcdna = revcomp(cdna) 
	allprot.append("".join([ code.get(revcdna[i:i+3], "X") for i in xrange(0,len(revcdna),3) ]) )
	allprot.append("".join([ code.get(revcdna[i:i+3], "X") for i in xrange(1,len(revcdna),3) ]) )
	allprot.append("".join([ code.get(revcdna[i:i+3], "X") for i in xrange(2,len(revcdna),3) ]) )
	return allprot

def findlongest(sequencelist):
	# finds the longest continuous ORF, and also the one with the fewest pieces
	# could add a metric for top two ORFs combined, or count the X's as well
	
	lenscores=[]
	fragscores=[]
	for seq in sequencelist:
		orfs=seq.rstrip().replace('X','*').split('*')
		# print orfs
		lenlist= [len(frag) for frag in orfs]
#		print lenlist
		lenscores.append(max(lenlist))
		fragscores.append(len(lenlist)-lenlist.count(0))
	# print fragscores
	# print lenscores

	maxlenlist = [i for i,x in enumerate(lenscores) if x==max(lenscores)]
#	maxlenind=lenscores.index(max(lenscores)) # max of the stored fragment lengths
	# print maxlenind
	minfraglist = [i for i,x in enumerate(fragscores) if x==min(fragscores)]

#	minfragind=fragscores.index(min(fragscores)) # min of the stored number of fragments
#	print maxlenlist,minfraglist
	# print lenscores, maxlenind
	# print fragscores, minfragind
	# make a master list with all the "winners"
	masterlist= list(set(minfraglist + maxlenlist))
	# print "masterlist",masterlist
	if (len(masterlist)==1):
		return [ [sequencelist[masterlist[0]]] , masterlist ]
		
	else:
		return [[sequencelist[x] for x in masterlist],masterlist]# return a count of 1 if just one each mismatch, 
		
# EMD FUNCTION



##########################
# START MAIN PROGRAM
skipcount=0
checkfile=0
minsize=1
blockwidth=71
noread=False
sequences={}
ordernames=[]
	
if len(sys.argv)<2:
	sys.stderr.write(usage)

else:
	# second parameter for printing all seqs
	# use zero to print all. otherwise it will be a minimum size
	# default is 1 (no size limit)
	if len(sys.argv)>2:
		minsize=int(sys.argv[2])
	if len(sys.argv)>3:
		checkfile=int(sys.argv[3])
		
	inputfilename = sys.argv[1]
	
	if not(os.path.exists(inputfilename)):
		sequences['userinput']=inputfilename.replace('\n','').replace(' ','').replace('\t','') # make a fake name
		noread=True  # Take the input from the command line instead of the file
		sys.stderr.write("\n# No file found, translating input as sequence:\n\n")
		ordernames=['userinput']
	
	else:
		# This function is in the mybio module
		if checkfile:
			sequences,ordernames = mb.readfasta(inputfilename,False,True,True) # don't append sequences, but *do* remove dupes, and in order
		else:
			sequences,ordernames = mb.readfasta_nocheck(inputfilename,False, False, False,True) # don't append sequences, but *do* remove dupes, and in order


#	for name in sequences.keys():
	for name in ordernames:
		sequence=sequences[name]
		newseq = re.sub(r'[\n\s]','',sequence).upper()
		#newseq = sequence.replace('\n','').replace(' ','').replace('\t','')

		# print newseq

		# returns all 6-frames of translation as a list
		seqlist=dna_translate(newseq)
		
		# print all six frames
		if minsize==0:
			for myseq in seqlist:
				print ">"+name[:20]+ "_"+str(seqlist.index(myseq))
				mb.printblock(myseq)
	
		else:
			# finds the lonest sequence and what frame it's in
			[longorf,frame] = findlongest(seqlist)
			# print "longorf,frame",longorf,frame
			
			for ind  in range(len(longorf)):
				print ">" + name + "_" + str(frame[ind])
				mb.printblock(longorf[ind],blockwidth)
				
		# print seqlist[0].split('*')
	
