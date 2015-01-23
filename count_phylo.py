import os
import sys, getopt
import csv
from collections import Counter
import re


def main(argv):
	inputfile = ''
	outputfile = ''
	
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
		print 'write in the form:script.py -i <infilfe> -o <outfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'write in the form:script.py -i <infilfe> -o <outfile>'
			sys.exit()
		elif opt in ("-i", "--infile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg
			fo = open(outputfile,"rw+")
	if FileCheck(inputfile):
		lol= list(csv.reader(open(inputfile,'rb'),delimiter='\t'))
		phylo_count = list()
		taxas=[]
		#for record in lol[-1][:]:
		taxonomy = [sub_list[-1] for sub_list in lol]
		for tax in taxonomy:
			
			record=re.split(r'\[\'|\'\, \'|\'\]', tax)
			taxas.append(record[3])
		taxa_dict=Counter(taxas)
		
		for x in taxa_dict:
			if outputfile != '':
				fo.write(x+"\t"+str(taxa_dict[x]))
			else:
				print x+"\t"+str(taxa_dict[x])
		
		


def FileCheck(fn):
	try:
		open(fn,"r")
		return 1
	except IOError:
		print "Error: File does not exist or cannot be opened."
		return 0
if __name__ == "__main__":
	main(sys.argv[1:])
		