import os
import argparse
import os.path
import logging
import math
import xml.dom.minidom
try:
	import numpy
except ImportError,e:
	print "Numerical Python (numpy) not installed properly! This is a dependency."
	sys.exit(0)
try:
	import Bio
except ImportError,e:
	print "Biopython not installled properly! This is a dependency."
	sys.exit(0)
try:
	import matplotlib
except ImportError,e:
	print "Matplotlib not installled properly! This is a dependency."
	sys.exit(0)	
try:
	from Bio.Phylo.PAML import yn00
except ImportError,e:
	print "Could not import PAML. This is needed to calculate Dn/Ds."
	
from Bio.Blast.Applications import NcbiblastpCommandline
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Alphabet import IUPAC
from Bio import AlignIO
from Bio import Align
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio.Align.Applications import ClustalwCommandline
import subprocess as sub
import time
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import re
import csv
import ast


class Step(object):
	name = ''
	output_files = [] # override with the output files
	
	def clean(self, state):
		for f in self.output_files:
			path_to_delete = state[f]
			if os.path.exists(path_to_delete):
				logging.debug('deleting %s' % (path_to_delete))
				os.remove(path_to_delete)

	def run(self, state, input):
		raise NotImplementedError()


class MsaStep(Step):
	input = None
	output = None
	#format changed from 'CLUSTAL'
	format = 'CLUSTAL'
	def run(self, state):
		POSSIBLE_PATHS= ["/Users/Wes/Desktop/clustalo","c:\\Program Files (x86)\\ClustalW2\\clustalw2.exe"]
		selected_path = None
		for path in POSSIBLE_PATHS:
			if os.path.exists(path):
				selected_path = path
				break
		cline = ClustalwCommandline(selected_path, infile=state[self.input], outfile=state[self.output], output=self.format)
		cline()
		
class Step1(MsaStep):
	name = 'Step1'
	output_files = ['msa']
	input = 'input'
	output = output_files[0]
	format= 'FASTA'
			
			
class Step2(Step):
	name= 'Step2'
	output_files=['position_data','temp', 'blast_xml']
	
	def makedb(self,state,dbpath):
		blastDB = "makeblastdb -in "+dbpath+" -input_type fasta -dbtype prot -out "+state.blastdb#+" -logfile "+state.log
		os.system(str(blastDB))
		
		return state.blastdb
	def column_from_residue_number(self,id, res_no, aln):
		align = AlignIO.read(aln,"fasta")
		
		#TODO: handle if alignment is not found
		print(align[0].id)
		rec = next((r for r in align if r.id == id), None)
		
		#TODO:handle if rec is not found in alignment
		j = 0
		for i, res in enumerate(rec.seq.tostring()):
			if res!='-':
				if j==int(res_no):
					return i
				j+=1
	def blastp(self,state,fas, db):
		blast_in = (fas)
		blast_path = os.path.join(state.working_dir, 'blast.xml')
		blastp_cline = NcbiblastpCommandline(query=blast_in, db=db, outfmt=5, out=blast_path, max_target_seqs=1)
		os.system(str(blastp_cline))
		print blastp_cline
		return blast_path
	
	def run(self,state):
		Entrez.email = "weslfield@gmail.com"
		with file(state.msa) as input_file:
			msa_object = Bio.AlignIO.read(input_file,'fasta')
		print 'Input the positions of the alignment you want to track seperated by spaces and then press Enter.'
		cmd_positions = raw_input('-->')
		positions=cmd_positions.split()
		#add check to make sure all are integers!
		print 'Thanks! Now enter the filepath of your local BLASTX database to identify the sequences.'
		db_path = raw_input('-->')
		db=self.makedb(state, db_path)
		print 'Thanks! Now enter a reference ID to track within the MSA (leave out ">").'
		refID = raw_input('-->')
		
		
		with file(state.msa, "r") as input_file:
			
			
			adjusted_positions=[]
			for position in positions:
				to_append=self.column_from_residue_number(refID,position,state.msa)
				adjusted_positions.append(to_append)
			print adjusted_positions
			with open(state.position_data, "w") as save_file:	
				for record in SeqIO.parse(input_file,'fasta'):
					record.seq.alphabet = IUPAC.protein
					with open(os.path.join(state.working_dir, 'temp.fasta'), "w+") as temp:
						temp.write(">"+record.id+"\n")
						temp.write(str(record.seq.ungap('-'))+"\n")
						#temp.write(record.format('fasta'))
						
					blast=self.blastp(state,temp.name,db)
					with open(blast) as to_read:
						blast_record = NCBIXML.read(to_read)
						#blast_record = next(blast_records)
						#blast_record = next(blast_records)
					foo=[record.id]
						
					for position in adjusted_positions:
						foo.append(record.seq[int(position)-1])
					#for alignment in blast_record.alignments:
						#print alignment.title+"\n"
					try:
						print "The title is "+blast_record.alignments[0].title
					except IndexError:
						print "can't find the record"
					else:
						name=re.split(r'[\[\]]+', blast_record.alignments[0].title)
						print name[1]
						handle = Entrez.esearch(db="Taxonomy", term=name[1])
						try:
							record = Entrez.read(handle)
						except RuntimeError:
							continue
						else:
							pass
						try:
							taxonomy = Entrez.efetch(db="Taxonomy", id=record["IdList"][0], retmode="xml")
						except IndexError:
							print "there was an error with Entrez efetch!"
						else:
							taxrecords = Entrez.read(taxonomy)
							tax_list=taxrecords[0]["Lineage"].split('; ')
							del tax_list[0:2]
						
							new_foo=foo+tax_list

							save_file.write('\t'.join(map(str,new_foo)))
							save_file.write("\n")

					os.remove(os.path.join(state.working_dir, 'temp.fasta'))
					os.remove(blast)
					


class Step3(Step):
	name='Step3'
	output_files=['mutations_to_krona','to_krona','from_krona','mutations_from_krona']
	def run(self,state):
		lol=list(csv.reader(open(state.position_data, 'rb'), delimiter='\t'))
		print len(lol)
		DICT={}
		for i in range(len(lol[1])-2):
			
			for j in range(len(lol)):
				if i == 0:
					DICT[lol[j][0]]= [lol[j][i+1]]
				else:
					DICT[lol[j][0]].append(lol[j][i+1])
				print DICT[lol[j][0]]
		with file(state.mutations_to_krona,"w") as mutations:
			for key in DICT:
				
				to_write='1\t'+'\t'.join(DICT[key])+'\n'
				mutations.write(to_write)
		
		execute_krona="perl "+state.base_dir+"/KronaTools-2.4/scripts/ImportText.pl -n 'Amino Acid Position Frequencies' -o "+state.mutations_from_krona+" "+state.mutations_to_krona
		try:
			p = sub.check_call(['perl',state.base_dir+'/KronaTools-2.4/scripts/ImportText.pl','-n','Amino Acid Position Frequencies','-o',state.mutations_from_krona,state.mutations_to_krona],stdout=sub.PIPE,stderr=sub.PIPE)
		except sub.CalledProcessError:
			print "There was a problem running Krona. Please make sure you have installed Krona and not altered its location within this folder."
			
		with file(state.to_krona,"w") as to_krona:
			for i in range(len(lol)):
				tax_list = ast.literal_eval(lol[i][-1])
				#tax_list = re.split(r'[(\[\')(\'\, \')(\'\])]+', lol[i][-1])
				tax_line = '\t'.join(tax_list)
				to_krona.write('1\t'+tax_line+'\n')
		execute_krona="perl "+state.base_dir+"/KronaTools-2.4/scripts/ImportText.pl -n 'Phylogeny' -o "+state.from_krona+" "+state.to_krona
		
		try:
			p2 = sub.check_call(['perl',state.base_dir+'/KronaTools-2.4/scripts/ImportText.pl','-n','Phylogeny','-o',state.from_krona,state.to_krona],stdout=sub.PIPE,stderr=sub.PIPE)
		
		except sub.CalledProcessError:
			print "There was a problem running Krona. Please make sure you have installed Krona and not altered its location within this folder."


STEPS = [Step1, Step2, Step3]



		
			
