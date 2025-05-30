import sys, argparse
# this is the function for running blast local using Python for RNA
# percent identity is set that it finds hits with minimum of 15bp off-targets
parser = argparse.ArgumentParser(description='get input and output filenames')
inputfile = ''
outputfile = ''
database = ''
parser.add_argument('-i','--inputfile',help='input file name')
parser.add_argument('-o','--outputfile',help='output file name')
parser.add_argument('-db','--database',help='database being BLASTed')
parser.add_argument('-plm','--probelengthmin',type=float,help='min length of probe')
parser.add_argument('-mhs','--homologysizemin',type=float,help='min homology size for blast hits')
parser.add_argument('-bpath','--blastpath',help='subdirectory location of blastn')
args = parser.parse_args()
sequence_data = args.inputfile # open(args.inputfile).read() 
#sequence_data = "blast_example_oneseq.fasta"
database = args.database 
# database = '/gpfs23/scratch/neuertg/Mouse_Genomic_plus_Transcript_plus_rRNA/mouse_genomic_transcript_PlusrRNA'
# database = 'D:\Mouse_Genomic_plus_Transcript_plus_rRNAv2\Mouse_Genomic_plus_Transcript_plus_rRNAv2'
outputfile = args.outputfile
pl = args.probelengthmin
mhs = args.homologysizemin
blastPath = args.blastpath
#outputfile = 'results3.xml'
import os
import time
import subprocess
start = time.time()
#Wthout quotes
#strand by default is both [
#strand="minus"  #strand="plus", strand="both"  [
#Limit output alignment percentage (perc_identity)  15/length(sequence)
per_id = 100*mhs/pl
current_dir = os.getcwd()
location = os.path.join(current_dir,blastPath)
outfmt = 5
reward = 1
penalty = -3
word_size = 7
gapopen = 5
gapextend = 2
evalue = 1000
num_alignments = 1000
strand = 'plus'
dust = 'no'
sub_command = location + ' -out ' + outputfile + ' -outfmt ' + str(outfmt) + ' -num_alignments ' + str(num_alignments) + ' -query ' + sequence_data + ' -db ' + database + ' -evalue ' + str(evalue) + ' -word_size ' + str(word_size) + ' -gapopen ' + str(gapopen) + ' -gapextend ' + str(gapextend) + ' -strand ' + strand + ' -penalty ' + str(penalty) + ' -reward ' + str(reward) + ' -dust '+ dust + ' -perc_identity ' + str(per_id)
print(sub_command)
result = subprocess.run(sub_command, shell=True,stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
combined_output = result.stdout.decode('utf-8')
print(combined_output)                     


#cline = NcbiblastnCommandline(cmd=location,db=database,query=sequence_data,out=outputfile,outfmt=5,reward=1,penalty=-3,word_size=7,gapopen=5,gapextend=2,strand='plus',num_alignments=1000,dust='no',evalue=1000,perc_identity = per_id)
#print(type(cline))
#print(cline)
#tdout, stderr = cline()
#blastn -out results.xml -outfmt 5 -num_alignments 1000 -query blast_example_oneseq.fasta -db /gpfs23/scratch/neuertg/Mouse_Genomic_plus_Transcript_plus_rRNAv2/mouse_genomic_transcript_PlusrRNA -evalue 1000 -word_size 7 -gapopen 5 -gapextend 2 -penalty -3 -reward 1 -dust no
#from Bio.Blast import NCBIXML
#E_VALUE_THRESH = 1e-2
#for record in NCBIXML.parse(open("results.xml")): 
#  if record.alignments: 
#    print("\n") 
#    print("query: %s" % record.query[:100]) 
#    for align in record.alignments: 
#      for hsp in align.hsps: 
#        if hsp.expect < E_VALUE_THRESH: 
#          print("match: %s " % align.title[:100])
#stop = time.time()
#print("The time of the run:", stop - start)

#def python_BLAST_local(fasta_file):

#return NCBIXML.parse(open("results.xml"))
