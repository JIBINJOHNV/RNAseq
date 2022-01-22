import pandas as pd
import glob
import argparse
import os


print("Pleasse activate conda activate RNAseq conda environment; otherwise may fail")

parser=argparse.ArgumentParser(description="It is for Incorporating additional information and Initial filtering of Annovar out put file")
parser.add_argument('-NumberOfThreads','--NumberOfThreads', help="Number of threads", required=True)
parser.add_argument('-genomeDir','--genomeDir', help="Directory where fasta file saved; should end with /", required=True)
parser.add_argument('-gtf','--gtf', help="GTF File", required=True)
parser.add_argument('-fasta','--fasta', help="Genome fasta file name", required=True)
parser.add_argument('-ReadLength','--ReadLength', help="Read lengthe", required=True)

args=parser.parse_args()


NumberOfThreads=args.NumberOfThreads
genomeDir=args.genomeDir
gtf=args.gtf
gfasta=args.fasta
ReadLength_1=int(args.ReadLength)-1
print(ReadLength_1)


cmd= "STAR --runThreadN "+ NumberOfThreads +" --runMode genomeGenerate " +  " --genomeDir "+ genomeDir + ''' \
        --genomeFastaFiles '''+ gfasta + " --sjdbGTFfile "+ gtf + " --sjdbOverhang "+ str(ReadLength_1)

os.system(cmd)
