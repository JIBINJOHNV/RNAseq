# RNAseq
#Create STAR INDEX

python Create_STAR_INdex.py -NumberOfThreads 5 -genomeDir DenomeDirPath/ -gtf Rattus_norvegicus.mRatBN7.2.105.gtf -fasta Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.fa -ReadLength 100



Command to DE Analysis

Rscript DESEQ2_V2.R --ReadCount TestData/counts_matrix_Test.txt --SampleInfo TestData/SampleInfoe.txt --contrast TestData/Contraxt.txt --OutpuName Test


