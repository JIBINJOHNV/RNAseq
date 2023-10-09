vi GSE135251_Part_200_216_Star_alignment.py



nohup python3 GSE135251_Part_200_216_Star_alignment.py 2>&1 > GSE135251_Part_200_216_Star_alignment.out &


import os
import glob

#fastq_path="/data/amrendra/Analysis/Cleanedfastqfiles/GSE126848/Fastp/" #47 samplesa single end, 75bp
#star_index="/data/amrendra/reference/TCGA/star-2.7.5c_GRCh38.d1.vd1_gencode.v36_75bp/"

fastq_path="/data/amrendra/Analysis/Cleanedfastqfiles/GSE135251/Fastp/" #
star_index="/data/amrendra/reference/TCGA/star-2.7.5c_GRCh38.d1.vd1_gencode.v36/"


sample_list=glob.glob(f'{fastq_path}*')

sample_names=[x.split("/")[-1].split("_")[0] for x in sample_list ]
sample_names=list(set(sample_names))
sample_names=sorted(sample_names)

n=1
for sample_name in sample_names[200:]:
    print(f' Now {n}th  {sample_name} is running \n\n')
    os.system(f'''/home/amrendra/software/STAR-2.7.11a/bin/Linux_x86_64/STAR \
        --readFilesIn {fastq_path}{sample_name}_fastp_R1.gz  {fastq_path}{sample_name}_fastp_R2.gz  \
        --outSAMattrRGline "ID:{sample_name} SM:{sample_name} PL:illumina LB:Lib1 DS:GSE126848" \
        --genomeDir {star_index} \
        --readFilesCommand zcat \
        --runThreadN 15 \
        --twopassMode Basic \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFilterType BySJout \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --limitSjdbInsertNsj 1200000 \
        --outFileNamePrefix {sample_name}_ \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs None \
        --alignSoftClipAtReferenceEnds Yes \
        --quantMode TranscriptomeSAM GeneCounts \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --genomeLoad NoSharedMemory \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
        --chimOutJunctionFormat 1 \
        --chimMainSegmentMultNmax 1 ''')
    n=n+1
        
