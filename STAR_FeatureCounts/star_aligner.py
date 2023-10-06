import os

sample_names=["sample1","sample2","sample3","sample4","sample5" ]

fasta_file="Pathto_fasta_file"
star_index="path_to_star_index"

for sample_name in sample_names:
    os.system(f'''STAR
        --readFilesIn {sample_name}.fastp.R1.gz {sample_name}.fastp.R1.gz \
        --outSAMattrRGline "ID:your_rg_id SM:{sample_name} PL:illumina LB:Lib1" \
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
        --outFileNamePrefix {sample_name} \
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
        --chimMainSegmentMultNmax 1 \
        --outSAMattributes NH HI AS nM NM ch''')

