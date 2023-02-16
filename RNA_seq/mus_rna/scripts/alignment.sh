#!/bin/bash


######################################################################################################################################
##################################################     BUILD INDEX    ################################################################
######################################################################################################################################
GENBUILD=~/Desktop/rna_seq/mus_rna/output/gencode_build

# Unzip the .gz file, needed as STAR only accepts .fa
gunzip -c ${GENBUILD}/GRCm39.primary_assembly.genome.fa.gz > ${GENBUILD}/GRCm39.primary_assembly.genome.fa
gunzip -c ${GENBUILD}/gencode.vM32.primary_assembly.annotation.gtf.gz > ${GENBUILD}/gencode.vM32.primary_assembly.annotation.gtf

STAR --runThreadN 4 --runMode genomeGenerate --genomeSAsparseD 3 --limitGenomeGenerateRAM 15000000000 --genomeSAindexNbases 12 \
        --genomeDir ${GENBUILD}/starIndex --genomeFastaFiles ${GENBUILD}/GRCm39.primary_assembly.genome.fa \
        --sjdbGTFfile ${GENBUILD}/gencode.vM32.primary_assembly.annotation.gtf


# removing unzipped files to save space
rm ${GENBUILD}/GRCm39.primary_assembly.genome.fa && rm ${GENBUILD}/gencode.vM32.primary_assembly.annotation.gtf # Removing uncompressed GTF file




######################################################################################################################################
#################################################     ALIGNING READS     #############################################################
######################################################################################################################################
# Getting the accession list to pass the paired reads
SRAPATH=~/Desktop/rna_seq/mus_rna/input/SraRunTable.txt
#tail -n +2 ${SRAPATH} | grep "[^-]B cells from spleen" | grep "RNA-Seq" |cut -d ',' -f 1 > SRR_Acc_List.txt
RDSPATH=~/Desktop/rna_seq/mus_rna/output/reads/trREADS
ALIGNOUT=~/Desktop/rna_seq/mus_rna/output/starAlign

# Pass the Acession list to align the trimmed paired sample reads in parallel
cat SRR_Acc_List.txt | \
        parallel -j 1 "STAR --runThreadN 4 --genomeLoad LoadAndKeep --genomeDir ${GENBUILD}/starIndex" \
        "--readFilesIn ${RDSPATH}/{}_1P.fastq ${RDSPATH}/{}_2P.fastq --outFilterIntronMotifs RemoveNoncanonicalUnannotated" \
        "--outFileNamePrefix ${ALIGNOUT}/{} --limitBAMsortRAM 5000000000 --outSAMtype BAM SortedByCoordinate" \
        "--outReadsUnmapped Fastx"





######################################################################################################################################
#################################################     MICROBIAL CONTAMINATION     ####################################################
######################################################################################################################################

KRAKENDB=https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v2_8GB_201904.tgz #MiniKraken link
DBOUT=~/Desktop/rna_seq/mus_rna/output/contamination/micro #output for classification and unzipping

# Download miniKraken if not enough computational resources
wget ${KRAKENDB} -P ${DBOUT}

# unzip the folder
tar zxvf ${DBOUT}/minikraken2_v2_8GB_201904.tgz -C ${DBOUT}

# Run kraken 2 for taxa classification on the unmapped reads from star alignment
cat SRR_Acc_List.txt | \
        parallel -j 1 "kraken2 --db ${DBOUT}/minikraken2_v2_8GB_201904_UPDATE --threads 4 --paired ${ALIGNOUT}/{}Unmapped.out.mate1  ${ALIGNOUT}/{}Unmapped.out.mate>        "--memory-mapping --use-names -output ${DBOUT}/{}.out --report ${DBOUT}/{}.report" \

rm  ${DBOUT}/minikraken2_v2_8GB_201904_UPDATE




######################################################################################################################################
######################################################     rRNA CONTAMINATION     ####################################################
######################################################################################################################################

# Download the entire rRNA seq data from ncbi for mice using : https://www.ncbi.nlm.nih.gov/nuccore/?term=txid10090%5Borganism%3Aexp%5D
# To Download > click "send to" below the search bar > Tick "Complete Record" > In Destion field select "File" > select "FASTA" as format
rRNAPATH=~/Desktop/rna_seq/mus_rna/output/contamination/rRNA

# Create an rRNA index using bwa and the downloaded rRNA fasta file
bwa index ${rRNAPATH}/mus_rRNA.fasta

# Alignment using MEM algorithm. Outputs sam, so convert it to bam and output some stats
cat SRR_Acc_List.txt | parallel -j 2 "bwa mem -t 4 -P ${rRNAPATH}/mus_rRNA.fasta ${RDSPATH}/{}_1P.fastq ${RDSPATH}/{}_2P.fastq |" \
        "samtools view -@ 4 -b -o ${rRNAPATH}/bwaAlign/{}_rRNA.bam && samtools flagstat -@ 6 ${rRNAPATH}/bwaAlign/{}_rRNA.bam > ${rRNAPATH}/bwaAlign/{}.out"
