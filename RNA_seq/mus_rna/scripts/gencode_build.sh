#! /bin/bash


#######################################################################################
##############################   BUILD  DOWNLOAD   ####################################
#######################################################################################

# Build the mouse genome from gencode build
REFLNK=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/GRCm39.primary_assembly.genome.fa.gz
GTFLNK=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.primary_assembly.annotation.gtf.gz
BUILDDIR=~/Desktop/nf_tut/rna_seq/mus_rna/output/gencode_build # path for output dir of assembly

# Use wget to get to download ref genome and annotation file for building index
#wget ${REFLNK} -P ${BUILDDIR}
#wget ${GTFLNK} -P ${BUILDDIR}





########################################################################################
###################################   READS  DOWNLOAD   ################################
########################################################################################

#Get the reads for p53+/- mice with radiation vs control

SRAPATH=~/Desktop/nf_tut/rna_seq/mus_rna/input/SraRunTable.txt
# GET the SRR list just for the samples needed for investigation
VAR=$(tail -n +2 ${SRAPATH} | grep "[^-]B cells from spleen" | grep "RNA-Seq" |cut -d ',' -f 1)
READSDIR=~/Desktop/nf_tut/rna_seq/mus_rna/output/reads



: '

# This is a loop for downloading the data using SRA-toolkit. Recommended to use prefetch to download
# The conditional statement checks if the file has already been downloaded in case of a crash

for i in ${VAR}
    do
        if [ -f ${i} ]
            then
                echo "${i} already downloaded"
        else
            echo "(o) Downloading SRA entry: ${i}" 
            # downloading SRA entry
            prefetch ${i}
	    fastq-dump --defline-qual '+' --split-files --gzip ${i}/${i}.sra -O ${READSDIR} && rm -r ${i}
            echo "(o) Done downloading ${i}"
        fi
    done

'

# Large dataset so this is just the first 250k reads for each file instead to accomedate for the low system resources
#cat ~/Desktop/nf_tut/rna_seq/mus_rna/input/SRA_list.txt | parallel wget {} -P ${READSDIR}




##########################################################################################################
#######################################   QC   ###########################################################
##########################################################################################################

#For looking at  pre trimmed quality check and filetering
# PRE QC check
QCOUT=~/Desktop/nf_tut/rna_seq/mus_rna/output/results/fastqc
#fastqc -t 8 ${READSDIR}/*fastq.gz -o ${QCOUT} && multiqc ${QCOUT} -o ${QCOUT}



# No adaptor contamination found but using it for demonstration purposes
#Using the provided illumina adaptor sequence to remove potential contaminations 
ADPATH=~/miniconda3/envs/nexflo/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa
printf '%s\n' "${VAR[@]}" | \
	parallel -j 2 "trimmomatic PE -threads 4 ${READSDIR}/{}_1.fastq.gz ${READSDIR}/{}_2.fastq.gz" \
	"-baseout ${READSDIR}/trREADS/{}.fastq ILLUMINACLIP:${ADPATH}:2:30:10:2:keepBothReads MINLEN:35 2>${READSDIR}/trREADS/{}trimming.log"

#POST QC check
multiqc ${READSDIR}/trREADS/ -o ${QCOUT}/trREADS/

