#!/bin/bash

GENBUILD=~/Desktop/nf_tut/rna_seq/mus_rna/output/gencode_build

## Unzip the .gz file, needed as STAR only accepts .fa
gunzip -c ${GENBUILD}/GRCm39.primary_assembly.genome.fa.gz > ${GENBUILD}/GRCm39.primary_assembly.genome.fa
gunzip -c ${GENBUILD}/gencode.vM32.primary_assembly.annotation.gtf.gz > ${GENBUILD}/gencode.vM32.primary_assembly.annotation.gtf


## --genomeSAsparseD, --limitGenomeGenerateRAM, --genomeSAindexNb flags for reducing RAM requirements but runs slower
STAR --runThreadN 4 --runMode genomeGenerate --genomeSAsparseD 3 --limitGenomeGenerateRAM 16000000000 --genomeSAindexNbases 12 \
	--genomeDir ${GENBUILD}/starIndex --genomeFastaFiles ${GENBUILD}/GRCm39.primary_assembly.genome.fa \
	--sjdbGTFfile ${GENBUILD}/gencode.vM32.primary_assembly.annotation.gtf


rm ${GENBUILD}/GRCm39.primary_assembly.genome.fa && rm ${GENBUILD}/gencode.vM32.primary_assembly.annotation.gtf
