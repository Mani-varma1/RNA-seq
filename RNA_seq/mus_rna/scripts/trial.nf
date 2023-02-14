//process_tuple_input.nf
nextflow.enable.dsl=2


params.loc = "~/Desktop/nf_tut/rna_seq/mus_rna/output/reads/*{1,2}.fastq.gz"
process FASTQC{
	input:
	tuple val(sample_id), path(reads)

	output:
	file "*_fastqc.{zip,html}" into fastqc_results
  
	script:
	"""
	fastqc ${reads}
	"""
}

/*process MULTIQC{
	input:
	file ('fastqc/*') from fastqc_results.collect().ifEmpty([])
	
	output:
	file "multiqc_report.html" into multiqc_report
	file "multiqc_data"
	
	script:
	"""
	multiqc .
	"""
}
*/
reads_ch = Channel.fromFilePairs(params.loc)
println reads_ch
workflow {
	FASTQC(reads_ch)
	//MULTIQC(fastqc_results)
}
