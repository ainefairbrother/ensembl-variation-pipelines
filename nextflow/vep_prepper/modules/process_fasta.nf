process PROCESS_FASTA {
	input:
	tuple path(fasta), val(outdir_suffix)

	output:

	shell:
	out_dir = params.outdir + "/" + outdir_suffix
	out_filename = fasta.getName().replace(".gz", ".bgz")

	'''
	gzip -dc !{fasta} | bgzip -c > !{out_filename}
	samtools faidx !{out_filename}

	# move
	mkdir -p !{out_dir}
    mv *.bgz* !{out_dir}/
	'''
}
