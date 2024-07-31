process PROCESS_INPUT {
	debug true

	input:

	output:
	path "input.csv"
	
	shell:
	base_dir = params.base_dir
	gff_filename = params.gff
	fasta_filename = params.fasta
	
	'''
	> input.csv

	for species_dir in !{base_dir}/*
	do
		species=$(basename ${species_dir})

		for assembly_dir in ${species_dir}/*
		do
			assembly=$(basename ${assembly_dir})

			for genebuild_provider_dir in ${assembly_dir}/*
			do
				genebuild_provider=$(basename ${genebuild_provider_dir})

				# get gff file
				gff_file=None
				if [[ -d ${genebuild_provider_dir}/geneset ]]
				then
					last_update_date=$(ls ${genebuild_provider_dir}/geneset | sort -r | head -n 1)
					file=${genebuild_provider_dir}/geneset/${last_update_date}/!{gff_filename}

					if [[ -e ${file} ]]
					then
						gff_file=${file}
					fi
				fi

				# get fasta file
				fasta_file=None
				if [[ -d ${genebuild_provider_dir}/genome ]]
				then
					file=${genebuild_provider_dir}/genome/!{fasta_filename}

					if [[ -e ${file} ]]
					then
							fasta_file=${file}
					fi
				fi
	
				# print csv line
				if [[ ${gff_file} != None && ${fasta_file} != None ]]
				then
					outdir_prefix=${species}/${assembly}/${genebuild_provider}
					echo "${gff_file},${fasta_file},${outdir_prefix}" >> input.csv
				fi
			done
		done
	done

	'''
}
