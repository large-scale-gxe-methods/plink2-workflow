task process_phenos {
	
	File phenofile
	String sample_id_header
	String outcome
	String exposure
	String? covar_names
	String? delimiter
	String? missing

	command {
		python3 /format_plink2_phenos.py ${phenofile} ${sample_id_header} ${outcome} ${exposure} "${covar_names}" "${delimiter}" ${missing}
	}

	runtime {
		docker: "quay.io/large-scale-gxe-methods/plink2-workflow"
		memory: "2 GB"
	}

        output {
                File pheno_fmt = "plink2_phenotypes.txt"
                File plink2_parameter_file = "plink2_parameters_string.txt"
	}
}


task run_interaction {
  
	File genofile_pgen
	File genofile_psam
	File genofile_pvar
	File phenofile
	String outcome
	Boolean binary_outcome
	String exposure
	String? covar_names
	File plink2_parameter_file
	Int? memory
	Int? disk
	Int threads
	Int monitoring_freq
	
	String covar_name_str = exposure + " " + covar_names
	String plink2_parameter_string = read_string(plink2_parameter_file)

        command {
		dstat -c -d -m --nocolor 1 > system_resource_usage.log &
		atop -x -P PRM 1 | grep '(GEM)' > process_resource_usage.log &

		/plink2 --pgen ${genofile_pgen} \
			--psam ${genofile_psam} \
			--pvar ${genofile_pvar} \
			--allow-extra-chr \
			--pheno-name ${outcome} \
			${true="--1" false="" binary_outcome} \
			--pheno ${phenofile} \
			--covar-name ${covar_name_str} \
			--glm interaction \
			--parameters ${plink2_parameter_string} \
			--threads ${threads} \
			--out plink2_res

			mv plink2_res.${outcome}.glm.${true='logistic' false='linear' binary_outcome} plink2_res
   		 }


	runtime {
		docker: "quay.io/large-scale-gxe-methods/plink2-workflow"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
	}

        output {
		File res = "plink2_res"
		File system_resource_usage = "system_resource_usage.log"
		File process_resource_usage = "process_resource_usage.log"
        }
}

task standardize_output {

	File resfile
	String exposure
	Boolean binary_outcome
	String outfile_base = basename(resfile)
	String outfile = "${outfile_base}.fmt"

	command {
		python3 /format_plink2_output.py ${resfile} ${exposure} ${binary_outcome} ${outfile}
	}

	runtime {
		docker: "quay.io/large-scale-gxe-methods/plink2-workflow"
		memory: "2 GB"
	}

        output {
                File res_fmt = "${outfile}"
	}
}

task cat_results {

	Array[File] results_array

	command {
		head -1 ${results_array[0]} > all_results.txt && \
			for res in ${sep=" " results_array}; do tail -n +2 $res >> all_results.txt; done
	}
	
	runtime {
		docker: "quay.io/large-scale-gxe-methods/plink2-workflow"
		disks: "local-disk 5 HDD"
	}
	output {
		File all_results = "all_results.txt"
	}
}
			

workflow run_plink2 {

	Array[File] genofiles_pgen
	Array[File] genofiles_psam
	Array[File] genofiles_pvar
	File phenofile
	String sample_id_header
	String outcome
	Boolean binary_outcome
	String exposure_names
	String? covar_names = ""
	String? delimiter = ","
	String? missing = "NA"
	Boolean? robust
	Int? memory = 10
	Int? disk = 20
	Int? threads = 1
	Int? monitoring_freq = 1

	call process_phenos {
		input:
			phenofile = phenofile,
			sample_id_header = sample_id_header,
			outcome = outcome,
			exposure = exposure_names,
			covar_names = covar_names,
			delimiter = delimiter,
			missing = missing
	}

	scatter (i in range(length(genofiles_pgen))) {
 		call run_interaction {
 			input:
 				genofile_pgen = genofiles_pgen[i],
 				genofile_psam = genofiles_psam[i],
 				genofile_pvar = genofiles_pvar[i],
 				phenofile = process_phenos.pheno_fmt,
 				outcome = outcome,
 				binary_outcome = binary_outcome,
				exposure = exposure_names,
 				covar_names = covar_names,
 				plink2_parameter_file = process_phenos.plink2_parameter_file,
 				memory = memory,	
 				disk = disk,
				threads = threads,
				monitoring_freq = monitoring_freq
 		}
 	}
 
 	scatter (resfile in run_interaction.res) {
 		call standardize_output {
 			input:
 				resfile = resfile,
 				exposure = exposure_names,
				binary_outcome = binary_outcome
 		}
 	}	
 
 	call cat_results {
 		input:
 			results_array = standardize_output.res_fmt
 	}

         output {
		File results = cat_results.all_results
 		Array[File] system_resource_usage = run_interaction.system_resource_usage
 		Array[File] process_resource_usage = run_interaction.process_resource_usage
 	}

	parameter_meta {
		genofiles_pgen: "Array of PLINK2 genotype (.pgen) filepaths."
		genofiles_psam: "Array of PLINK2 sample (.psam) filepaths."
		genofiles_pvar: "Array of PLINK2 variant (.pvar) filepaths."
		phenofile: "Phenotype filepath. Does not need to be in PLINK format (will be processed as part of the workflow)."	
		sample_id_header: "Optional column header name of sample ID in phenotype file."
		outcome: "Column header name of phenotype data in phenotype file."
                binary_outcome: "Boolean: is the outcome binary? Otherwise, quantitative is assumed."
		exposure_names: "Column header name(s) of the exposures for genotype interaction testing (space-delimited). Only one exposures is currently allowed."
		covar_names: "Column header name(s) of any covariates for which only main effects should be included (space-delimited). This set should not overlap with exposures or int_covar_names."
		delimiter: "Delimiter used in the phenotype file."
		missing: "Missing value key of phenotype file."
		cpu: "Minimum number of requested cores."
		disk: "Requested disk space (in GB)."
		monitoring_freq: "Delay between each output for process monitoring (in seconds). Default is 1 second."
	}

	meta {
		author: "Kenny Westerman"
		email: "kewesterman@mgh.harvard.edu"
		description: "Run interaction tests using PLINK2 and return summary statistics for 1-DF and 2-DF tests."
	}
}
