task process_phenos {
	
	File phenofile
	String sample_id_header
	String outcome
	String covar_headers
	String exposure
	String? delimiter = ","
	String? missing = "NA"

	command {
		python3 /format_plink2_phenos.py ${phenofile} ${sample_id_header} ${outcome} "${covar_headers}" ${exposure} "${delimiter}" ${missing}
	}

	runtime {
		docker: "quay.io/large-scale-gxe-methods/plink2-workflow"
		memory: "2 GB"
	}

        output {
                File pheno_fmt = "plink2_phenotypes.csv"
	}
}


task run_interaction {
  
    File genofile_pgen
    File genofile_psam
    File genofile_pvar
    File phenofile
    String outcome
	Boolean binary_outcome
	String covar_headers
	Int? memory = 10
	Int? disk = 20

        command {
        	echo "" > resource_usage.log
        	dstat -c -d -m --nocolor 10 1>>resource_usage.log &
			/plink2 --pgen ${genofile_pgen} \
				--psam ${genofile_psam} \
				--pvar ${genofile_pvar} \
				--allow-extra-chr \
				--pheno-name ${outcome} \
				--pheno ${phenofile} \
				--covar-name ${covar_headers} \
				--glm interaction \
				--parameters 1-7 \
				--out plink2_res
   		 }


	runtime {
		docker: "quay.io/large-scale-gxe-methods/plink2-workflow"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
	}

        output {
            File res = "plink2_res.${outcome}.glm.linear"
			File resource_usage = "resource_usage.log"
        }
}

task standardize_output {

	File resfile
	String exposure
	String outfile_base = basename(resfile)
	String outfile = "${outfile_base}.fmt"

	command {
		python3 /format_plink2_output.py ${resfile} ${exposure} ${outfile}
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
	String covar_headers
	String exposures
	String? delimiter
	String? missing
	Boolean? robust
	Int? memory
	Int? disk

	call process_phenos {
		input:
			phenofile = phenofile,
			sample_id_header = sample_id_header,
			outcome = outcome,
			covar_headers = covar_headers,
			exposure = exposures,
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
 				covar_headers = covar_headers,
 				memory = memory,	
 				disk = disk
 		}
 	}
 
 	scatter (resfile in run_interaction.res) {
 		call standardize_output {
 			input:
 				resfile = resfile,
 				exposure = exposures
 		}
 	}	
 
 	call cat_results {
 		input:
 			results_array = standardize_output.res_fmt
 	}

         output {
                 File results = cat_results.all_results
 		Array[File] resource_usage = run_interaction.resource_usage
 	}

	parameter_meta {
 	genofiles: "Array of genotype filepaths in .bgen format. NEEDS UPDATING RE: PLINK FORMAT"
 	phenofile: "Phenotype filepath."
 	sample_id_header: "Column header name of sample ID in phenotype file."
 	outcome: "Column header name of phenotype data in phenotype file."
 	binary_outcome: "Boolean: is the outcome binary? Otherwise, quantitative is assumed."
 	covar_headers: "Column header names of the selected covariates in the pheno data file (space-delimited)."
 	exposures: "Column header name(s) of the covariates to use as exposures for genotype interaction testing (space-delimited). All exposures must also be provided as covariates."
 	delimiter: "Delimiter used in the phenotype file."
 	missing: "Missing value key of phenotype file."
 	memory: "Requested memory (in GB)."
 	disk: "Requested disk space (in GB)."
 	}

	meta {
		author: "Kenny Westerman"
		email: "kewesterman@mgh.harvard.edu"
		description: "Run interaction tests using the ProbABEL package and return summary statistics for 1-DF and 2-DF tests."
	}
}
