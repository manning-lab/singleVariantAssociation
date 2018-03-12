task conditionalPhenotype {
	Array[File] genotype_files
	File? phenotype_file
	String? id_col
	File? sample_file
	String? snps
	String label

	Int disk

	command {
		echo "Input files" > conditionalPhenotype.log
		echo "genotype_files: ${sep="," genotype_files}" >> conditionalPhenotype.log
		echo "phenotype_file: ${phenotype_file}" >> conditionalPhenotype.log
		echo "id_col: ${id_col}" >> conditionalPhenotype.log
		echo "sample_file: ${sample_file}" >> conditionalPhenotype.log
		echo "snps: ${snps}" >> conditionalPhenotype.log
		echo "label: ${label}" >> conditionalPhenotype.log
		echo "disk: ${disk}" >> conditionalPhenotype.log
		echo "" >> conditionalPhenotype.log
		dstat -c -d -m --nocolor 10 1>>conditionalPhenotype.log &
		R --vanilla --args ${sep="," genotype_files} ${phenotype_file} ${id_col} ${default="NA" sample_file} ${snps} ${label} < /singleVariantAssociation/preprocess_conditional.R
	}

	runtime {
		docker: "tmajarian/singlevariantassociation:latest"
		disks: "local-disk ${disk} SSD"
		memory: "10G"
	}

	output {
		File new_phenotype_file = "${label}_phenotypes.csv"
		File alt_ref = "${label}_alleles.txt"
		File log_file = "conditionalPhenotype.log"
	}
}

task fitNull {
	File genotype_file
	File? phenotype_file
	String? outcome_name
	String? outcome_type
	String? covariates_string
	String? conditional_string
	String? ivars_string
	String? group_var
	File? sample_file
	String label
	File? kinship_matrix
	String? id_col

	Int memory
	Int disk

	command {
		echo "Input files" > fitNull.log
		echo "genotype_file: ${genotype_file}" >> fitNull.log
		echo "phenotype_file: ${phenotype_file}" >> fitNull.log
		echo "outcome_name: ${outcome_name}" >> fitNull.log
		echo "outcome_type: ${outcome_type}" >> fitNull.log
		echo "covariates_string: ${covariates_string}" >> fitNull.log
		echo "conditional_string: ${conditional_string}" >> fitNull.log
		echo "ivars_string: ${ivars_string}" >> fitNull.log
		echo "group_var: ${group_var}" >> fitNull.log
		echo "sample_file: ${sample_file}" >> fitNull.log
		echo "label: ${label}" >> fitNull.log
		echo "kinship_matrix: ${kinship_matrix}" >> fitNull.log
		echo "id_col: ${id_col}" >> fitNull.log
		echo "memory: ${memory}" >> fitNull.log
		echo "disk: ${disk}" >> fitNull.log
		echo "" >> fitNull.log
		dstat -c -d -m --nocolor 10 1>>fitNull.log &
		R --vanilla --args ${genotype_file} ${phenotype_file} ${outcome_name} ${outcome_type} ${default="NA" covariates_string} ${default="NA" conditional_string} ${default="NA" ivars_string} ${default="NA" group_var} ${default="NA" sample_file} ${label} ${kinship_matrix} ${id_col} < /singleVariantAssociation/genesis_nullmodel.R
	}

	runtime {
		docker: "tmajarian/singlevariantassociation:latest"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		File model = "${label}_null.RDa"
		File log_file = "fitNull.log"
	}
}

task assocTest {

	File gds_file
	File? null_file
	String label
	String? test
	String? ivars_string
	Int? mac
	String? variant_range

	Int memory
	Int disk

	command {
		echo "Input files" > assocTest.log
		echo "gds_file: ${gds_file}" >> assocTest.log
		echo "null_file: ${null_file}" >> assocTest.log
		echo "label: ${label}" >> assocTest.log
		echo "test: ${test}" >> assocTest.log
		echo "ivars_string: ${ivars_string}" >> assocTest.log
		echo "mac: ${mac}" >> assocTest.log
		echo "variant_range: ${variant_range}" >> assocTest.log
		echo "memory: ${memory}" >> assocTest.log
		echo "disk: ${disk}" >> assocTest.log
		echo "" >> assocTest.log
		dstat -c -d -m --nocolor 10 1>>assocTest.log &
		R --vanilla --args ${gds_file} ${null_file} ${label} ${default="Score" test} ${default="5" mac} ${default="NA" ivars_string} ${default="NA" variant_range} < /singleVariantAssociation/association.R
	}

	meta {
		author: "Tim Majarian"
		email: "tmajaria@braodinstitute.org"	
	}
	
	runtime {
		docker: "tmajarian/singlevariantassociation:latest"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		File assoc = "${label}.assoc.RData"
		File log_file = "assocTest.log"
	}
}

task summary {
	String pval
	Float? pval_threshold
	String label
	Array[File] assoc

	Int memory
	Int disk

	command {
		echo "Input files" > summary.log
		echo "pval: ${pval}" >> summary.log
		echo "label: ${label}" >> summary.log
		echo "assoc: ${sep = ',' assoc}" >> summary.log
		echo "memory: ${memory}" >> summary.log
		echo "disk: ${disk}" >> summary.log
		echo "" >> summary.log
		dstat -c -d -m --nocolor 10 1>>summary.log &
		R --vanilla --args ${pval} ${default="0.0001" pval_threshold} ${label} ${sep = ',' assoc} < /singleVariantAssociation/summary.R
	}
	
	runtime {
		docker: "tmajarian/singlevariantassociation:latest"
  	    disks: "local-disk ${disk} SSD"
        memory: "${memory}G"
	}

	output {
		File allassoccsv = "${label}.assoc.csv"
		File topassoccsv = "${label}.topassoc.csv"
		File plots = "${label}_association_plots.png"
		File log_file = "summary.log"
	}
}

workflow w_assocTest {
	# conditionalPhenotype inputs
	String? these_snps


	# fitNull inputs
	Array[File] these_genotype_files
	File? this_phenotype_file
	String? this_outcome_name
	String? this_outcome_type
	String? this_covariates_string
	String? this_ivars_string
	String? this_group_var
	File? this_sample_file
	String this_label
	File? this_kinship_matrix
	String? this_id_col
	
	# assocTest inputs
	File? this_null_file
	String? this_test
	Int? this_mac
	String? this_variant_range

	# summary inputs
	String this_pval
	Float? this_pval_threshold	

	# inputs to all
	Int this_memory
	Int this_disk

	File null_genotype_file = these_genotype_files[0]

	Boolean need_null = defined(this_null_file)

	if(defined(these_snps)) {

		call conditionalPhenotype {
			input: genotype_files = these_genotype_files, phenotype_file = this_phenotype_file, id_col = this_id_col, sample_file = this_sample_file, snps = these_snps, label = this_label, disk = this_disk
		}
		
		call fitNull as fitNullConditional {
			input: genotype_file = null_genotype_file, phenotype_file = conditionalPhenotype.new_phenotype_file, outcome_name = this_outcome_name, outcome_type = this_outcome_type, covariates_string = this_covariates_string, conditional_string = these_snps, ivars_string = this_ivars_string, group_var = this_group_var, sample_file = this_sample_file, label = this_label, kinship_matrix = this_kinship_matrix, id_col = this_id_col, memory = this_memory, disk = this_disk
		}

		scatter(this_genotype_file in these_genotype_files) {
		
			call assocTest as assocTestConditional {
				input: gds_file = this_genotype_file, null_file = fitNullConditional.model, label = this_label, test = this_test, mac = this_mac, ivars_string = this_ivars_string, variant_range = this_variant_range, memory = this_memory, disk = this_disk
			}
		}

		call summary as summaryConditional {
			input: pval = this_pval, pval_threshold = this_pval_threshold, label = this_label, assoc = assocTestConditional.assoc, memory = this_memory, disk = this_disk
		}
	}

	if (!defined(these_snps)) {
	
		if(!need_null) {
			
			call fitNull {
				input: genotype_file = null_genotype_file, phenotype_file = this_phenotype_file, outcome_name = this_outcome_name, outcome_type = this_outcome_type, covariates_string = this_covariates_string, conditional_string = these_snps, ivars_string = this_ivars_string, group_var = this_group_var, sample_file = this_sample_file, label = this_label, kinship_matrix = this_kinship_matrix, id_col = this_id_col, memory = this_memory, disk = this_disk
			}

			scatter(this_genotype_file in these_genotype_files) {
			
				call assocTest {
					input: gds_file = this_genotype_file, null_file = fitNull.model, label = this_label, test = this_test, mac = this_mac, ivars_string = this_ivars_string, variant_range = this_variant_range, memory = this_memory, disk = this_disk
				}
			}

			call summary {
				input: pval = this_pval, pval_threshold = this_pval_threshold, label = this_label, assoc = assocTest.assoc, memory = this_memory, disk = this_disk
			}

		} 

		if(need_null) {

			scatter(this_genotype_file in these_genotype_files) {
			
				call assocTest as assocNull {
					input: gds_file = this_genotype_file, null_file = this_null_file, label = this_label, test = this_test, mac = this_mac, ivars_string = this_ivars_string, variant_range = this_variant_range, memory = this_memory, disk = this_disk
				}
			}

			call summary as summaryNull {
				input: pval = this_pval, pval_threshold = this_pval_threshold, label = this_label, assoc = assocNull.assoc, memory = this_memory, disk = this_disk
			}
		}
	}
}
