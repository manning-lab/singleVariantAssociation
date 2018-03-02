task conditionalPhenotype {
	Array[File] genotype_files
	File? phenotype_file
	String? id_col
	File? sample_file
	String? snps
	String label

	Int disk

	command {
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
		R --vanilla --args ${genotype_file} ${phenotype_file} ${outcome_name} ${outcome_type} ${default="NA" covariates_string} ${default="NA" conditional_string} ${default="NA" ivars_string} ${default="NA" group_var} ${default="NA" sample_file} ${label} ${kinship_matrix} ${id_col} < /singleVariantAssociation/genesis_nullmodel.R
	}

	runtime {
		docker: "tmajarian/singlevariantassociation:latest"
		disks: "local-disk ${disk} SSD"
		memory: "${memory}G"
	}

	output {
		File model = "${label}_null.RDa"
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
