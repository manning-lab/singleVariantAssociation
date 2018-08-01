task summaryCSV {
	String pval
	Float? pval_threshold
	String label
	File assoc

	Int memory
	Int disk

	command {
		R --vanilla --args ${pval} ${default="0.0001" pval_threshold} ${label} ${assoc} < /singleVariantAssociation/summaryCSV.R
	}
	
	runtime {
		docker: "manninglab/singlevariantassociation:latest"
  	    disks: "local-disk ${disk} SSD"
        memory: "${memory} GB"
        bootDiskSizeGb: 20
	}

	output {
		File topassoccsv = "${label}.topassoc.csv"
		File plots = "${label}_association_plots.png"
	}
}

workflow w_summaryCSV {
	String this_pval
	Float? this_pval_threshold
	String this_label
	File this_assoc

	Int this_memory
	Int this_disk

	call summaryCSV {
		input: pval = this_pval, pval_threshold = this_pval_threshold, label = this_label, assoc = this_assoc, memory = this_memory, disk = this_disk
	}
}