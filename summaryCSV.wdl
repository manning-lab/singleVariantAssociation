task summaryCSV {
	String? pval
	Float? pval_threshold
	String label
	Array[File] assoc

	Int memory
	Int disk

	command {
		R --vanilla --args ${default="NA" pval} ${default="0.0001" pval_threshold} ${label} ${sep="," assoc} < /singleVariantAssociation/summaryCSV.R
	}
	
	runtime {
		docker: "manninglab/singlevariantassociation:genesis2"
  	    disks: "local-disk ${disk} SSD"
        memory: "${memory} GB"
        bootDiskSizeGb: 20
	}

	output {
		File assoccsv = "${label}.all.assoc.csv"
		File topassoccsv = "${label}.top.assoc.csv"
		File plots = "${label}.association.plots.png"
	}
}

workflow w_summaryCSV {
	String? this_pval
	Float? this_pval_threshold
	String this_label
	Array[File] these_assoc

	Int this_memory
	Int this_disk

	call summaryCSV {
		input: pval = this_pval, pval_threshold = this_pval_threshold, label = this_label, assoc = these_assoc, memory = this_memory, disk = this_disk
	}
}