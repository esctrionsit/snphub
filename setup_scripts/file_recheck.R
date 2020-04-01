file_recheck <- function(){
	cat(blue$bgWhite$bold("[ ] Begining the second round file checking") %+% "\n")
	bcf_header <- readLines(pipe(paste(path_bcftools, "view", path_vcf, "-h", sep = " ")))
	bcf_header <- bcf_header[length(bcf_header)]
	bcf_header <- strsplit(bcf_header, split = "\t")[[1]]
	# metadata
	fra_glo_metadata <- read.table(path_metadata, header = T, as.is = T, fill=T, sep="\t")
	for(ID in fra_glo_metadata$vcfID){
		if(!ID %in%  bcf_header){
			cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Detect file ") %+% black$bgWhite$bold(path_metadata) %+% black$bgWhite$bold(" is not correct.") %+% "\n"
				%+% black$bgWhite$bold("    ") %+% black$bgWhite$bold("vcfID ") %+% red$bgWhite$bold(ID) %+% black$bgWhite$bold(" does not appear in the vcf files.") %+% "\n"
				)
		}
	}

	cat(green$bgWhite$bold("[+] Second round file checking finished") %+% "\n")
}