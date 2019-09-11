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

	# groupdata
	fra_glo_groupdata <- read.table(path_groupdata, col.names = c("Group", "Name"), header = F, as.is = T)
	tryCatch(
		{
			for(i in 1:nrow(fra_glo_groupdata)){
				group_name <- fra_glo_groupdata[i,1]
				group_cont <- fra_glo_groupdata[i,2]
				group_cont <- strsplit(group_cont, split = ",")[[1]]
				for(accession in group_cont){
					ID <- as.character(fra_glo_metadata[which(fra_glo_metadata$AccessionName == accession),]$vcfID)
					if(!ID %in%  bcf_header){
						cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Detect file ") %+% black$bgWhite$bold(path_groupdata) %+% black$bgWhite$bold(" is not correct.") %+% "\n"
							%+% black$bgWhite$bold("    ") %+% black$bgWhite$bold("AccessionName ") %+% red$bgWhite$bold(accession) %+% black$bgWhite$bold(" does not appear in the vcf files.") %+% "\n"
							)
					}
				}
			}
		}, error = function(e) {
				cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Something wrong while rechecking ") %+% black$bgWhite$bold(path_groupdata) %+% black$bgWhite$bold(".") %+% "\n"
					%+% black$bgWhite$bold("    ") %+% black$bgWhite$bold(e) %+% "\n"
					)
		}
	)
	cat(green$bgWhite$bold("[+] Second round file checking finished") %+% "\n")
}