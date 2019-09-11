file_check <- function(status){
	cat(blue$bgWhite$bold("[ ] Begining the first round file checking") %+% "\n")

	file_sub_list <- c(path_fasta, path_raw_gff3, path_metadata, path_groupdata, path_sam_location)
	for(i in file_sub_list){
		if(!file.exists(i)){
			cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Detect file ") %+% red$bgWhite$bold(i) %+% black$bgWhite$bold(" is not exist or its path is wrong.") %+% "\n")
			status <- 1
		}
	}

	# metadata
	is_correct <- tryCatch(
		{
			fra_glo_metadata <- read.table(path_metadata, header = T, as.is = T, fill=T, sep="\t")
			TRUE
		}, error = function(e) {
			cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Detect file ") %+% black$bgWhite$bold(path_metadata) %+% black$bgWhite$bold(" is not correct.") %+% "\n"
				%+% black$bgWhite$bold("    ") %+% black$bgWhite$bold(e) %+% "\n"
				)
			status <- 1
			FALSE
		}
	)
	if(is_correct){
		if(length(unique(fra_glo_metadata$DisplayName)) != length(fra_glo_metadata$DisplayName)){
			cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Detect file ") %+% black$bgWhite$bold(path_metadata) %+% black$bgWhite$bold(" is not correct.") %+% "\n"
				%+% black$bgWhite$bold("    ") %+% black$bgWhite$bold("DisplayName column (third column) should be unique.") %+% "\n"
				)
			status <- 1
		}
		if(length(unique(fra_glo_metadata$AccessionName)) != length(fra_glo_metadata$AccessionName)){
			cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Detect file ") %+% black$bgWhite$bold(path_metadata) %+% black$bgWhite$bold(" is not correct.") %+% "\n"
				%+% black$bgWhite$bold("    ") %+% black$bgWhite$bold("AccessionName column (second column) should be unique.") %+% "\n"
				)
			status <- 1
		}
		if(length(unique(fra_glo_metadata$vcfID)) != length(fra_glo_metadata$vcfID)){
			cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Detect file ") %+% black$bgWhite$bold(path_metadata) %+% black$bgWhite$bold(" is not correct.") %+% "\n"
				%+% black$bgWhite$bold("    ") %+% black$bgWhite$bold("vcfID column (first column) should be unique.") %+% "\n"
				)
			status <- 1
		}
	}
	# groupdata
	is_correct <- tryCatch(
		{
			fra_glo_groupdata <- read.table(path_groupdata, col.names = c("Group", "Name"), header = F, as.is = T)
			TRUE
		}, error = function(e) {
			cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Detect file ") %+% black$bgWhite$bold(path_groupdata) %+% black$bgWhite$bold(" is not correct.") %+% "\n"
				%+% black$bgWhite$bold("    ") %+% black$bgWhite$bold(e) %+% "\n"
				)
			status <- 1
			FALSE
		}
	)
	if(is_correct){
		if(length(unique(fra_glo_groupdata$Group)) != length(fra_glo_groupdata$Group)){
			cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Detect file ") %+% black$bgWhite$bold(path_groupdata) %+% black$bgWhite$bold(" is not correct.") %+% "\n"
				%+% black$bgWhite$bold("    ") %+% black$bgWhite$bold("Group names in first column should be unique.") %+% "\n"
				)
			status <- 1
		}
		for(i in 1:nrow(fra_glo_groupdata)){
			s <- fra_glo_groupdata[i, 2]
			s <- strsplit(s, split = ",")[[1]]
			for(j in s){
				if(!(j %in% fra_glo_metadata$AccessionName)){
					cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Detect file ") %+% black$bgWhite$bold(path_groupdata) %+% black$bgWhite$bold(" is not correct.") %+% "\n"
						%+% black$bgWhite$bold("    ") %+% black$bgWhite$bold("Wrong sample label in group ") %+% red$bgWhite$bold(fra_glo_groupdata[i, 1]) %+% "\n" 
						)
					status <- 1
				}
				break
			}
		}
	}

	###############################################################################################
	# Finish
	cat(green$bgWhite$bold("[+] File checking finished") %+% "\n")

	#return
	status
}