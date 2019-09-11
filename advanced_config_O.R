# filepaths
path_fa_index <- paste(path_datafolder, "/", name_fasta, ".fai", sep="")
path_geneindex <- paste(path_datafolder, "/", "geneinfo.txt", sep="")
path_vcf <- paste(path_datafolder, "/", "output.ann.bcf.gz", sep="")

path_gff3 <- paste(path_datafolder, "/", name_gff3, ".gz", sep="")
path_raw_gff3 <- paste(path_datafolder, "/", name_gff3, sep="")
path_fasta <- paste(path_datafolder, "/", name_fasta, sep="")

path_metadata <- paste(path_datafolder, "/", name_metadata, sep="")
path_groupdata <- paste(path_datafolder, "/", name_groupdata, sep="")
path_sam_location <- paste(path_datafolder, "/", name_sam_location, sep="")

if(is.na(name_sysinfo)){
	path_sysinfo <- NA
}else{
	path_sysinfo <- paste(path_datafolder, "/", name_sysinfo, sep="")
}


if(is.na(name_UIsetting)){
	path_UIsetting <- NA
}else{
	path_UIsetting <- paste(path_datafolder, "/", name_UIsetting, sep="")
}
