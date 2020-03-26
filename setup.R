source("setup.conf")
source("advanced_config.R")
source("./setup_scripts/envir_check.R")
source("./setup_scripts/file_check.R")
source("./setup_scripts/file_make.R")
source("./setup_scripts/file_recheck.R")

status <- 0

if(data_type != "VCF" && data_type != "vcf" && data_type != "hapmap"){
	status = 1
	print("Wrong data_type in \"setup.conf\"!")
}

envir_check(status)

if(status != 1 && data_type == "hapmap"){
	if(!file.exists(path_hapmap)){
		cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Detect file ") %+% red$bgWhite$bold(path_hapmap) %+% black$bgWhite$bold(" is not exist or its path is wrong.") %+% "\n")
		status <- 1
	}else{
		system(paste("python ./data_transfer/hapmap2vcf.py -i ", path_hapmap, " -o ", path_datafolder, "/snphub_convert_output.vcf", sep=""))
		path_vcfolder <- path_datafolder
	}

}

if(status != 1){ setup <- file_check(status) }

if(status != 1){ setup <- file_make(status) }

if(status != 1){ file_recheck() }

system(paste("/bin/bash ./setup_scripts/ClearTmpFiles.sh ", path_datafolder, sep=""))