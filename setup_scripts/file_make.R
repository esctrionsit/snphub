file_make <- function(status){
	cat(blue$bgWhite$bold("[ ] Begining mew file making") %+% "\n")

	tryCatch(
		{
			system(paste("/bin/bash setup_scripts/IndexFaGff.sh", path_datafolder, path_samtools, name_fasta, path_tabix, name_gff3, sep=" "))

			system(paste("/bin/bash setup_scripts/GenomeDatabaseBuild.sh", path_datafolder, name_gff3, name_fasta, sep=" "))

			system(paste("/bin/bash setup_scripts/AnnoWithSnpEff.sh", path_datafolder, path_vcfolder, name_gff3, path_bcftools, sep=" "))

			system(paste("/bin/bash setup_scripts/MakeGeneIndex.sh", path_datafolder, name_gff3, sep=" "))

			system(paste("/bin/bash setup_scripts/ClearTmpFiles.sh", path_datafolder, sep=" "))
		}, error = function(e) {
			cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Some things are wrong during file building. Please check the log for more details.") %+% "\n"
				%+% black$bgWhite$bold("    ") %+% black$bgWhite$bold(e) %+% "\n"
				)
			status <- 1
		}
	)
	cat(green$bgWhite$bold("[+] new file making finished") %+% "\n")
	#return
	status
}