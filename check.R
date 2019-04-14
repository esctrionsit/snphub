source("./config.R")

###############################################################################################
# Package checking
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

if(!is.installed("crayon")){
	warning("Detect package \"crayon\" is not installed in your R enviroment.")
	warning("Trying to install the \"crayon\" package.")
	warning("If failed, please try to install it mannually.")
	install.packages("crayon")
}
suppressPackageStartupMessages(library("crayon"))

packages_check_list <- c("RColorBrewer", "ggplot2", "ggtree", "shiny", "pegas", "vcfR", "ape", "DT")

for(i in packages_check_list){
	if(!is.installed(i)){
		cat(red("[-]") %+% " Detect package " %+% red(i) %+% " is not installed in your R enviroment.\n"
			%+% "    Trying to install the " %+% red(i) %+%" package.\n"
			%+% "    If failed, please try to install it mannually.\n"
			)
		install.packages(i)
	}
}
###############################################################################################
# Config checking
tool_check_list <- c(path_bcftools, path_samtools)
for(i in tool_check_list){
	tryCatch(
		{
			t <- read.table(pipe(paste(i, "--version", sep = " ")), sep = "~")
		}, error = function(e) {
			cat(red("[-]") %+% " Detect file " %+% red(i) %+% " is not exist.\n")
		}
	)
}

tryCatch(
	{
		t <- read.table(pipe(path_seqkit), sep = "~")
	}, error = function(e) {
		cat(red("[-]") %+% " Detect file " %+% red(path_seqkit) %+% " is not exist.\n")
	}
)


file_sub_list <- c(path_fasta, path_fa_index, path_vcf)
for(i in file_sub_list){
	if(!file.exists(i)){
		cat(red("[-]") %+% " Detect file " %+% red(i) %+% " is not exist.\n")
	}
}

# gene index
is_correct <- tryCatch(
	{
		fra_glo_geneindex <- read.table(path_geneindex, col.names = c("CHROM", "begin", "end", "gene_name"), header = F, as.is = T)
		TRUE
	}, error = function(e) {
		cat(red("[-]") %+% " Detect file " %+% path_geneindex %+% " is not correct.\n"
			%+% "    " %+% e
			)
		FALSE
	}
)
if(is_correct){
	if(sum(is.na(as.numeric(fra_glo_geneindex$begin))) > 0){
		cat(red("[-]") %+% " Detect file " %+% path_geneindex %+% " is not correct.\n"
			%+% "    " %+% "The second column should be numeric (beginning position of gene).\n"
			)
	}
	if(sum(is.na(as.numeric(fra_glo_geneindex$end))) > 0){
		cat(red("[-]") %+% " Detect file " %+% path_geneindex %+% " is not correct.\n"
			%+% "    " %+% "The third column should be numeric (ending position of gene).\n"
			)
	}
	if(length(unique(fra_glo_geneindex$gene_name)) != length(fra_glo_geneindex$gene_name)){
		cat(red("[-]") %+% " Detect file " %+% path_geneindex %+% " is not correct.\n"
			%+% "    " %+% "Gene name in forth column should be unique.\n"
			)
	}
}
# metadata
is_correct <- tryCatch(
	{
		fra_glo_metadata <- read.table(path_metadata, col.names = c("Accession", "Label", "Name"), header = F, as.is = T)
		TRUE
	}, error = function(e) {
		cat(red("[-]") %+% " Detect file " %+% path_metadata %+% " is not correct.\n"
			%+% "    " %+% e
			)
		FALSE
	}
)
if(is_correct){
	if(length(unique(fra_glo_metadata$Accession)) != length(fra_glo_metadata$Accession)){
		cat(red("[-]") %+% " Detect file " %+% path_metadata %+% " is not correct.\n"
			%+% "    " %+% "Accession in first column should be unique.\n"
			)
	}
	if(length(unique(fra_glo_metadata$Label)) != length(fra_glo_metadata$Label)){
		cat(red("[-]") %+% " Detect file " %+% path_metadata %+% " is not correct.\n"
			%+% "    " %+% "Label in second column should be unique.\n"
			)
	}else{
		tryCatch(
			{
				t <- read.table(pipe(paste(path_bcftools, "view", path_vcf, "-h", "-r", paste(fra_glo_metadata$Accession, collapse = ","), sep=" ")), comment.char = "", sep = "~")
			}, error = function(e) {
				cat(red("[-]") %+% " Detect file " %+% path_groupdata %+% " is not correct.\n"
					%+% "    " %+% "Accession in first column should be in the vcf/bcf file.\n"
					)
			}
		)
	}
}
# groupdata
is_correct <- tryCatch(
	{
		fra_glo_groupdata <- read.table(path_groupdata, col.names = c("Group", "Name"), header = F, as.is = T)
		TRUE
	}, error = function(e) {
		cat(red("[-]") %+% " Detect file " %+% path_groupdata %+% " is not correct.\n"
			%+% "    " %+% e %+% "\n"
			)
		FALSE
	}
)
if(is_correct){
	if(length(unique(fra_glo_groupdata$Group)) != length(fra_glo_groupdata$Group)){
		cat(red("[-]") %+% " Detect file " %+% path_groupdata %+% " is not correct.\n"
			%+% "    " %+% "Group name in first column should be unique.\n"
			)
	}
	for(i in 1:nrow(fra_glo_groupdata)){
		s <- fra_glo_groupdata[i, 2]
		s <- strsplit(s, split = ",")[[1]]
		for(j in s){
			if(!(j %in% fra_glo_metadata$Label)){
				cat(red("[-]") %+% " Detect file " %+% path_groupdata %+% " is not correct.\n"
					%+% "    " %+% "Wrong sample label in group " %+% red(fra_glo_groupdata[i, 1]) %+% "\n" 
					)
			}
			break
		}
	}
}

###############################################################################################
# Finish
cat(green("Finish\n"))