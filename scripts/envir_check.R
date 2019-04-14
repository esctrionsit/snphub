source("../config.R")

###############################################################################################
# Package checking
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

if(!is.installed("crayon")){
	stop("Detect package \"crayon\" is not installed in your R enviroment.\nPlease use \ninstall.packages(\"crayon\")\nto install it FIRST.")
}
suppressPackageStartupMessages(library("crayon"))

packages_check_list <- c("RColorBrewer", "ggplot2", "ggtree", "shiny", "pegas", "vcfR", "ape", "DT")

for(i in packages_check_list){
	if(!is.installed(i)){
		cat(red("[-]") %+% " Detect package " %+% red(i) %+% " is not installed in your R enviroment.\n"
			%+% "    Please use " %+% red("install.packages(\"" %+% i %+% "\")") %+%" to install it.\n"
			)
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

###############################################################################################
# Finish
cat(green("[ ] Environment Checking Finished\n"))