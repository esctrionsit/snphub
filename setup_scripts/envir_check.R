envir_check <- function(status){

	###############################################################################################
	# Package checking
	is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

	if(!is.installed("crayon")){
		stop("Detect package \"crayon\" is not installed in your R enviroment.\nPlease use \ninstall.packages(\"crayon\")\nto install it FIRST.")
	}
	suppressPackageStartupMessages(library("crayon"))

	packages_check_list <- c("ggplot2", "ggmap", "shiny", "dplyr", "rjson", "pegas", "vcfR", "ape", "DT")

	for(i in packages_check_list){
		if(!is.installed(i)){
			cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Detect package ") %+% red$bgWhite$bold(i) %+% black$bgWhite$bold(" is not installed in your R enviroment.") %+% "\n"
				%+% black$bgWhite$bold("    Maybe you can use ") %+% red$bgWhite$bold("install.packages(\"") %+% red$bgWhite$bold(i) %+% red$bgWhite$bold("\")") %+% black$bgWhite$bold(" to install it.") %+% "\n"
				)
			status <- 1
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
				cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Detect program ") %+% red$bgWhite$bold(i) %+% black$bgWhite$bold(" is not exist.") %+% "\n")
				status <- 1
			}
		)
	}

	tryCatch(
		{
			t <- read.table(pipe(path_seqkit), sep = "~")
		}, error = function(e) {
			cat(red$bgWhite$bold("[-]") %+% black$bgWhite$bold(" Detect program ") %+% red$bgWhite$bold(path_seqkit) %+% black$bgWhite$bold(" is not exist.") %+% "\n")
			status <- 1
		}
	)

	###############################################################################################
	# Finish
	cat(green$bgWhite$bold("[+] Environment Checking Finished") %+% "\n")

	#return
	status
}