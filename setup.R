source("setup.conf")
source("advanced_config.R")
source("./setup_scripts/envir_check.R")
source("./setup_scripts/file_check.R")
source("./setup_scripts/file_make.R")
source("./setup_scripts/file_recheck.R")

status <- 0

envir_check(status)

if(status != 1){ setup <- file_check(status) }

if(status != 1){ setup <- file_make(status) }

if(status != 1){ file_recheck() }

system(paste("/bin/bash ./setup_scripts/ClearTmpFiles.sh ", path_datafolder, sep=""))