ns_main <- function(ns_ty1, ns_ty2, ns_co, ns_ro){
    ns_co <- paste(strsplit(ns_co, split = " ")[[1]], collapse="")
	fra_ns_finalres <<- data.frame()
	ns_chr <- pub_check_chr_value(ns_ro)
	ns_beg <- pub_check_pos_begin(ns_ro)
    ns_end <- pub_check_pos_end(ns_ro)
    ns_sam <- pub_check_sample_name(ns_co, T)

    if(sum(is.na(ns_chr), is.na(ns_beg), is.na(ns_end)) > 0 || length(ns_chr) == 0) { return(ns_error_message(1)) }
    if(sum(is.na(ns_sam)) > 0 || length(ns_sam) == 0) { return(ns_error_message(2)) }

    ns_chr <- ns_chr[1]
	ns_beg <- ns_beg[1]
    ns_end <- ns_end[1]

    ns_readdata(ns_ty1, ns_ty2, ns_sam, ns_chr, paste(ns_beg, ns_end, sep = "-"))
    names(fra_ns_finalres) <- c(paste(ns_ty1, ns_ty2, ns_co, ns_ro, sep=" | "))
    #return
    fra_ns_finalres
}

ns_error_message <- function(code){
	if(code == 1){
        data.frame(Info=c("Error 0001:Invalid region."))
    }else if(code == 2){
        data.frame(Info=c("Error 0002:Invalid accession detected."))
    }
}

ns_readdata <- function(ty1, ty2, samp, ns_ge, ro){
    for(i in samp) {
        shell <- paste(path_samtools, " faidx ", path_fasta," ", ns_ge, ":", ro, sep = "")
        if(i != "#RAW"){
            shell <- paste(shell, " |", path_bcftools, " consensus ", path_vcf, " -s ", i, sep = "")
            if(ty1 != "both"){
                shell <- paste(shell, " -e \'TYPE!=\"", ty1, "\"\' ", sep = "")
            }
            if(ty2 != "all") { shell <- paste(shell, " -H R ", sep = "") } 
            shell <- paste(shell, " |", path_seqkit, " seq", sep = "")
        }
        #fra_ns_finalres <<- rbind(fra_ns_finalres, data.frame(c(shell)))
       # next
        strdata <- read.table(pipe(shell), as.is = T)
        head <- strdata[1,1]
        head <- strsplit(head, ">")[[1]]
        if(i == "#RAW"){
            tmp_samp <- "RAW"
        }else{
            tmp_samp <- as.character(fra_glo_metadata[which(fra_glo_metadata$Accession == i),]$Label)
        }
        head <- paste(">", tmp_samp, ":", head[2], sep="")
        strdata[1,1] <- head
        fra_ns_finalres <<- rbind(fra_ns_finalres, strdata)
    }
    names(fra_ns_finalres) <<- c("")
}
