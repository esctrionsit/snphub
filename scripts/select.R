snp_main <- function(ty, oi, co_t, co_f, co_e, ro, ro_ext, maf ="0" , bso= "No", mlr = "1"){
    if(is.null(bso)){
         maf = "0"
         bso = "No"
         mlr = "1"
    }

    text_snp_currpara <<- ""

	sel_chr <- pub_check_chr_value(ro)

    if(is.na(as.numeric(ro_ext)) || as.numeric(ro_ext) < 0) { return(snp_error_message(1)) }

    sel_beg <- pub_check_pos_begin(ro, ro_ext)
    sel_end <- pub_check_pos_end(ro, ro_ext)

    if(sum(is.na(sel_chr), is.na(sel_beg), is.na(sel_end)) > 0 || length(sel_chr) == 0) { return(snp_error_message(2)) }
    if(snp_range_is_too_long(sel_beg, sel_end, ty)) { return(snp_error_message(3)) }
    if(paste(co_t, co_f, co_e, sep = "") == "") {return(snp_error_message(10))}
    
    sel_sam <- pub_check_sample_name(c(co_t, co_f, co_e))

    if(sum(is.na(sel_sam)) > 0 || length(sel_sam) == 0) { return(snp_error_message(4)) }
    if(sample_is_unique(sel_sam) == F) {return(snp_error_message(8))}
    if(is.na(as.numeric(maf))) { return(snp_error_message(5)) }
    if(is.na(as.numeric(mlr))) { return(snp_error_message(9)) }

    shell <- snp_fetch_data(oi, sel_chr, sel_beg, sel_end, ty, sel_sam, maf, bso, mlr)
    #return(data.frame(c(shell)))

    if(nrow(fra_snp_orivcf) == 0) { return(snp_error_message(6)) }

    sel_rule <- snp_get_sam_rule(co_t, co_f, co_e)
    sel_code <- snp_core_select(sel_rule,oi)

    if(sel_code > 0) { return(snp_error_message(sel_code + 100)) }
    if(nrow(fra_snp_finalres) == 0) { return(snp_error_message(7)) }

    snp_finalres_header_mapping()

    text_snp_currpara <<- paste("Parameter: ", ty, " ; ", ro, " ; flask length ", ro_ext, sep="")

    fra_snp_finalres
}


snp_range_is_too_long <- function(sel_beg, sel_end, ty) { F }

snp_fetch_data <- function(oi, sel_chr, sel_beg, sel_end, ty, fet_sam, maf="0", bso="No", mlr="1") {
    shell <- paste(path_bcftools, " view ", path_vcf, sep = "")
    if(ty == "snp"){
        shell <- paste(shell, " -v snps ", sep = "")
    }else if(ty == "indel"){
        shell <- paste(shell, " -v indels ", sep = "")
    }
    if(bso == "Yes"){
        shell <- paste(shell, " -M 2 ", sep = "")
    }
    shell <- paste(shell, " -r ", sep = "")
    for(i in 1:length(sel_chr)){
        if(i == length(sel_chr)){
            shell <- paste(shell, sel_chr[i], ":", sel_beg[i], "-", sel_end[i], " ", sep = "")
        }else{
            shell <- paste(shell, sel_chr[i], ":", sel_beg[i], "-", sel_end[i], ",", sep = "")
        }
    }
    shell <- paste(shell, " -s ", paste(fet_sam, collapse=","), sep = "")
    if(maf != "0"){
        shell <- paste(shell, " -e 'MAF<", maf, "' ", sep="")
    }
    if(mlr != "1"){
        shell <- paste(shell, " -e 'F_MISSING>", mlr, "' ", sep="")
    }
    shell <- paste(shell, " | ", path_bcftools, " query ", sep="")
    shell <- paste(shell, "-f '%CHROM\\t%POS", sep="")
    if("ANN" %in% oi){
        shell <- paste(shell, "\\t%ANN", sep="")
    }
    shell <- paste(shell, "\\t%REF\\t%ALT[\\t%GT", sep="")
    if("DP" %in% oi){
        shell <- paste(shell, ":%DP", sep="")
    }
    if("GQ" %in% oi){
        shell <- paste(shell, ":%GQ", sep="")
    }
    shell <- paste(shell, "]\\n' ", sep="")

    bcftools_leng <- read.table(pipe(paste(shell, " | wc -l")), header = F, comment.char = "", as.is = T)
    if(bcftools_leng[1,1] != "0"){
        fra_snp_orivcf <<- read.table(pipe(shell), header = F, comment.char = "#", as.is = T)
        if("ANN" %in% oi){
            names(fra_snp_orivcf) <<- c("CHROM", "POS", "ANN", "REF", "ALT", fet_sam)
        }else{
            names(fra_snp_orivcf) <<- c("CHROM", "POS", "REF", "ALT", fet_sam)
        }
    }else{
        fra_snp_orivcf <<- data.frame()
    }
    #return
    shell
}

snp_get_sam_rule <- function(co_t, co_f, co_e) {
    rule <- c()
    rule <- c(rule, rep("+", time=length(strsplit(co_t, split = ",")[[1]])))
    rule <- c(rule, rep("-", time=length(strsplit(co_f, split = ",")[[1]])))
    rule <- c(rule, rep(" ", time=length(strsplit(co_e, split = ",")[[1]])))
    #return
    rule
}

sample_is_unique <- function(sel_sam){
    for(i in 1:length(sel_sam)){
        if(sel_sam[i] %in% sel_sam[-i]){
            return(F)
        }
    }
    T
}

snp_core_select <- function(rule, oi, s_b=5) {
    swi <- T
    res <- c()
    if("ANN" %in% oi){
        s_b <- 6
    }
    for(i in 1:nrow(fra_snp_orivcf)){
        swi <- T
        for(j in s_b:ncol(fra_snp_orivcf)){
            tmp <- fra_snp_orivcf[i,j]
            tmp <- substr(tmp,1,3)
            if(rule[j-s_b+1] == "+"){
                if(tmp == "0/0" || tmp == "./."){
                    swi <- F
                    break
                }
            }else if(rule[j-s_b+1] == "-"){
                if(tmp != "0/0" && tmp != "./."){
                    swi <- F
                    break
                }
            }else if(rule[j-s_b+1] == " "){
                break
            }
        }
        if(swi){
            res <- c(res, as.numeric(fra_snp_orivcf[i,2]))
        }
    }    
    fra_snp_finalres <<- subset(fra_snp_orivcf, POS %in% res)

    #return
    0
}

snp_finalres_header_mapping <- function(){
    tmp <- names(fra_snp_finalres)
    for(i in 1:length(tmp)){
        tmp2 <- as.character(fra_glo_metadata[which(fra_glo_metadata$Accession == tmp[i]),]$Label)
        if(length(tmp2) != 1){
            next
        }else{
            tmp[i] <- tmp2
        }
    }
    names(fra_snp_finalres) <<- tmp
    #return
    0
}

snp_error_message <- function(code){
    if(code == 1){
        data.frame(Info=c("Error 0001:Flanking length must be integer."))
    }else if(code == 2){
        data.frame(Info=c("Error 0002:Invalid region."))
    }else if(code == 3){
        data.frame(Info=c("Error 0003:Region is too long."))
    }else if(code == 4){
        data.frame(Info=c("Error 0004:Invalid accession detected."))
    }else if(code == 5){
        data.frame(Info=c("Error 0005:MAF must be numeric."))
    }else if(code == 6){
        data.frame(Info=c("Error 0006:No variation found in current regions."))
    }else if(code == 7){
        data.frame(Info=c("Error 0007:No variation found under current filters."))
    }else if(code == 8){
        data.frame(Info=c("Error 0008:Duplicate sample name is not allowed."))
    }else if(code == 9){
        data.frame(Info=c("Error 0009:Maxium missing rate must be numeric"))
    }else if(code == 10){
        data.frame(Info=c("Error 0010:No sample name found."))
    }
}