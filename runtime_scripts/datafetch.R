dt_fetch_data <- function(dt_chr, dt_beg, dt_end, dt_sam) {
    shell <- paste(path_bcftools, " view ", sep = "")
    shell <- paste(shell, path_vcf, " -s ", paste(dt_sam, collapse=","), " --min-ac=1 --threads 4 -m2 -M2 -r ", sep = "")
    for(i in 1:length(dt_chr)){
        if(i == length(dt_chr)){
            shell <- paste(shell, dt_chr[i], ":", dt_beg[i], "-", dt_end[i], " ", sep = "")
        }else{
            shell <- paste(shell, dt_chr[i], ":", dt_beg[i], "-", dt_end[i], ",", sep = "")
        }
    }
    rand <- read.table(pipe("cat /dev/urandom | head -n 10 | md5sum | head -c 10"))
    rand <- as.character(rand[1,1])

    shell <- paste(shell, " > ", getwd(), "/tmp/tmp.", rand, sep = "")

    print(shell)
    
    system(shell)

    code = tryCatch({
        tmp_data <- read.table(paste(getwd(), "/tmp/tmp.", rand, sep = ""), comment.char = "#")
        0
    }, error = function(e) {
        return(1)
    })

    if(code == 0){
        vcfR_dt_orivcf <<- read.vcfR(paste(getwd(), "/tmp/tmp.", rand, sep = ""))
        system(paste("rm -f ", getwd(), "/tmp/tmp.", rand, sep = ""))
        #return
        0
    }else{
        vcfR_dt_orivcf <<- data.frame()
        system(paste("rm -f ", getwd(), "/tmp/tmp.", rand, sep = ""))
        #return
        code
    }    
    
}

hp_fetch_data <- function(hp_chr, hp_beg, hp_end, hp_sam, maf) {
    shell <- paste(path_bcftools, " view ", path_vcf, sep = "")
    if(maf != "0"){
        shell <- paste(shell, " -e 'MAF<", maf, "' ", sep="")
    }
    shell <- paste(shell, " -r ", sep = "")
    for(i in 1:length(hp_chr)){
        if(i == length(hp_chr)){
            shell <- paste(shell, hp_chr[i], ":", hp_beg[i], "-", hp_end[i], " ", sep = "")
        }else{
            shell <- paste(shell, hp_chr[i], ":", hp_beg[i], "-", hp_end[i], ",", sep = "")
        }
    }
    shell <- paste(shell, " -s ", paste(hp_sam, collapse=","), sep = "")
    shell <- paste(shell, " | ", path_bcftools, " query ", sep="")
    shell <- paste(shell, "-f '%CHROM\\t%POS\\t%ANN[\\t%GT]\\n' ", sep = "")

    bcftools_leng <- read.table(pipe(paste(shell, " | wc -l")), header = F, comment.char = "", as.is = T)
    if(bcftools_leng[1,1] != "0"){
        fra_hp_orivcf <<- read.table(pipe(shell), header = F, comment.char = "#", as.is = T)
        names(fra_hp_orivcf) <<- c("CHROM", "POS","ANN", hp_sam)
    }else{
        fra_hp_orivcf <<- data.frame()
    }
    #return
    0
}

hm_fetch_data <- function(chr, beg, end, sam){
    shell <- paste(path_bcftools, " view ", path_vcf, sep = "")
    shell <- paste(shell, " -r ", sep = "")
    for(i in 1:length(chr)){
        if(i == length(chr)){
            shell <- paste(shell, chr[i], ":", beg[i], "-", end[i], " ", sep = "")
        }else{
            shell <- paste(shell, chr[i], ":", beg[i], "-", end[i], ",", sep = "")
        }
    }
    shell <- paste(shell, " -s ", paste(sam, collapse=","), sep = "")
    shell <- paste(shell, " | ", path_bcftools, " query ", sep="")
    shell <- paste(shell, "-f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' ", sep = "")

    bcftools_leng <- read.table(pipe(paste(shell, " | wc -l")), header = F, comment.char = "", as.is = T)
    if(bcftools_leng[1,1] != "0"){
        fra_hm_orivcf <<- read.table(pipe(shell), header = F, comment.char = "#", as.is = T)
        names(fra_hm_orivcf) <<- c("CHROM", "POS", "REF", "ALT", sam)
    }else{
        fra_hm_orivcf <<- data.frame()
        return(1)
    }
    #return
    0
}

hn_fetch_data <- function(hn_chr, hn_beg, hn_end, hn_sam) {
    shell <- paste(path_bcftools, " view ", sep = "")
    shell <- paste(shell, path_vcf, " -s ", paste(hn_sam, collapse=","), " --min-ac=1 --threads 4 -m2 -M2 -v snps -r ", sep = "")
    for(i in 1:length(hn_chr)){
        if(i == length(hn_chr)){
            shell <- paste(shell, hn_chr[i], ":", hn_beg[i], "-", hn_end[i], " ", sep = "")
        }else{
            shell <- paste(shell, hn_chr[i], ":", hn_beg[i], "-", hn_end[i], ",", sep = "")
        }
    }
    shell <- paste(shell, " | sed 's%1/1%1|1%g; s%0/0%0|0%g; s%./.%.|.%g' ", sep = "")
    rand <- read.table(pipe("cat /dev/urandom | head -n 10 | md5sum | head -c 10"))
    rand <- as.character(rand[1,1])

    shell <- paste(shell, " > ", getwd(), "/tmp/tmp.", rand, sep = "")
    
    system(shell)

    code = tryCatch({
        tmp_data <- read.table(paste(getwd(), "/tmp/tmp.", rand, sep = ""), comment.char = "#")
        0
    }, error = function(e) {
        return(1)
    })

    if(code == 0){
        vcfR_hn_orivcf <<- read.vcfR(paste(getwd(), "/tmp/tmp.", rand, sep = ""))
        system(paste("rm -f ", getwd(), "/tmp/tmp.", rand, sep = ""))
        #return
        0
    }else{
        vcfR_hn_orivcf <<- data.frame()
        system(paste("rm -f ", getwd(), "/tmp/tmp.", rand, sep = ""))
        #return
        code
    }    
}

lp_fetch_data <- function(chr, beg, end, sam){
    shell <- paste(path_bcftools, " view ", path_vcf, sep = "")
    shell <- paste(shell, " -r ", sep = "")
    for(i in 1:length(chr)){
        if(i == length(chr)){
            shell <- paste(shell, chr[i], ":", beg[i], "-", end[i], " ", sep = "")
        }else{
            shell <- paste(shell, chr[i], ":", beg[i], "-", end[i], ",", sep = "")
        }
    }
    shell <- paste(shell, " -s ", paste(sam, collapse=","), sep = "")
    shell <- paste(shell, " | ", path_bcftools, " query ", sep="")
    shell <- paste(shell, "-f '%CHROM\\t%POS\\t%AC\\t%AN\\t%ANN\\t[\\t%GT]\\n' ", sep = "")
    print(shell)
    bcftools_leng <- read.table(pipe(paste(shell, " | wc -l")), header = F, comment.char = "", as.is = T)
    if(bcftools_leng[1,1] != "0"){
        fra_lp_orivcf <<- read.table(pipe(shell), header = F, comment.char = "#", as.is = T)
        names(fra_lp_orivcf) <<- c("CHROM", "POS", "AC", "AN", "ANN", sam)
    }else{
        fra_lp_orivcf <<- data.frame()
        return(1)
    }
    #return
    0
}

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

    print(shell)

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
