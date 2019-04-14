hn_main <- function(co, ro, ext) {
    withProgress(message = 'Drawing', detail = "  Drawing Haplotype Network plot...", value = 5, {
        text_hn_currpara <<- ""

        hn_chr <- pub_check_chr_value(ro)

        if(is.na(as.numeric(ext)) || as.numeric(ext) < 0) { return(hn_error_message(1)) }

        hn_beg <- pub_check_pos_begin(ro, ext)
        hn_end <- pub_check_pos_end(ro, ext)

        if(sum(is.na(hn_chr), is.na(hn_beg), is.na(hn_end)) > 0 || length(hn_chr) == 0) { return(hn_error_message(2)) }
        if(hn_range_is_too_long(hn_beg, hn_end)) { return(hn_error_message(3)) }

        hn_sam <- hn_check_sample_name(co)
        hn_gro <- pub_sub_group(co)

        if(sum(is.na(hn_sam)) > 0 || length(hn_sam) == 0) { return(hn_error_message(4)) }
        if(sum(is.na(hn_gro)) > 0 || nrow(hn_gro) == 0) { return(hn_error_message(6)) }

        fet_code <- hn_fetch_data(hn_chr, hn_beg, hn_end, hn_sam)

        if(fet_code != 0) { return(hn_error_message(fet_code + 100)) }

        p <- hn_draw_plot(hn_gro, co, ro, ext)

        text_hn_currpara <<- paste("Parameter: ", ro, " ; flask length ", ext, sep="")

        p
    })
}

hn_check_sample_name <- function(co) {
    res <- c()
    if(grepl("\\}", co)){
        tmp1 <- strsplit(co, split = "\\}")[[1]]
        for(i in tmp1){
            samps <- strsplit(i, split = "\\{")[[1]][2]
            samps <- strsplit(samps, split = ",")[[1]]
            for(j in samps){
                tmp2 <- as.character(fra_glo_metadata[which(fra_glo_metadata$Label == j),]$Accession)
                if(length(tmp2) != 1) { tmp2 <- NA }
                res <- c(res, tmp2)
            }
        }
    }else{
        tmp1 <- strsplit(co, split = ",")[[1]]
        for(i in tmp1){
            tmp2 <- as.character(fra_glo_groupdata[which(fra_glo_groupdata$Group == i),]$Name)
            if(length(tmp2) != 1) {
                tmp2 <- NA
            }else{
                tmp2 <- strsplit(tmp2, split = ",")[[1]]
            }
            for(j in tmp2){
                tmp3 <- as.character(fra_glo_metadata[which(fra_glo_metadata$Label == j),]$Accession)
                if(length(tmp3) != 1) { tmp3 <- NA }
                res <- c(res, tmp3)
            }
        }
    }
    #return
    res
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

hn_draw_plot <- function(groupf, co, ro, ext){
    dna <- vcfR2DNAbin(vcfR_hn_orivcf, consensus = T, extract.haps = F)
    hap <- haplotype(dna)
    hap <- sort(hap, what = "labels")
    hap.net <- haploNet(hap, getProb = FALSE)

    # group file
    groupf <- groupf[match(labels(dna), groupf$Sample),]

    # do some magic with the row names and categories that I figured out a while ago 
    # but no longer entirely understand how it works. Essentially, this gives a table
    # with haplotypes and their frequencies by category so you can do the color-coded
    # pies for each haplotype in the plot.
    hap.pies <- with(
      stack(setNames(attr(hap,'index'),1:length(attr(hap,'index')))),
      table(hap=as.numeric(as.character(ind)),groupf[values,"Group"])
    )
    rownames(hap.pies) <- rownames(hap)
    pal <- brewer.pal(ncol(hap.pies), "Set3")
    plot(hap.net, size=attr(hap.net, "freq")^0.5, bg=pal, pie=hap.pies, scale.ratio = 2, fast = F, show.mutation = 3)
    mtext(co, line = 2, adj = 0)
    mtext(paste(ro, " with extension ", ext, sep = ""), line = 1, adj = 0)
    #plot(hap.net, bg=pal, pie=hap.pies, scale.ratio = 1)
    categories <- colnames(hap.pies)
    legend("topright", legend=categories,fill=pal,bty="n",cex=1.2,ncol=1)
}

hn_range_is_too_long <- function(hn_beg, hn_end) { F }


hn_error_message <- function(code){
    code
}