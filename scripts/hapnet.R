hn_main <- function(co, ro, anno, ext) {
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
        if(FALSE %in% (hn_sam == pub_unique_sample(hn_sam))){ return(hn_error_message(5)) }

        fet_code <- hn_fetch_data(hn_chr, hn_beg, hn_end, hn_sam)

        if(fet_code != 0) { return(hn_error_message(fet_code + 100)) }

        hap_chk = hn_check_hap()

        if(length(hap_chk)==1 && !is.na(as.numeric(hap_chk))){ return(hn_error_message(201)) }

        p <- hn_draw_plot(hn_gro, co, ro, ext, anno)

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

hn_check_hap <- function(){
    dna <- vcfR2DNAbin(vcfR_hn_orivcf, consensus = T, extract.haps = F)
    hap <- haplotype(dna)
    hap <- sort(hap, what = "labels")
    withCallingHandlers(
        warning = function(cnd){
            return(-1)
        },  
        dna
    )
}

hn_draw_plot <- function(groupf, co, ro, ext, anno){
    meta_data <- fra_glo_metadata_raw
    LIST <- meta_data[,3]
    names(LIST) <- meta_data[,1]

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
    ## layout
    layout(mat = matrix(c(2,1,3),ncol = 1), heights = c(1,5,1))

    # main plot
    par(mar=c(0,3,0,3))
    color <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD")
    pal <- color[1:ncol(hap.pies)]
    plotObj <- plot(hap.net, size=attr(hap.net, "freq")^0.5, bg=pal, pie=hap.pies, scale.ratio = 2, fast = F, show.mutation = 3)
    if(as.character(anno) == "Yes"){
        # xleft <- min(plotObj$xx)
        # xright <- max(plotObj$xx)
        # yupper <- max(plotObj$yy)
        # ybuttom <- min(plotObj$yy)
        xleft <- par("usr")[1]
        xright <- par("usr")[2]
        ybuttom <- par("usr")[3]
        yupper <- par("usr")[4]

        xrange <- xright - xleft
        yrange <- yupper - ybuttom


        sample_order <- LIST[attr(dna, "dimnames")[[1]]]
        hap_order <- attr(hap, "dimnames")[[1]]

        if (xrange >= yrange){
            textx <- xleft + xrange/10
            texty <- yupper + (xrange - yrange)/2 - yrange/3
        } else {
            textx <- xleft - (yrange - xrange)/2 + xrange/10
            texty <- yupper - yrange/3
        }

        #toprint <- ""
        offl <- 0
        for (i in seq_along(hap_order)){
            tmpname <- paste0("Samples in haplotype ", hap_order[i], ":", paste(sample_order[attr(hap, "index")[[i]]], collapse = ", "))
            #toprint <- paste0(toprint, tmpname, "\n")
            toprint <- stringi::stri_wrap(tmpname, exdent = 4, width = 130)
            text(x = textx, y = texty - offl, labels = paste(toprint, collapse = "\n"), pos = 4)
            offl <- offl + length(toprint)*yrange*0.020
        }
        # text(x = textx, y = texty, labels = stringi::stri_wrap(toprint, exdent = 2, width = 150, normalize = F), pos = 4)
    }

    ## plot2: title and legend
    par( mar = c(0,3,3,3) )
    plot(x=0, type="n", bty="n", xaxt="n", yaxt="n", 
        xlab="", ylab="", 
        xlim=c(0, 1), ylim=c(0, 1), xaxs="i", yaxs="i")
    categories <- colnames(hap.pies)
    text(0.5, 1, cex = 3, "HapNet", pos = 1)
    legend("bottomright", fill = pal, legend = categories, border = F, box.col = "white", horiz=TRUE, cex = 1.5, x.intersp = 0.5)

    ## plot3: tag and pars
    # ========================
    # plot3: tag and paras
    #
    t = read.table(pipe("date +%Y%m%d%H%M%S"))
    d = as.character(t[1,1])
    filenametag <- paste("SnpHub_HapNet_", d, sep="")
    parameter <- paste("Parameter: ", co, "; ", ro,"; flanking ", ext, ";", sep="")
    #
    par( mar = c(3,3,0,3) )
    plot(NULL, NULL, type="n", bty="n", xaxt="n", yaxt="n", 
         xlab="", ylab="", 
         xlim=c(0, 1), ylim=c(0, 1), xaxs="i", yaxs="i")
    text(x = 1, y = 0.5, paste0(filenametag, "\n", parameter), pos = 2)
}

hn_range_is_too_long <- function(hn_beg, hn_end) { F }


hn_error_message <- function(code){
    code
}