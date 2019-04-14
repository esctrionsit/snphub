dt_main <- function(ty, ty_tree, dbr, co, ro, ext){
    withProgress(message = 'Drawing', detail = "  Drawing PhyloTree plot...", value = 5, {
        text_dt_currpara <<- ""

        dt_chr <- pub_check_chr_value(ro)

        if(is.na(as.numeric(ext)) || as.numeric(ext) < 0) { return(dt_error_message(1)) }

        dt_beg <- pub_check_pos_begin(ro, ext)
        dt_end <- pub_check_pos_end(ro, ext)

        if(sum(is.na(dt_chr), is.na(dt_beg), is.na(dt_end)) > 0 || length(dt_chr) == 0) { return(dt_error_message(2)) }
        if(dt_range_is_too_long(dt_beg, dt_end)) { return(dt_error_message(3)) }
        
        dt_sam <- dt_check_sample_name(co)
        dt_gro <- pub_sub_group(co)

        if(sum(is.na(dt_sam)) > 0 || length(dt_sam) == 0) { return(dt_error_message(4)) }
        if(sum(is.na(dt_gro)) > 0 || nrow(dt_gro) == 0) { return(dt_error_message(6)) }
        
        code <- dt_fetch_data(dt_chr, dt_beg, dt_end, dt_sam)

        if(code != 0) { return(code+100) }

        p <- dt_draw_plot(dt_gro, ty, ty_tree, dbr, co, ro, ext)

        if(ty == "MDS"){
            tmp <- "MDS"
        }else{
            tmp <- ty_tree
        }
        text_dt_currpara <<- paste("Parameter: ", ro, " ; flask length ", ext, " ; ", tmp, sep="")
        #return
        p
    })
}

dt_range_is_too_long <- function(dt_beg, dt_end) { F }

dt_check_sample_name <- function(co) {
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

dt_error_message <- function(code) {
    code
}

dt_draw_alpha <- function (colour, alpha = NA){
    col <- grDevices::col2rgb(colour, TRUE)/255
    if (length(colour) != length(alpha)) {
        if (length(colour) > 1 && length(alpha) > 1) {
            stop("Only one of colour and alpha can be vectorised")
        }
        if (length(colour) > 1) {
            alpha <- rep(alpha, length.out = length(colour))
        }else if (length(alpha) > 1) {
            col <- col[, rep(1, length(alpha)), drop = FALSE]
        }
    }
    alpha[is.na(alpha)] <- col[4, ][is.na(alpha)]
    new_col <- grDevices::rgb(col[1, ], col[2, ], col[3, ], alpha)
    new_col[is.na(colour)] <- NA
    new_col
}

dt_draw_plot <- function(group_info, dt_pl_ty, treetype, dbr, co, ro, ext) {
    if(dbr=="Yes"){
        edgelen = T
    }else{
        edgelen = F
    }
    names(group_info) <- c("Accession", "Group")
    dna <- vcfR2DNAbin(vcfR_dt_orivcf, consensus = T, extract.haps = F)
    distdf <- dist.dna(dna, model = "raw")
    LIST  <- fra_glo_metadata$Name
    names(LIST) <- fra_glo_metadata$Accession
    attr(distdf, "Labels") <- LIST[attr(distdf, "Labels")]
    info_df <- merge(merge(data.frame(Name = attr(distdf, "Labels")), fra_glo_metadata, by = "Name"), group_info, by = "Accession")
    names(info_df) <- c("Accession", "label", "Name", "Group")
    p <- tryCatch({
        if("NJ-tree"==dt_pl_ty){
            if (sum(is.na(distdf)) > 0){
                tree <- njs(distdf)
                bool_dt_warning <<- T
            } else {
                tree <- nj(distdf)
            }
            color <- rgb(c(67, 252, 141, 231, 166, 255, 229, 179), c(153, 141, 160, 138, 216, 217, 196, 179), c(125, 98, 203, 195, 84, 47, 148, 179), maxColorValue = 255)
            color <- color[1: length(unique(info_df$Group))]
            names(color) <- unique(info_df$Group)
            col.vector <- vector(mode="character",length=nrow(tree$edge))
            n.tips <- length(tree$tip.label)
            col.vector[tree$edge[,2]>n.tips] <- "#333333"
            edge.data <- as.data.frame(tree$edge)
            for(i in seq_along(tree$tip.label)){
              edge.row <- as.numeric(rownames(edge.data[edge.data$V2==i,]))
              col.vector[edge.row] <- color[info_df$Group[match(tree$tip.label[i], info_df$label)]]
            }
            tip.color <- color[info_df$Group]
            names(tip.color) <- as.character(info_df$label)
            layout(mat = matrix(c(2,1,3),ncol = 1), heights = c(1,5,1))

            # plot1: main plot
            par(mar=c(0,3,0,3))
            plot(tree, edge.color = col.vector, tip.col = tip.color[tree$tip.label], type = treetype, use.edge.length = edgelen)
            
            # ========================
            # plot2: header and legend
            #
            par( mar = c(0,3,3,3) )
            plot(x=0, type="n", bty="n", xaxt="n", yaxt="n", 
                 xlab="", ylab="", 
                 xlim=c(0, 1), ylim=c(0, 1), xaxs="i", yaxs="i")

            # Text ex:  "chr2: 90,936-91,653"
            # text( 0.5, 1, cex = 3,
            #      paste(chr, " : ", 
            #            format(LeftMost, big.mark = ",", scientific = F), " - ", 
            #            format(RightMost, big.mark = ",", scientific = F), sep = "") )
            text(0.5, 1, cex = 3, "Phylogenetic tree", pos = 1)

            # Draw the color legends
            legend("bottomright", fill = color, legend = names(color), border = F, box.col = "white", horiz=TRUE, cex = 1.5, x.intersp = 0.5)

            # ========================
            # plot3: tag and paras
            #
            t = read.table(pipe("date +%Y%m%d%H%M%S"))
            d = as.character(t[1,1])
            filenametag <- paste("SnpHub_PhyloTree_", d, sep="")
            parameter <- paste("Parameter: ", co, "; ", ro,"; flanking ", ext, ";", sep="")
            #
            par( mar = c(3,3,0,3) )
            plot(NULL, NULL, type="n", bty="n", xaxt="n", yaxt="n", 
                 xlab="", ylab="", 
                 xlim=c(0, 1), ylim=c(0, 1), xaxs="i", yaxs="i")
            text(x = 1, y = 0.5, paste0(filenametag, "\n", parameter), pos = 2)
        }else{
            # MDS
            if (sum(is.na(distdf)) > 0){
                return(201)
            }
            tree <- nj(distdf)
            loc <- cmdscale(distdf)
            loc_df <- merge(data.frame(label = row.names(loc), loc), info_df, by = "label")
            colnames(loc_df)[c(2,3)] <- c("Dim1", "Dim2")
            color <- rainbow(length(unique(info_df$Group)))
            names(color) <- unique(info_df$Group)
            #plot1: main plot
            par(mar=c(5,6,0,6))
            plot(loc_df$Dim1, loc_df$Dim2, type = "n", xlab = "Dim1", ylab = "Dim2", frame.plot = FALSE, cex.axis = 1.5, cex.lab = 2)
            text(loc_df$Dim1, loc_df$Dim2, col="gray", labels = loc_df$label)
            points(loc_df$Dim1, loc_df$Dim2, col=dt_draw_alpha(color[loc_df$Group],0.6), pch = 16)
            #========================
            # plot2: header and legend
            #
            par( mar = c(0,3,3,3) )
            plot(x=0, type="n", bty="n", xaxt="n", yaxt="n", 
                 xlab="", ylab="", 
                 xlim=c(0, 1), ylim=c(0, 1), xaxs="i", yaxs="i")

            # Text ex:  "chr2: 90,936-91,653"
            # text( 0.5, 1, cex = 3,
            #      paste(chr, " : ", 
            #            format(LeftMost, big.mark = ",", scientific = F), " - ", 
            #            format(RightMost, big.mark = ",", scientific = F), sep = "") )
            text(0.5, 1, cex = 3, "Multidimensional scaling", pos = 1)

            # Draw the color legends
            legend("bottomright", fill = color, legend = names(color), border = F, box.col = "white", horiz=TRUE, cex = 1.5, x.intersp = 0.5)

            #
            # ========================
            # plot3: tag and paras
            #
            t = read.table(pipe("date +%Y%m%d%H%M%S"))
            d = as.character(t[1,1])
            filenametag <- paste("SnpHub_PhyloTree_", d, sep="")
            parameter <- paste("Parameter: ", co, "; ", ro,"; flanking ", ext, ";", sep="")
            #
            par( mar = c(3,3,0,3) )
            plot(NULL, NULL, type="n", bty="n", xaxt="n", yaxt="n", 
                 xlab="", ylab="", 
                 xlim=c(0, 1), ylim=c(0, 1), xaxs="i", yaxs="i")
            text(x = 1, y = 0.5, paste0(filenametag, "\n", parameter), pos = 2)
        }
    }, error = function(e) {
        paste("202: Unknown draw time errror." , e, collapse = "\n")
    })
    
    
    #return
    p
}
