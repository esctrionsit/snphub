lp_main <- function(co, ro, lp_ntrans, ext = "0"){
	withProgress(message = 'Drawing', detail = "  Drawing SnpFreq...", value = 5, {
		text_lp_currpara <<- ""

	    lp_chr <- pub_check_chr_value(ro)

	    if(is.na(as.numeric(ext)) || as.numeric(ext) < 0) { return(lp_error_message(1)) }
	    if(is.na(as.numeric(lp_ntrans)) || as.numeric(lp_ntrans) < 0) { return(lp_error_message(5)) }

	    lp_beg <- pub_check_pos_begin(ro, ext)
	    lp_end <- pub_check_pos_end(ro, ext)

	    if(sum(is.na(lp_chr), is.na(lp_beg), is.na(lp_end)) > 0 || length(lp_chr) == 0) { return(lp_error_message(2)) }
	    if(lp_range_is_too_long(lp_beg, lp_end)) { return(lp_error_message(3)) }

	    lp_sam <- pub_check_gro_sample_name(c(co))
	    if(sum(is.na(lp_sam)) > 0 || length(lp_sam) == 0) { return(lp_error_message(4)) }

        lp_gro <- pub_sub_group(co)
        if(sum(is.na(lp_gro)) > 0 || nrow(lp_gro) == 0) { return(lp_error_message(6)) }

	    lp_chr <- lp_chr[1]
	    lp_beg <- lp_beg[1]
	    lp_end <- lp_end[1]

	    #fet_code <- lp_fetch_data(lp_chr, lp_beg, lp_end, lp_sam)

	    #if(fet_code != 0) { return(lp_error_message(fet_code + 100)) }

	    tran_code <- lp_trans_data(lp_chr, lp_beg, lp_end, lp_gro)
	    if(tran_code != 0) { return(lp_error_message(tran_code + 100)) }

	    lp_ntrans <- as.numeric(lp_ntrans)
	    p <- lp_draw_plot(lp_chr, lp_ntrans, lp_beg, lp_end, co, ro)

	    text_lp_currpara <<- paste("Parameter: ", co, " ; ", ro, sep="")

	    p
	})
}

lp_trans_data <- function(lp_chr, lp_beg, lp_end, lp_gro){
	groups <- pub_unique_sample(as.character(lp_gro[,2]))
	res <- data.frame()
	count <- 0
	for(group in groups){
		#fetch data
		sam <- lp_gro[which(lp_gro$Group == group),]$Sample
		fet_code <- lp_fetch_data(lp_chr, lp_beg, lp_end, sam)
		if(fet_code != 0) { return(fet_code) }

		inf <- fra_lp_orivcf$ANN
    	ac <- fra_lp_orivcf$AC

		#init result data frame
		if(count == 0){
			leng <- nrow(fra_lp_orivcf)
    		res <- data.frame(fra_lp_orivcf$CHROM, fra_lp_orivcf$POS, rep(0, time=leng))
    		for(i in 1:length(inf)){
				ANN <- inf[i]
				if(grepl("missense_variant", ANN)){
					res[i,3] <- 1
				}else if(grepl("synonymous_variant", ANN)){
					res[i,3] <- 2
				}else if(grepl("frameshift_variant", ANN)){
					res[i,3] <- 3
				}else if(grepl("stop_gained", ANN) || grepl("stop_lost", ANN)){
					res[i,3] <- 4
				}else if(grepl("splice_region_variant", ANN)){
					res[i,3] <- 5
				}else{
					res[i,3] <- 6
				}
				names(res) <- c("CHROM", "POS", "TYPE")
			}
		}

		#adding group freq
		tmp_col <- c()
		tmp_an <- as.numeric(fra_lp_orivcf$AN)
		for(i in 1:length(inf)){
			if(grepl(",", ac[i])){
				lis <- strsplit(ac[i], split = ",")[[1]]
				tmp_col <- c(tmp_col, as.numeric(lis[1])/tmp_an[i])
			}else{
				tmp_col <- c(tmp_col, as.numeric(ac[i])/tmp_an[i])
			}
		}

		#combine to data frame
		tmp_name <- c(names(res), group)
		res <- cbind(res, tmp_col)
		names(res) <- tmp_name

		#
		count <- count+1
	}

    fra_lp_drawdata <<- res
    #return
    0
}

lp_draw_plot <- function(chr, ntrans, beg, end, co, ro){
    # Read the input file
	mainDF <- fra_lp_drawdata
	n_group <- ncol(mainDF)-3
	group_names <- colnames(mainDF)[-c(1,2,3)]
	#
	
	LeftMost = as.numeric(beg)
	RightMost = as.numeric(end)
	#
	# Read the annotation file
	# tabix iwgsc_refseqv1.0_HighConf_2017Mar13.gff3.gz chr7A:1-100000
	tmpannoDF <- system2(path_tabix, args = c(path_gff3, paste0(chr, ":", LeftMost, "-", RightMost)), stdout = T)
	annoDF <- data.frame(do.call(rbind, strsplit(tmpannoDF, "\t", fixed=TRUE)), stringsAsFactors=F)
	if(ncol(annoDF)>4){
	    annoDF[,4] <- as.numeric(annoDF[,4])
	    annoDF[,5] <- as.numeric(annoDF[,5])
	    colnames(annoDF) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
	    annoDF$ID <- lp_getAttributeField(annoDF$attributes, "ID")
	    annoDF$Par <- lp_getAttributeField(annoDF$attributes, "Parent")
	    N_anno <- sum(annoDF$feature == "mRNA", na.rm = TRUE)
	    if(ntrans < N_anno) {N_anno <- ntrans}
	}else{
	    N_anno <- 0
	}
	#
	# ========================
	# plot region definition
	# plot1: main plot; plot2: header and legend
	layout(mat = matrix(c(2,1,3),ncol = 1), heights = c(1,5,1))
	# xlim = c(LeftMost, RightMost)
	# ylim = c(1,10)
	#
	color = c("#00D9FFFF", "#00FF19FF", "#FF0000FF", "#FFB300FF","#0026FFFF", "gray")
	names(color) <- c("missense_variant", "synonymous_variant", "frameshift_variant", "stop_gained/stop_lost", "splice_region_variant", "others")
	#
	# ========================
	# plot1: main plot
	#
	par(mar=c(0,6,0,3))
	xlength = RightMost - LeftMost
	# Set the xlim and ylim
	plot(NULL, NULL, 
	    xlim = c(LeftMost, RightMost), 
	    ylim = c(1 - N_anno*1.5 - 1, 6 + (n_group-1)*6),
	    xlab = "", ylab = "",
	    yaxt = "n", xaxt = "n", frame.plot = FALSE,
    	cex.axis = 1.5, cex = 1.5)
	# -------------------------
	# plot axis
	axis(1, pos = 1, lwd = 4, cex.axis = 2, xaxs = "i")
	axis(2, 1:5, c("0%", "25%", "50%", "75%", "100%"), las=2, cex.axis = 2, lwd = 4)
	# Draw the DNA strands
	#abline( h = 1, lwd = 2 )
	# plot axis
	# xaxs does not work as expected
	for (j in seq_len(n_group)){
		abline( h = 1 + (j-1)*6, lwd = 4 )
		axis(1, pos = 1, lwd = 4, cex.axis = 2, xaxs = "i")
		axis(2, at = c(1 + (j-1)*6, 5 + (j-1)*6), labels = c("0","1"), las=2, cex.axis = 2, lwd = 4)
		text(x = LeftMost, y = j*6, group_names[j], pos = 1, cex=2.5)
		# Draw the loll
		for( i in seq_len(nrow(mainDF))){
			# lines( c(pos[i], pos[i]), c(0, mC[i,k])+k*yscale , lwd=5)
			segments( mainDF[i,2], 1 + (j-1)*6, mainDF[i,2], 1 + mainDF[i,3+j]*4 + (j-1)*6, lwd = 3, col = color[mainDF[i,3]] )
			points( mainDF[i,2], 1 + mainDF[i,3+j]*4 + (j-1)*6, cex = 2, pch = 16, col = color[mainDF[i,3]] )
		}
	}
	#
	# ==========================
	# Draw the gene annotation
	if(N_anno > 0){
	    Ttrans <- annoDF[annoDF$feature == "mRNA", "ID"]
	    for (i in c(1:N_anno)){
	        trans <- Ttrans[i]
	        annoDFtmp <- rbind(annoDF[annoDF$feature == "mRNA" & annoDF$ID == trans,], 
	                                             annoDF[annoDF$Par == trans, ])
	        annoDFtmp <- annoDFtmp[!is.na(annoDFtmp$seqname),]
	        
	        y_center = 0.5 - N_anno*1.5 + i*1.5 - 1.5
	        # gene name on upmost
	        for (j in seq_len(nrow(annoDFtmp)) ) {
	            annotype = annoDFtmp[j,3]
	            Left = annoDFtmp[j,4]
	            Right = annoDFtmp[j,5]
	            x_center = Left + (Right-Left)/2
	            
	            if (annotype == "exon"){
	                rect(Left, y_center-0.2, Right, y_center+0.2, col = "white")
	            } else if (annotype == "CDS"){
	                rect(Left, y_center-0.2, Right, y_center+0.2, col = "black")
	            } else if (annotype == "mRNA"){
	                strand = annoDFtmp[j,7]
	                StrandSign=(c(1,-1)[c("+","-")==strand])
	                x0 = c(Left, Right)[c("+","-")==strand]
	                x1 = x0 - (RightMost-LeftMost)/30*StrandSign
	                x2 = c(Right, Left)[c("+","-")==strand]
	                segments(Left, y_center, Right, y_center, col = "black")
			        arrows(x0, y_center, x1, y_center, length = 0.08 )
			        text(x0, y_center+0.4, labels = annoDFtmp[j,10], cex = 1, pos = 4)
	            }
	        }
	    }
	}
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
	text(0.5, 1, cex = 3, "SnpFreq", pos = 1)
	# Draw the color legends
	legend("bottomleft", fill = color, legend = names(color), border = F, box.col = "white", horiz=T, cex = 2, x.intersp = 0.5)
	
	#
	# ========================
	# plot3: tag and paras
	#
	t = read.table(pipe("date +%Y%m%d%H%M%S"))
    d = as.character(t[1,1])
    filenametag <- paste("SnpHub_SnpFreq_", d, sep="")
    parameter <- paste("Parameter: ", co, "; ", ro,";", sep="")
    #
	par( mar = c(3,3,0,3) )
	plot(NULL, NULL, type="n", bty="n", xaxt="n", yaxt="n", 
	     xlab="", ylab="", 
	     xlim=c(0, 1), ylim=c(0, 1), xaxs="i", yaxs="i")
	text(x = 1, y = 0.5, paste0(filenametag, "\n", parameter), pos = 2, cex = 2)

}

lp_getAttributeField <- function (x, field, attrsep = ";") {
	s = strsplit(x, split = attrsep, fixed = TRUE)
	sapply(s, function(atts) {
		a = strsplit(atts, split = "=", fixed = TRUE)
    	m = match(field, sapply(a, "[", 1))
    	if (!is.na(m)) {
    		rv = a[[m]][2]
    	}else {
    		rv = as.character(NA)
    	}
    	return(rv)
  	})
}

lp_range_is_too_long <- function(lp_beg, lp_end){ F }

lp_error_message <- function(code){
    code
}
