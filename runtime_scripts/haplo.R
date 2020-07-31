hp_main <- function(co, ro, ext="0", cluster="Yes", flip="Yes", maf="0", maf_igm="Yes") {
    withProgress(message = 'Drawing', detail = "  Drawing Haplotype Heatmap plot...", value = 5, {
        library(ggplot2)

        co <- paste(strsplit(co, split = " ")[[1]], collapse="")
        text_hp_currpara <<- ""

        hp_chr <- pub_check_chr_value(ro)

        if(is.na(as.numeric(ext)) || as.numeric(ext) < 0) { return(hp_error_message(1)) }
        if(is.na(as.numeric(maf))) { return(hp_error_message(7)) }

        hp_beg <- pub_check_pos_begin(ro, ext)
        hp_end <- pub_check_pos_end(ro, ext)

        if(sum(is.na(hp_chr), is.na(hp_beg), is.na(hp_end)) > 0 || length(hp_chr) == 0) { return(hp_error_message(2)) }
        if(hp_range_is_too_long(hp_beg, hp_end)) { return(hp_error_message(3)) }
        
        hp_sam <- pub_check_gro_sample_name(co)
        hp_gro <- pub_sub_group(co)

        if(sum(is.na(hp_sam)) > 0) { return(hp_error_message(4)) }
        if(sum(is.na(hp_gro)) > 0 || nrow(hp_gro) == 0) { return(hp_error_message(6)) }

        hp_sam <- as.character(hp_gro[,1])

        unique_sam <- pub_unique_sample(hp_sam)
        
        hp_fetch_data(hp_chr, hp_beg, hp_end, unique_sam, maf, maf_igm)

        if(nrow(fra_hp_orivcf) == 0) { return(hp_error_message(5)) }

        fra_hp_orivcf <<- subset(fra_hp_orivcf, select=c("CHROM", "POS", "POS", "ANN", hp_sam))

        tmp_df <- data.frame(c("GENE","MUTYPE", "_"),c("TYPE","TYPE", "_"))
        names(tmp_df) <- c("Sample", "Group")
        hp_gro <- rbind(tmp_df, hp_gro)

        hp_trans_data(hp_sam, hp_gro)

        text_hp_currpara <<- paste("Parameter: ", ro, " ; flask length ", ext, sep="")

        p <- hp_draw_plot(co, ro, ext, cluster, flip)  
        #return
        p
    })
}

hp_range_is_too_long <- function(hp_beg, hp_end) {
    leng <- 0;
    for(i in 1:length(hp_beg)){
        leng <- leng + (as.numeric(hp_end[i]) - as.numeric(hp_beg[i]))
    }
    if(leng > 5500000){
        T
    }else{
        F
    }
}

hp_error_message <- function(code) {
    code
}

hp_trans_data <- function(hp_sam, hp_gro) {
	new_row_name <- paste(fra_hp_orivcf$CHROM, fra_hp_orivcf$POS, sep=":")

    tmp_name <- names(fra_hp_orivcf)
    tmp_name[3] <- "GENE"
    tmp_name[4] <- "MUTYPE"
    names(fra_hp_orivcf) <<- tmp_name

    assist_line <- rep("-2", nrow(fra_hp_orivcf))
    assist_val <- 1
    for(i in 1:nrow(fra_hp_orivcf)){
        if(assist_val > 9){ break }
        if(i == 1){ next }
        assist_line[i] <- paste("-2.", as.character(assist_val), sep="")
        assist_val <- assist_val+1
    }
    fra_hp_orivcf <<- cbind(fra_hp_orivcf[,1:4], assist_line, fra_hp_orivcf[,5:ncol(fra_hp_orivcf)])
    tmp_name <- c(tmp_name[1:4], "_", tmp_name[5:length(tmp_name)])
    names(fra_hp_orivcf) <<- tmp_name

	fra_hp_orivcf <<- as.matrix(fra_hp_orivcf)

	fra_hp_orivcf[which(fra_hp_orivcf %in% lis_hp_code_missing)] <<- "-1"
	fra_hp_orivcf[which(fra_hp_orivcf %in% lis_hp_code_nosnp)] <<- "0"
	fra_hp_orivcf[which(fra_hp_orivcf %in% lis_hp_code_hete)] <<- "0.5"
	fra_hp_orivcf[which(fra_hp_orivcf %in% lis_hp_code_homo)] <<- "1"

    tmp_beg <- 1
    tmp_ind <- subset(fra_glo_geneindex, CHROM == fra_hp_orivcf[1,1])
    for(i in 1:nrow(fra_hp_orivcf)){
        for(j in tmp_beg:nrow(tmp_ind)){
            if(as.numeric(tmp_ind[j,3]) < as.numeric(fra_hp_orivcf[i,3])){
                tmp_beg <- tmp_beg + 1
                next
            }
            if(as.numeric(tmp_ind[j,2]) < as.numeric(fra_hp_orivcf[i,3]) && as.numeric(tmp_ind[j,3]) > as.numeric(fra_hp_orivcf[i,3])){
                fra_hp_orivcf[i,3] <<- "3.0"
                break
            }
            if(as.numeric(tmp_ind[j,2]) > as.numeric(fra_hp_orivcf[i,3])){
                fra_hp_orivcf[i,3] <<- "3.1"
                break
            }
        }
    }

    for(i in 1:nrow(fra_hp_orivcf)){
        ANN <- fra_hp_orivcf[i, 4]
        if(grepl("missense_variant", ANN)){
            fra_hp_orivcf[i,4] <<- "2.1"
        }else if(grepl("synonymous_variant", ANN)){
            fra_hp_orivcf[i,4] <<- "2.2"
        }else if(grepl("frameshift_variant", ANN)){
            fra_hp_orivcf[i,4] <<- "2.3"
        }else if(grepl("stop_gained", ANN) || grepl("stop_lost", ANN)){
            fra_hp_orivcf[i,4] <<- "2.4"
        }else if(grepl("splice_region_variant", ANN)){
            fra_hp_orivcf[i,4] <<- "2.5"
        }else{
            fra_hp_orivcf[i,4] <<- "2.0"
        }
    }

    if(nrow(fra_hp_orivcf) != 1){
        fra_hp_orivcf <<- as.data.frame(t(fra_hp_orivcf[,-(1:2)]))
    }else{
        fra_hp_orivcf <<- as.data.frame(t(t(fra_hp_orivcf[,-(1:2)])))
    }
	colnames(fra_hp_orivcf) <<- new_row_name

	fra_hp_orivcf$Sample <<- c("GENE", "MUTYPE", "_", hp_sam)
    fra_hp_orivcf$Group <<- hp_gro$Group
}

hp_draw_plot <- function(co, ro, ext, cluster, flip) {
	# add info
	#sample_present <- fra_glo_metadata[fra_glo_metadata$Label %in% fra_hp_orivcf$Sample,]
    #add_sam <- data.frame(c("GENE","MUTYPE"), c("GENE","MUTYPE"), c("GENE","MUTYPE"))
    #names(add_sam) <- c("Accession", "Label", "Name")
    #sample_present <- rbind(sample_present, add_sam)
	#LIST <- sample_present$Label
	#names(LIST) <- sample_present$Name
    new_sam <- c()
    for(j in fra_hp_orivcf$Sample){
        if(as.character(j)=="GENE" || as.character(j)=="MUTYPE" || as.character(j)=="_"){
            new_sam <- c(new_sam, as.character(j))
            next
        }
        tmp3 <- as.character(fra_glo_metadata[which(fra_glo_metadata$Accession == j),]$Name)
        if(length(tmp3) != 1) { tmp3 <- NA }
        if(tmp3 %in% new_sam){
            assist_index <- 1
            while(T){
                if(paste(tmp3, ".", as.character(assist_index), sep="") %in% new_sam){
                    assist_index <- assist_index + 1
                    next
                }
                tmp3 <- paste(tmp3, ".", as.character(assist_index), sep="")
                break
            }
        }
        new_sam <- c(new_sam, tmp3)
    }
    fra_hp_orivcf$Sample <<- new_sam
	#fra_hp_orivcf$Sample <<- LIST[fra_hp_orivcf$Sample]

    if(flip == "Yes"){
        fra_hp_orivcf$Sample <<- factor(fra_hp_orivcf$Sample, levels = rev(fra_hp_orivcf$Sample))
    } else{
        fra_hp_orivcf$Sample <<- factor(fra_hp_orivcf$Sample, levels = fra_hp_orivcf$Sample)
    }

	if(cluster == "Yes"){
        # re-order
        tryCatch(
            error = function(cnd) {
                fra_hp_orivcf$Sample <<- factor(fra_hp_orivcf$Sample)
            },
            {
                suppressWarnings(row.order <- hclust(dist(fra_hp_orivcf[,-1]))$order)
                fra_hp_orivcf$Sample <<- factor(fra_hp_orivcf$Sample, levels = fra_hp_orivcf$Sample[row.order])
            }
        )
    }

	# melt
	df <- reshape2::melt(fra_hp_orivcf, id.vars = c("Sample", "Group"))
	colnames(df) <- c("Sample", "Group", "Mutation", "Type")
    df$Group <- factor(df$Group, levels = unique(df$Group))

	# turn discrete
	df$Type <- as.factor(df$Type)
    color_lis <- c()
    count_L <- 0
    count_M <- 0
    count_R <- 0
    color_lis_L <- c()
    color_lis_M <- c()
    color_lis_R <- c()
    if("-1" %in% df$Type){
        levels(df$Type)[levels(df$Type)=="-1"] <- "Missing"
        count_L <- count_L+1
        color_lis_L <- c(color_lis_L, "Missing"="gray95")
    }
	if("0" %in% df$Type){
	   levels(df$Type)[levels(df$Type)=="0"] <- "None"
       count_L <- count_L+1
       color_lis_L <- c(color_lis_L, "None"="grey")
    }
    if("0.5" %in% df$Type){
	   levels(df$Type)[levels(df$Type)=="0.5"] <- "Heter"
       count_L <- count_L+1
       color_lis_L <- c(color_lis_L, "Heter"="deepskyblue")
    }
    if("1" %in% df$Type){
	   levels(df$Type)[levels(df$Type)=="1"] <- "Homo"
       count_L <- count_L+1
       color_lis_L <- c(color_lis_L, "Homo"="blue")
    }
    if("2.0" %in% df$Type){
        levels(df$Type)[levels(df$Type)=="2.0"] <- "others"
        count_R <- count_R+1
        color_lis_R <- c(color_lis_R, "others"="#484848")
    }
    if("2.1" %in% df$Type){
        levels(df$Type)[levels(df$Type)=="2.1"] <- "missense_variant"
        count_R <- count_R+1
        color_lis_R <- c(color_lis_R, "missense_variant"="#00D9FFFF")
    }
    if("2.2" %in% df$Type){
        levels(df$Type)[levels(df$Type)=="2.2"] <- "synonymous_variant"
        count_R <- count_R+1
        color_lis_R <- c(color_lis_R, "synonymous_variant"="#00FF19FF")
    }
    if("2.3" %in% df$Type){
        levels(df$Type)[levels(df$Type)=="2.3"] <- "frameshift_variant"
        count_R <- count_R+1
        color_lis_R <- c(color_lis_R, "frameshift_variant"="#FF0000FF")
    }
    if("2.4" %in% df$Type){
        levels(df$Type)[levels(df$Type)=="2.4"] <- "stop_gained/stop_lost"
        count_R <- count_R+1
        color_lis_R <- c(color_lis_R, "stop_gained/stop_lost"="#FFB300FF")
    }
    if("2.5" %in% df$Type){
        levels(df$Type)[levels(df$Type)=="2.5"] <- "splice_region_variant"
        count_R <- count_R+1
        color_lis_R <- c(color_lis_R, "splice_region_variant"="#0026FFFF")
    }
    if("3.0" %in% df$Type){
        levels(df$Type)[levels(df$Type)=="3.0"] <- "gene_body"
        count_M <- count_M+1
        color_lis_M <- c(color_lis_M, "gene_body"="#dfbb27")
    }
    if("3.1" %in% df$Type){
        levels(df$Type)[levels(df$Type)=="3.1"] <- "intergenic"
        count_M <- count_M+1
        color_lis_M <- c(color_lis_M, "intergenic"="black")
    }

    color_lis <- c()
    assist_factor <- c(" ", "  ", "   ", "    ", "     ", "      ", "       ", "        ", "         ", "          ")
    if(count_L == count_M && count_L == count_R){
        count_assist <- 0
        if("-2" %in% df$Type){
            if(count_assist < 3){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2"] <- assist_factor[count_assist]
        }
        if("-2.1" %in% df$Type){
            if(count_assist < 3){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.1"] <- assist_factor[count_assist]
        }
        if("-2.2" %in% df$Type){
            if(count_assist < 3){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.2"] <- assist_factor[count_assist]
        }
        if("-2.3" %in% df$Type){
            if(count_assist < 3){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.3"] <- assist_factor[count_assist]
        }
        if("-2.4" %in% df$Type){
            if(count_assist < 3){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.4"] <- assist_factor[count_assist]
        }
        if("-2.5" %in% df$Type){
            if(count_assist < 3){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.5"] <- assist_factor[count_assist]
        }
        if("-2.6" %in% df$Type){
            if(count_assist < 3){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.6"] <- assist_factor[count_assist]
        }
        if("-2.7" %in% df$Type){
            if(count_assist < 3){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.7"] <- assist_factor[count_assist]
        }
        if("-2.8" %in% df$Type){
            if(count_assist < 3){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.8"] <- assist_factor[count_assist]
        }
        if("-2.9" %in% df$Type){
            if(count_assist < 3){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.9"] <- assist_factor[count_assist]
        }
        if(count_assist>=1){
            color_lis <- c(color_lis, "white")
            tmp_lis_name <- names(color_lis)
            tmp_lis_name[length(tmp_lis_name)] <- assist_factor[1]
            names(color_lis) <- tmp_lis_name
        }
        color_lis <- c(color_lis, color_lis_L)
        if(count_assist>=2){
            color_lis <- c(color_lis, "white")
            tmp_lis_name <- names(color_lis)
            tmp_lis_name[length(tmp_lis_name)] <- assist_factor[2]
            names(color_lis) <- tmp_lis_name
        }
        color_lis <- c(color_lis, color_lis_M)
        if(count_assist>=3){
            color_lis <- c(color_lis, "white")
            tmp_lis_name <- names(color_lis)
            tmp_lis_name[length(tmp_lis_name)] <- assist_factor[3]
            names(color_lis) <- tmp_lis_name
        }
        color_lis <- c(color_lis, color_lis_R)
    }else if(count_L >= count_M && count_L >= count_R){
        count_1 <- count_L-count_M
        count_2 <- count_L-count_R
        count_assist <- 0
        if("-2" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2"] <- assist_factor[count_assist]
        }
        if("-2.1" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.1"] <- assist_factor[count_assist]
        }
        if("-2.2" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.2"] <- assist_factor[count_assist]
        }
        if("-2.3" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.3"] <- assist_factor[count_assist]
        }
        if("-2.4" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.4"] <- assist_factor[count_assist]
        }
        if("-2.5" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.5"] <- assist_factor[count_assist]
        }
        if("-2.6" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.6"] <- assist_factor[count_assist]
        }
        if("-2.7" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.7"] <- assist_factor[count_assist]
        }
        if("-2.8" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.8"] <- assist_factor[count_assist]
        }
        if("-2.9" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.9"] <- assist_factor[count_assist]
        }
        print(c(count_1, count_2))
        print(c(count_L, count_M, count_R, count_assist))
        print(color_lis)
        color_lis <- c(color_lis, color_lis_L)
        print(color_lis)
        if(count_assist > count_1){
            if(count_1>0){
                    for(i in 1:count_1){
                    color_lis <- c(color_lis, "white")
                    tmp_lis_name <- names(color_lis)
                    tmp_lis_name[length(tmp_lis_name)] <- assist_factor[i]
                    names(color_lis) <- tmp_lis_name
                }
            }
        }else{
            for(i in 1:count_assist){
                color_lis <- c(color_lis, "white")
                tmp_lis_name <- names(color_lis)
                tmp_lis_name[length(tmp_lis_name)] <- assist_factor[i]
                names(color_lis) <- tmp_lis_name
            }
        }
        print(color_lis)
        color_lis <- c(color_lis, color_lis_M)
        print(color_lis)
        print(c(count_1, count_assist))
        if(count_assist > count_1){
            for(i in (count_1+1):count_assist){
                print(i)
                color_lis <- c(color_lis, "white")
                tmp_lis_name <- names(color_lis)
                tmp_lis_name[length(tmp_lis_name)] <- assist_factor[i]
                names(color_lis) <- tmp_lis_name
            }
        }
        print(color_lis)
        color_lis <- c(color_lis, color_lis_R)
        print(color_lis)
    }else if(count_M >= count_L && count_M >= count_R){
        count_1 <- count_M-count_L
        count_2 <- count_M-count_R
        count_assist <- 0
        if("-2" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2"] <- assist_factor[count_assist]
        }
        if("-2.1" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.1"] <- assist_factor[count_assist]
        }
        if("-2.2" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.2"] <- assist_factor[count_assist]
        }
        if("-2.3" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.3"] <- assist_factor[count_assist]
        }
        if("-2.4" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.4"] <- assist_factor[count_assist]
        }
        if("-2.5" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.5"] <- assist_factor[count_assist]
        }
        if("-2.6" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.6"] <- assist_factor[count_assist]
        }
        if("-2.7" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.7"] <- assist_factor[count_assist]
        }
        if("-2.8" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.8"] <- assist_factor[count_assist]
        }
        if("-2.9" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.9"] <- assist_factor[count_assist]
        }

        if(count_assist > count_1){
            if(count_1 > 0){
                    for(i in 1:count_1){
                    color_lis <- c(color_lis, "white")
                    tmp_lis_name <- names(color_lis)
                    tmp_lis_name[length(tmp_lis_name)] <- assist_factor[i]
                    names(color_lis) <- tmp_lis_name
                }
            }
        }else{
            for(i in 1:count_assist){
                color_lis <- c(color_lis, "white")
                tmp_lis_name <- names(color_lis)
                tmp_lis_name[length(tmp_lis_name)] <- assist_factor[i]
                names(color_lis) <- tmp_lis_name
            }
        }
        color_lis <- c(color_lis, color_lis_L)
        color_lis <- c(color_lis, color_lis_M)
        if(count_assist > count_1){
            for(i in (count_1+1):count_assist){
                color_lis <- c(color_lis, "white")
                tmp_lis_name <- names(color_lis)
                tmp_lis_name[length(tmp_lis_name)] <- assist_factor[i]
                names(color_lis) <- tmp_lis_name
            }
        }
        color_lis <- c(color_lis, color_lis_R)
    }else{
        count_1 <- count_R-count_L
        count_2 <- count_R-count_M
        count_assist <- 0
        if("-2" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2"] <- assist_factor[count_assist]
        }
        if("-2.1" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.1"] <- assist_factor[count_assist]
        }
        if("-2.2" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.2"] <- assist_factor[count_assist]
        }
        if("-2.3" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.3"] <- assist_factor[count_assist]
        }
        if("-2.4" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.4"] <- assist_factor[count_assist]
        }
        if("-2.5" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.5"] <- assist_factor[count_assist]
        }
        if("-2.6" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.6"] <- assist_factor[count_assist]
        }
        if("-2.7" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.7"] <- assist_factor[count_assist]
        }
        if("-2.8" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.8"] <- assist_factor[count_assist]
        }
        if("-2.9" %in% df$Type){
            if(count_assist < count_1+count_2){ count_assist <- count_assist + 1 }
            levels(df$Type)[levels(df$Type)=="-2.9"] <- assist_factor[count_assist]
        }

        if(count_assist > count_1){
            if(count_1 > 0){
                    for(i in 1:count_1){
                    color_lis <- c(color_lis, "white")
                    tmp_lis_name <- names(color_lis)
                    tmp_lis_name[length(tmp_lis_name)] <- assist_factor[i]
                    names(color_lis) <- tmp_lis_name
                }
            }
        }else{
            for(i in 1:count_assist){
                color_lis <- c(color_lis, "white")
                tmp_lis_name <- names(color_lis)
                tmp_lis_name[length(tmp_lis_name)] <- assist_factor[i]
                names(color_lis) <- tmp_lis_name
            }
        }
        color_lis <- c(color_lis, color_lis_L)
        if(count_assist > count_1){
            for(i in (count_1+1):count_assist){
                color_lis <- c(color_lis, "white")
                tmp_lis_name <- names(color_lis)
                tmp_lis_name[length(tmp_lis_name)] <- assist_factor[i]
                names(color_lis) <- tmp_lis_name
            }
        }
        color_lis <- c(color_lis, color_lis_M)
        color_lis <- c(color_lis, color_lis_R)
    }
    df$Type <- ordered(as.character(df$Type), levels = names(color_lis))

    print(df)

    #get lab
    t = read.table(pipe("date +%Y%m%d%H%M%S"))
    d = as.character(t[1,1])
    filenametag <- paste("SnpHub_Heatmap_", d, sep="")
    parameter <- paste("Parameter: ", co, "; ", ro,"; flanking ", ext, ";", sep="")
    print(color_lis)
	# plot
    if(flip == "Yes"){
        p <- ggplot(df,aes(Mutation,Sample,fill=Type))+
            geom_tile() +
            scale_fill_manual(values = color_lis, drop=F) +
            #scale_color_discrete(limits = names(color_lis))+
            guides(fill=guide_legend(ncol=3)) +
            #scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30)) +
            scale_x_discrete(limits = rev(levels(df$Mutation))) +
            theme_minimal()+
            labs(title = "Genotype heatmap", caption = paste0(filenametag, "\n", parameter)) +
            facet_grid(rows = vars(Group), scales = "free", space = "free") +
            theme(legend.position = "top",
                legend.direction = "horizontal",
                legend.title = element_blank(),
                legend.justification = "right",
                plot.title = element_text(size = 30,hjust=0.5),
                axis.text.x = element_text(size = 17, angle = 70, hjust = 1),
                axis.text.y = element_text(size=14),
                axis.title = element_text(size=25),
                strip.text = element_text(size=25),
                aspect.ratio = 1) +
            NULL
    }else{
        p <- ggplot(df,aes(Sample,Mutation,fill=Type))+
            geom_tile() +
            scale_fill_manual(values = color_lis, drop=F) +
            guides(fill=guide_legend(ncol=3)) +
            #scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30)) +
            scale_y_discrete(limits = levels(df$Mutation)) +
            theme_minimal()+
            labs(title = "Haplotype heatmap", caption = paste0(filenametag, "\n", parameter)) +
            facet_grid(cols = vars(Group), scales = "free", space = "free") +
            theme(legend.position = "top",
                  legend.direction = "horizontal",
                  legend.title = element_blank(),
                  legend.justification = "right",
                  plot.title = element_text(size = 30,hjust=0.5),
                  axis.text.x = element_text(size=14, angle = 70, hjust = 1),
                  axis.text.y = element_text(size=14),
                  axis.title = element_text(size=25),
                  strip.text = element_text(size=25),
                  aspect.ratio = 1) +
            NULL
    }
    #return
    p
}

