hp_main <- function(co, ro, ext, cluster, flip, maf="0") {
    withProgress(message = 'Drawing', detail = "  Drawing Haplotype Heatmap plot...", value = 5, {
        text_hp_currpara <<- ""

        hp_chr <- pub_check_chr_value(ro)

        if(is.na(as.numeric(ext)) || as.numeric(ext) < 0) { return(hp_error_message(1)) }
        if(is.na(as.numeric(maf))) { return(hp_error_message(7)) }

        hp_beg <- pub_check_pos_begin(ro, ext)
        hp_end <- pub_check_pos_end(ro, ext)

        if(sum(is.na(hp_chr), is.na(hp_beg), is.na(hp_end)) > 0 || length(hp_chr) == 0) { return(hp_error_message(2)) }
        if(hp_range_is_too_long(hp_beg, hp_end)) { return(hp_error_message(3)) }
        
        hp_sam <- hp_check_sample_name(co)
        hp_gro <- pub_sub_group(co)

        if(sum(is.na(hp_sam)) > 0 || length(hp_sam) == 0) { return(hp_error_message(4)) }
        if(sum(is.na(hp_gro)) > 0 || nrow(hp_gro) == 0) { return(hp_error_message(6)) }

        unique_sam <- pub_unique_sample(hp_sam)
        
        hp_fetch_data(hp_chr, hp_beg, hp_end, unique_sam, maf)

        if(nrow(fra_hp_orivcf) == 0) { return(hp_error_message(5)) }

        fra_hp_orivcf <<- subset(fra_hp_orivcf, select=c("CHROM", "POS", hp_sam))

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

hp_check_sample_name <- function(co) {
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

hp_error_message <- function(code) {
    code
}

hp_trans_data <- function(hp_sam, hp_gro) {
	new_row_name <- paste(fra_hp_orivcf$CHROM, fra_hp_orivcf$POS, sep=":")

	fra_hp_orivcf <<- as.matrix(fra_hp_orivcf)

	fra_hp_orivcf[which(fra_hp_orivcf %in% lis_hp_code_missing)] <<- "-1"
	fra_hp_orivcf[which(fra_hp_orivcf %in% lis_hp_code_nosnp)] <<- "0"
	fra_hp_orivcf[which(fra_hp_orivcf %in% lis_hp_code_hete)] <<- "0.5"
	fra_hp_orivcf[which(fra_hp_orivcf %in% lis_hp_code_homo)] <<- "1"

	fra_hp_orivcf <<- as.data.frame(t(fra_hp_orivcf[,-(1:2)]))
	colnames(fra_hp_orivcf) <<- new_row_name
    
	fra_hp_orivcf$Sample <<- hp_sam
    fra_hp_orivcf$Group <<- hp_gro$Group
}

hp_draw_plot <- function(co, ro, ext, cluster, flip) {
	# add info
	sample_present <- fra_glo_metadata[fra_glo_metadata$Label %in% fra_hp_orivcf$Sample,]
	LIST <- sample_present$Label
	names(LIST) <- sample_present$Name
    new_sam <- c()
    for(j in fra_hp_orivcf$Sample){
        tmp3 <- as.character(fra_glo_metadata[which(fra_glo_metadata$Accession == j),]$Name)
        if(length(tmp3) != 1) { tmp3 <- NA }
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
	levels(df$Type)[levels(df$Type)=="-1"] <- "Missing"
	levels(df$Type)[levels(df$Type)=="0"] <- "None"
	levels(df$Type)[levels(df$Type)=="0.5"] <- "Heter"
	levels(df$Type)[levels(df$Type)=="1"] <- "Homo"

    #get lab
    t = read.table(pipe("date +%Y%m%d%H%M%S"))
    d = as.character(t[1,1])
    filenametag <- paste("SnpHub_Heatmap_", d, sep="")
    parameter <- paste("Parameter: ", co, "; ", ro,"; flanking ", ext, ";", sep="")
    
	# plot
    if(flip == "Yes"){
        p <- ggplot(df,aes(Mutation,Sample,fill=Type))+
            geom_tile() +
            scale_fill_manual(values = c("gray95", "grey", "deepskyblue", "blue")) +
            #scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30)) +
            scale_x_discrete(limits = rev(levels(df$Mutation))) +
            theme_minimal()+
            labs(title = "Haplotype heatmap", caption = paste0(filenametag, "\n", parameter)) +
            facet_grid(rows = vars(Group), scales = "free", space = "free") +
            theme(legend.position = "top",
                legend.direction = "horizontal",
                legend.title = element_blank(),
                legend.justification = "right",
                plot.title = element_text(size = 17,hjust=0.5),
                axis.text.x = element_text(angle = 70, hjust = 1),
                axis.title = element_text(size=25),
                strip.text = element_text(size=25),
                aspect.ratio = 1) +
            NULL
    }else{
        p <- ggplot(df,aes(Sample,Mutation,fill=Type))+
            geom_tile() +
            scale_fill_manual(values = c("gray95", "grey", "deepskyblue", "blue")) +
            #scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30)) +
            scale_y_discrete(limits = levels(df$Mutation)) +
            theme_minimal()+
            labs(title = "Haplotype heatmap", caption = paste0(filenametag, "\n", parameter)) +
            facet_grid(cols = vars(Group), scales = "free", space = "free") +
            theme(legend.position = "top",
                  legend.direction = "horizontal",
                  legend.title = element_blank(),
                  legend.justification = "right",
                  plot.title = element_text(size = 17,hjust=0.5),
                  axis.text.x = element_text(angle = 70, hjust = 1),
                  axis.title = element_text(size=25),
                  strip.text = element_text(size=25),
                  aspect.ratio = 1) +
            NULL
    }
    #return
    p
}

