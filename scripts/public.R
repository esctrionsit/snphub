#packages
library(ggplot2)
library(ggtree)
library(rjson)
library(shiny)
library(pegas)
library(vcfR)
library(ape)
library(DT)
##########################################################################################################################
##########################################################################################################################
#init function
pub_read_sysinfo <- function(path_sysinfo){
    res <- c()

    con <- file(path_sysinfo, "r")
    line <- readLines(con,n=1)
    while( length(line) != 0 ) {
         res <- c(res, line)
         line <- readLines(con,n=1)
    }
    close(con)

    res_fra <- data.frame(res)
    names(res_fra) <- c("")

    res_fra
}
##########################################################################################################################
##########################################################################################################################
#global varities
fra_glo_geneindex <- read.table(path_geneindex, col.names = c("CHROM", "begin", "end", "gene_name"), header = F, as.is = T)
fra_glo_metadata_raw <- read.table(path_metadata, header = T, as.is = T, fill=T)
fra_glo_metadata <- fra_glo_metadata_raw[,1:3]
names(fra_glo_metadata) <-  c("Accession", "Label", "Name")
fra_glo_groupdata <- read.table(path_groupdata, col.names = c("Group", "Name"), header = F, as.is = T)
fra_glo_avalichrnames <- read.table(path_fa_index, col.names = c("CHROM", "Length", "V1", "V2", "V3"), as.is = T)
lis_glo_avalichrnames <- fra_glo_avalichrnames$CHROM
fra_glo_samloc <- read.table(path_sam_location, header = F, stringsAsFactors = F, col.names = c("Acession", "location", "lon", "lat"))
if(!is.na(path_sysinfo)){
    fra_glo_sysdata <- pub_read_sysinfo(path_sysinfo)
}else{
    fra_glo_sysdata <- data.frame()
}

if(is.na(path_UIsetting)){
    path_UIsetting <- "./scripts/empty.json"
    json_glo_UIsetting <- fromJSON(paste(readLines(path_UIsetting), collapse=""))
}else{
    json_glo_UIsetting <- fromJSON(paste(readLines(path_UIsetting), collapse=""))
}


fra_snp_orivcf <- data.frame(c("fra_snp_orivcf"))
fra_snp_finalres <- data.frame(c("fra_snp_finalres"))
text_snp_currpara <- ""
#newstring
fra_ns_finalres <- data.frame(c('fra_ns_finalres'))
#haplotype plot
fra_hp_orivcf <- data.frame(c("fra_hp_orivcf"))
text_hp_currpara <- ""
#hapNet
vcfR_hn_orivcf <- "vcfR file"
text_hn_currpara <- ""
#distance tree
vcfR_dt_orivcf <- "vcfR file"
bool_dt_warning <- F
text_dt_currpara <- ""
#hapmap
fra_hm_orivcf <- data.frame(c("fra_hm_orivcf"))
fra_hm_drawdata <- data.frame(c("fra_hm_drawdata"))
int_hm_warning <- 0
text_hm_NAwarning <- ""
text_hm_drawsite <- ""
#lollipop
fra_lp_orivcf <- data.frame(c("fra_lp_orivcf"))
fra_lp_drawdata <- data.frame(c("fra_lp_drawdata"))
text_lp_currpara <- ""
#general
fra_glo_samshow <- fra_glo_metadata_raw
names(fra_glo_samshow) -> gtmp
gtmp[1] <- "vcfID"
gtmp[2] <- "Accession name"
gtmp[3] <- "Display name"
names(fra_glo_samshow) <- gtmp
fra_glo_groshow <- fra_glo_groupdata
names(fra_glo_groshow) <- c("Group name", "Contained samples")
fra_glo_chrshow <- fra_glo_avalichrnames[,1:2]
names(fra_glo_chrshow) <- c("Chromosome", "Length in bp")
##########################################################################################################################
##########################################################################################################################
#In haplotype plot, mutation from raw vcf is changed into 4 types: missing, nosnp, homozygous and heterozygous.
lis_hp_code_missing <- c("./.")
lis_hp_code_nosnp <- c("0/0")
lis_hp_code_homo <- c("1/1", "2/2", "3/3", "4/4", "5/5", "6/6", "7/7", "8/8", "9/9")
lis_hp_code_hete <- c("0/1", "0/2", "0/3", "0/4", "0/5", "0/6", "0/7", "0/8", "0/9",
    "1/0", "1/2", "1/3", "1/4", "1/5", "1/6", "1/7", "1/8", "1/9",
    "2/0", "2/1", "2/3", "2/4", "2/5", "2/6", "2/7", "2/8", "2/9",
    "3/0", "3/1", "3/2", "3/4", "3/5", "3/6", "3/7", "3/8", "3/9",
    "4/0", "4/1", "4/2", "4/3", "4/5", "4/6", "4/7", "4/8", "4/9",
    "5/0", "5/1", "5/2", "5/3", "5/4", "5/6", "5/7", "5/8", "5/9",
    "6/0", "6/1", "6/2", "6/3", "6/4", "6/5", "6/7", "6/8", "6/9",
    "7/0", "7/1", "7/2", "7/3", "7/4", "7/5", "7/6", "7/8", "7/9",
    "8/0", "8/1", "8/2", "8/3", "8/4", "8/5", "8/6", "8/7", "8/9",
    "9/0", "9/1", "9/2", "9/3", "9/4", "9/5", "9/6", "9/7", "9/8")


##########################################################################################################################
##########################################################################################################################
#global functions
pub_check_chr_value <- function(ro) {
    inputs <- strsplit(ro, split = ",")[[1]]
    res <- c()
    for(i in inputs){
        if (grepl(":", i)){
            tmp <- strsplit(i, split = ":")[[1]][1]
            if(tmp %in% lis_glo_avalichrnames) {
                res <- c(res, tmp)
            }else {
                res <- c(res, NA)
            }
        }else{
            tmp = as.character(fra_glo_geneindex[which(fra_glo_geneindex$gene_name == i),]$CHROM)
            if(length(tmp) != 1){
                res <- c(res, NA)
            }else{
                res <- c(res, tmp)
            }
        }
    }
    #return
    res
}

pub_check_pos_begin <- function(ro, ro_ext=0){
    inputs <- strsplit(ro, split = ",")[[1]]
    res <- c()
    for(i in inputs){
        if(grepl(":", i)){
            tmp <- strsplit(i, split = ":")[[1]][2]
            if(! grepl("-", tmp)) {
                res <- c(res, NA)
                next
            }
            tmp <- strsplit(tmp, split = "-")[[1]][1]
            tmp <- as.numeric(tmp) - as.numeric(ro_ext)
            res <- c(res, as.character(ifelse(tmp<0, 0, tmp)))
        }else{
            tmp = as.character(fra_glo_geneindex[which(fra_glo_geneindex$gene_name == i),]$begin)
            if(length(tmp) != 1){
                res <- c(res, NA)
            }else{
                tmp <- as.numeric(tmp) - as.numeric(ro_ext)
                res <- c(res, as.character(ifelse(tmp<0, 0, tmp)))
            }
        }
    }
    #return
    res
}

pub_check_pos_end <- function(ro, ro_ext=0) {
    inputs <- strsplit(ro, split = ",")[[1]]
    res <- c()
    for(i in inputs){
        if (grepl(":", i)){
            tmp <- strsplit(i, split = ":")[[1]][2]
            if(! grepl("-", tmp)) {
                res <- c(res, NA)
                next
            }
            tmp <- strsplit(tmp, split = "-")[[1]][2]
            tmp <- as.numeric(tmp) + as.numeric(ro_ext)
            res <- c(res, as.character(tmp))
        }else{
            tmp = as.character(fra_glo_geneindex[which(fra_glo_geneindex$gene_name == i),]$end)
            if(length(tmp) != 1){
                res <- c(res, NA)
            }else{
                tmp <- as.numeric(tmp) + as.numeric(ro_ext)
                res <- c(res, as.character(tmp))
            }
        }
    }
    #return
    res
}

pub_check_sample_name <- function(user_sam) {
    tmp_res <- c()
    res <- c()
    for(i in user_sam){
        tmp_res <- c(tmp_res, strsplit(i, split = ",")[[1]])
    }
    for(i in 1:length(tmp_res)){
        if(grepl("#", tmp_res[i])){
            tmp <- as.character(fra_glo_groupdata[which(fra_glo_groupdata$Group == substr(tmp_res[i], 2, nchar(tmp_res[i]))),]$Name)
            if(length(tmp) != 1) {
                tmp <- NA
            }else{
                tmp <- strsplit(tmp, split = ",")[[1]]
            }
            res <- c(res, tmp)
        }else{
            res <- c(res, tmp_res[i])
        }
    }

    for(i in 1:length(res)){
        tmp <- as.character(fra_glo_metadata[which(fra_glo_metadata$Label == res[i]),]$Accession)
        if(length(tmp) != 1){
            res[i] <- NA
        }else{
            res[i] <- tmp
        }
    }
    #res[! res %in% lis_glo_vcfsamples] <- NA
    #return
    res
}

pub_sub_group <- function(co){
    res <- data.frame()
    if(grepl("\\}", co)){
        tmp1 <- strsplit(co, split = "\\}")[[1]]
        for(i in tmp1){
            group_name <- strsplit(strsplit(i, split = "\\{")[[1]][1], split = ",")[[1]]
            group_name <- group_name[length(group_name)]
            samps <- strsplit(i, split = "\\{")[[1]][2]
            samps <- strsplit(samps, split = ",")[[1]]
            for(j in samps){
                tmp2 <- as.character(fra_glo_metadata[which(fra_glo_metadata$Label == j),]$Accession)
                if(length(tmp2) != 1) { tmp2 <- NA }
                res <- rbind(res, data.frame(c(tmp2), c(group_name)))
            }
        }
        names(res) <- c("Sample", "Group")
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
                res <- rbind(res, data.frame(c(tmp3), c(i)))
            }
        }
        names(res) <- c("Sample", "Group")
    }

    #return
    res
}