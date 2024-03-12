#packages
library(rjson)
library(shiny)
library(dplyr)
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
fra_glo_metadata_raw <- read.table(path_metadata, header = T, as.is = T, fill=T, sep="\t")
fra_glo_metadata <- fra_glo_metadata_raw[,1:3]
names(fra_glo_metadata) <-  c("Accession", "Label", "Name")
fra_glo_groupdata <- read.table(path_groupdata, col.names = c("Group", "Name"), header = F, as.is = T, sep="\t")
fra_glo_avalichrnames <- read.table(path_fa_index, col.names = c("CHROM", "Length", "V1", "V2", "V3"), as.is = T)
lis_glo_avalichrnames <- fra_glo_avalichrnames$CHROM
fra_glo_samloc <- read.table(path_sam_location, header = F, stringsAsFactors = F, col.names = c("Acession", "location", "lon", "lat"), sep="\t")
if(!is.na(path_sysinfo)){
    fra_glo_sysdata <- pub_read_sysinfo(path_sysinfo)
}else{
    fra_glo_sysdata <- data.frame()
}

if(is.na(path_UIsetting)){
    path_UIsetting <- "./runtime_scripts/empty.json"
    json_glo_UIsetting <- fromJSON(paste(readLines(path_UIsetting), collapse=""))
}else{
    json_glo_UIsetting <- fromJSON(paste(readLines(path_UIsetting), collapse=""))
}

if("ALL" %in% fra_glo_groupdata[,1]){
    fra_glo_groupdata[which(fra_glo_groupdata$Group == "ALL"),2] <- paste(fra_glo_metadata_raw$Accession, collapse=",")
}else{
    fra_glo_groupdata <- rbind(fra_glo_groupdata, c("ALL", paste(fra_glo_metadata_raw$Accession, collapse=",")))
}

text_pub_err_samples <- ""
text_pub_err_groups <- ""
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
fra_hn_detail <- data.frame()
#distance tree
vcfR_dt_orivcf <- "vcfR file"
bool_dt_warning <- F
text_dt_currpara <- ""
text_dt_detail_mess <- ""
text_dt_detail_expl <- ""
fra_dt_detail <- data.frame()
#hapmap
fra_hm_orivcf <- data.frame(c("fra_hm_orivcf"))
fra_hm_drawdata <- data.frame(c("fra_hm_drawdata"))
int_hm_warning <- 0
text_hm_NAwarning <- ""
text_hm_drawsite <- ""
fra_hm_detail <- data.frame()
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

pub_check_sample_name <- function(user_sam, support_raw=F) {
    text_pub_err_samples <<- ""
    text_pub_err_groups <<- ""
    tmp_res <- c()
    res <- c()
    for(i in user_sam){
        tmp_res <- c(tmp_res, strsplit(i, split = ",")[[1]])
    }
    if(length(tmp_res) == 0){
        return(c())
    }
    for(i in 1:length(tmp_res)){
        if(support_raw && "#RAW"==tmp_res[i]){
            tmp <- "#RAW"
            res <- c(res, tmp)
            next
        }
        if(grepl("#", tmp_res[i])){
            tmp <- as.character(fra_glo_groupdata[which(fra_glo_groupdata$Group == substr(tmp_res[i], 2, nchar(tmp_res[i]))),]$Name)
            if(length(tmp) != 1) {
                text_pub_err_groups <<- paste(text_pub_err_groups, tmp_res[i], sep= " ")
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
        if(support_raw && !is.na(res[i]) && res[i]=="#RAW"){
            next
        }
        tmp <- as.character(fra_glo_metadata[which(fra_glo_metadata$Label == res[i]),]$Accession)
        if(length(tmp) != 1){
            text_pub_err_samples <<- paste(text_pub_err_samples, res[i], sep= " ")
            res[i] <- NA
        }else{
            res[i] <- tmp
        }
    }
    #res[! res %in% lis_glo_vcfsamples] <- NA
    #return
    res
}

pub_check_gro_sample_name <- function(co){
    grolis <- pub_sub_rawgroup(co)
    text_pub_err_samples <<- ""
    text_pub_err_groups <<- ""
    if(length(is.na(grolis))==1){
        if(is.na(grolis)){
            return(NA)
        }
    }

    res <- c()
    for(groele in grolis){
        if(grepl("\\}", groele)){
            tmp1 <- strsplit(groele, split = "\\}")[[1]][1]
            samps <- strsplit(tmp1, split = "\\{")[[1]][2]
            samps <- strsplit(samps, split = ",")[[1]]
            for(j in samps){
                tmp2 <- as.character(fra_glo_metadata[which(fra_glo_metadata$Label == j),]$Accession)
                if(length(tmp2) != 1) {
                    tmp2 <- NA
                    text_pub_err_samples <<- paste(text_pub_err_samples, j, sep= " ")
                }
                res <- c(res, tmp2)
            }
        }else{
            if(groele %in% fra_glo_groupdata$Group){
                res <- c(res, groele)
            }else{
                res <- c(res, NA)
                text_pub_err_groups <<- paste(text_pub_err_groups, groele, sep= " ")
            }
        }
    }
    #return
    res
}

pub_sub_group <- function(co){
    grolis <- pub_sub_rawgroup(co)
    if(length(grolis)==1 && is.na(grolis)){
        return(NA)
    }
    res <- data.frame()
    for(groele in grolis){
        if(grepl("\\}", groele)){
            tmp1 <- strsplit(groele, split = "\\}")[[1]][1]
            group_name <- strsplit(strsplit(tmp1, split = "\\{")[[1]][1], split = ",")[[1]][1]
            samps <- strsplit(tmp1, split = "\\{")[[1]][2]
            samps <- strsplit(samps, split = ",")[[1]]
            for(j in samps){
                tmp2 <- as.character(fra_glo_metadata[which(fra_glo_metadata$Label == j),]$Accession)
                if(length(tmp2) != 1) { tmp2 <- NA }
                tmp_fra <- data.frame(c(tmp2), c(group_name))
                names(tmp_fra) <- c("Sample", "Group")
                res <- rbind(res, tmp_fra)
            }
        }else{
            tmp1 <- as.character(fra_glo_groupdata[which(fra_glo_groupdata$Group == groele),]$Name)
            if(length(tmp1) != 1){
                tmp1 <- NA
            }else{
                tmp1 <- strsplit(tmp1, split = ",")[[1]]
            }
            if(length(tmp1)>1 || !is.na(tmp1)){
                for(j in tmp1){
                    tmp2 <- as.character(fra_glo_metadata[which(fra_glo_metadata$Label == as.character(j)),]$Accession)
                    if(length(tmp2) != 1) { tmp2 <- NA }
                    tmp_fra <- data.frame(c(tmp2), c(groele))
                    names(tmp_fra) <- c("Sample", "Group")
                    res <- rbind(res, tmp_fra)
                }
            }else{
                return(NA)
            }
        }
    }

    #return
    res
}

#functions for pub_sub_rawgroup
pub_contain_left_brace <- function(str1){
  length(grep("{",str1,fixed = T)) == 1
}
pub_contain_right_brace <- function(str1){
  length(grep("}",str1,fixed = T)) == 1
}
`%+=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))
`%-=%` = function(e1,e2) eval.parent(substitute(e1 <- e1 - e2))

pub_sub_rawgroup <- function(co){
    outlist <- list()
    j <- 1
    in_middle <- 0
    str_splited <- strsplit(co, ",")[[1]]
    for (i in seq_len(length(str_splited))){
        if (pub_contain_left_brace(str_splited[i])){
            outlist[j] <- str_splited[i]
            if (pub_contain_right_brace(str_splited[i])) j %+=% 1
            else in_middle %+=% 1
        }else if (pub_contain_right_brace(str_splited[i])){
            outlist[j] <- paste0(outlist[j], ",", str_splited[i])
            in_middle %-=% 1
            j %+=% 1
        }else{
            if(in_middle == 1){
                outlist[j] <- paste0(outlist[j], ",", str_splited[i])
            }else{
                outlist[j] <- str_splited[i]
                j %+=% 1
            }
        }
        if (in_middle != 0 & in_middle != 1){
            outlist <- NA
            break
        }
    }
    outlist
}

pub_unique_sample <- function(sam){
    res <- c()
    for(i in 1:length(sam)){
        if(sam[i] %in% res){
            next
        }
        res <- c(res, sam[i])
    }
    res
}
