# filepaths
path_fasta <- "/data2/rawdata2/shinyZavitan/151210_zavitan_v2_pseudomolecules.fa"
path_fa_index <- "/data2/rawdata2/shinyZavitan/151210_zavitan_v2_pseudomolecules.fa.fai"
path_gff3 <- "/data2/rawdata2/shinyZavitan/Zavitan.sort.gff3.gz"
path_metadata <- "/data2/rawdata2/shinyZavitan/sample_name.txt"
path_geneindex <- "/data2/rawdata2/shinyZavitan/geneinfo.txt"
path_vcf <- "/data2/rawdata2/shinyZavitan/Zavitan.ann.bcf.gz"
path_groupdata <- "/data2/rawdata2/shinyZavitan//group_info.txt"
path_sam_location <- "/data2/rawdata2/shinyZavitan/location_info.txt"

path_sysinfo <- "/data2/rawdata2/shinyZavitan/sys_info.txt" #The ONLY table that HAVE header, also the ONLY optional table with default value NA

path_UIsetting <- "./sample/Zavitan.json"
#tool application paths
path_bcftools <- "/data/user/shinyug/bin/bcftools"
path_samtools <- "samtools"
path_seqkit <- "/data/user/shinyug/bin/seqkit"
path_tabix <- "/home/wangzh/bin/htslib/bin/tabix"

#system
.libPaths(c(.libPaths(), "/home/wangzh/R/x86_64-redhat-linux-gnu-library/3.4"))