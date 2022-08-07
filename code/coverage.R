library(optparse)

option_list <- list(
                make_option(c("--input"), type = "character", help = "Input File"),
                make_option(c("--mean_output"), type = "character", help = "mean_coverage.txt"),
                make_option(c("--median_output"), type = "character", help = "median_coverage.txt"),)

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
args <- commandArgs(TRUE)
options(stringsAsFactors=FALSE, width=160, scipen=999)

input <- opt$input
mean_output <- opt$mean_output
median_output <- opt$median_output


df = read.table(input,skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\t',nrows=1)
write.table(floor(df[,"MEAN_COVERAGE"]), mean_output, quote=F, col.names=F, row.names=F)
write.table(df[,"MEDIAN_COVERAGE"], median_output, quote=F, col.names=F, row.names=F)
