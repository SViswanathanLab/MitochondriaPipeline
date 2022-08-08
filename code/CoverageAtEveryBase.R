library(optparse)

option_list <- list(
                make_option(c("--control_region_shifted"), type = "character", help = "Input File"),
                make_option(c("--non_control_region"), type = "character", help = "mean_coverage.txt"),
                make_option(c("--output"), type = "character", help = "per_base_coverage.tsv"))

parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
args <- commandArgs(TRUE)
options(stringsAsFactors=FALSE, width=160, scipen=999)

control_region_shifted <- opt$control_region_shifted
non_control_region <- opt$non_control_region
per_base_coverage <- opt$output


shift_back = function(x) {
   if (x < 8570) {
  return(x + 8000)
} else {
  return (x - 8569)}
}

control_region_shifted = read.table(control_region_shifted, header=T,fill=T)
shifted_back = sapply(control_region_shifted[,"pos"], shift_back)
control_region_shifted[,"pos"] = shifted_back

beginning = subset(control_region_shifted, control_region_shifted[,'pos']<8000)
end = subset(control_region_shifted, control_region_shifted[,'pos']>8000)

non_control_region = read.table(non_control_region, header=T)
combined_table = rbind(beginning, non_control_region, end)
write.table(combined_table, per_base_coverage, row.names=F, col.names=T, quote=F, sep="\t")
