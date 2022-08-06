import os
import re

configfile: "config/config.yaml"
configfile: "config/samples.yaml"

#include: "AlignAndCall.snakefile"
#include: "AlignAndMarkDuplicates.snakefile"

rule all:
    input:
	expand("results/SubsetBamtoChrM/{sample}.bam", sample=config["samples"]),
	expand("results/SubsetBamtoChrM/{sample}.bai", sample=config["samples"])

rule SubsetBamtoChrM:
    input:
    	tumor_filepath = lambda wildcards: config["samples"][wildcards.tumor],
	normal_filepath = lambda wildcards: config["samples"][config["pairings"][wildcards.tumor]]
    output:
    	bam = protected("results/SubsetBamtoChrM/{tumor}_chrM.bam"),
	bai = protected("results/SubsetBamtoChrM/{tumor}_chrM.bai")
    params:
        gatk = config["gatk_path"],
	contig_name = config["contig_name"]
    log:
    	"logs/SubsetBamtoChrM/{tumor}.txt"
    shell:
	"({params.gatk} PrintReads \
      	-L {params.contig_name} \
      	--read-filter MateOnSameContigOrNoMappedMateReadFilter \
      	--read-filter MateUnmappedAndUnmappedReadFilter \
      	-I {input.tumor_filepath} \
      	-O {output.bam}) 2> {log}"
	
