configfile: "config/config.yaml"
configfile: "config/samples.yaml"

#include: "AlignAndCall.snakefile"
#include: "AlignAndMarkDuplicates.snakefile"

rule all:
	input:
		expand("results/SubsetBamtoChrM/{samples}_chrM.bam", tumor=config["pairings"]),
		expand("results/SubsetBamtoChrM/{samples}_chrM.bai", tumor=config["pairings"])

rule SubsetBamtoChrM:
    input:
    	tumor_filepath = lambda wildcards: config["samples"][wildcards.tumor],
	normal_filepath = lambda wildcards: config["samples"][config["pairings"][wildcards.tumor]]
    output:
    	bam = protected("results/SubsetBamtoChrM/{samples}_chrM.bam"),
	bai = protected("results/SubsetBamtoChrM/{samples}_chrM.bai")
    params:
        gatk = config["gatk_path"],
	contig_name = config["contig_name"]
    log:
    	"logs/SubsetBamtoChrM/{tumors}.txt"
    shell:
	"
	({params.gatk} PrintReads \
      	-L {params.contig_name} \
      	--read-filter MateOnSameContigOrNoMappedMateReadFilter \
      	--read-filter MateUnmappedAndUnmappedReadFilter \
      	-I {input.tumor_filepath} \
      	-output {output.bam}) 2> {log}"
