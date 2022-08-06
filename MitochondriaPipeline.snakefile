configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
    input:
	expand("results/SubsetBamtoChrM/{sample}.bam", sample=config["samples"]),
	expand("results/SubsetBamtoChrM/{sample}.bai", sample=config["samples"])

rule SubsetBamtoChrM:
    input:
    	tumor_filepath = lambda wildcards: config["samples"][wildcards.tumor],
	normal_filepath = lambda wildcards: config["samples"][config["pairings"][wildcards.tumor]]
    output:
    	bam = protected(expand("results/SubsetBamtoChrM/{sample}.bam", sample=config["samples"])),
	bai = protected(expand("results/SubsetBamtoChrM/{sample}.bai", sample=config["samples"]))
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
	
