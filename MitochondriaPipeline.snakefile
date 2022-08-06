configfile: "config/config.yaml"
configfile: "config/samples.yaml"

include: "AlignAndCall.snakefile"
include: "AlignAndMarkDuplicates.snakefile"


rule all:
	input: 
  
rule SubsetBamtoChrM:
    input:
    	tumor_filepath = lambda wildcards: config["samples"][wildcards.tumor],
	normal_filepath = lambda wildcards: config["samples"][config["pairings"][wildcards.tumor]]
    output:
    	
    params:
        gatk = config["gatk_path"],
	contig_name = config["contig_name"]
    logs:
    	"logs/SubsetBamtoChrM/{tumors}.txt"
    shell:
    	"set -e
	
	({params.gatk} PrintReads \
      	-L {params.contig_name} \
      	--read-filter MateOnSameContigOrNoMappedMateReadFilter \
      	--read-filter MateUnmappedAndUnmappedReadFilter \
      	-I {input.tumor_filepath} \
      	--read-index {input.input_bai} \
      	-output {output.}) 2> {log}"
