configfile: "config/config.yaml"
configfile: "config/samples.yaml"

include: "AlignAndCall.snakefile"
include: "AlignAndMarkDuplicates.snakefile"


rule all:
	input: 
  
rule SubsetBamtoChrM:
    input:
    output:
    params:
        gatk = config["gatk_path"],
	contig_name = config["contig_name"]
    logs:
    	"logs/SubsetBamtoChrM/{tumors}.txt"
    shell:
    	"set -e
	({params.gatk} PrintReads \
      	~{"-R " + ref_fasta} \
      	-L {params.contig_name} \
      	--read-filter MateOnSameContigOrNoMappedMateReadFilter \
      	--read-filter MateUnmappedAndUnmappedReadFilter \
        ~{"--gcs-project-for-requester-pays " + requester_pays_project} \
      	-I {input.input_bam} \
      	--read-index ~{input_bai} \
      	-output ~{basename}.bam) 2> {log}"
