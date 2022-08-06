configfile: "config/config.yaml"

include: "AlignAndCall.snakefile"
include: "AlignAndMarkDuplicates.snakefile"


rule all:
	input: 
  
