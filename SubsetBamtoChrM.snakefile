import os

configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
    input:
        expand("results/SubsetBamtoChrM/{tumor}/{tumor}.bam", tumor=config["pairings"]),
        expand("results/RevertSam/{tumor}/{tumor}.bam", tumor=config["pairings"]),
        expand("results/AlignAndMarkDuplicates/{tumor}/{tumor}_mba.bam", tumor=config["pairings"]),
        expand("results/AlignAndMarkDuplicates/{tumor}/{tumor}_md.bam", tumor=config["pairings"]),
        expand("results/AlignAndMarkDuplicates/{tumor}/{tumor}.bam", tumor=config["pairings"]),
        expand("results/AlignAndMarkDuplicates/{tumor}/{tumor}.metrics", tumor=config["pairings"])
        
rule SubsetBamtoChrM:
    input:
        tumor_filepath = lambda wildcards: config["samples"][wildcards.tumor],
        #normal_filepath = lambda wildcards: config["samples"][config["pairings"][wildcards.tumor]]
    output:
        bam = "results/SubsetBamtoChrM/{tumor}/{tumor}.bam",
        bai = "results/SubsetBamtoChrM/{tumor}/{tumor}.bai"
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

rule RevertSam:
    input:
        bam = "results/SubsetBamtoChrM/{tumor}/{tumor}.bam"
    output:
        bam = "results/RevertSam/{tumor}/{tumor}.bam",
    params:
        java = config["java"],
        picard_jar = config["picard_jar"]
    log:
        "logs/RevertSam/{tumor}.txt"
    shell:
        "({params.java} -Xmx1000m -jar {params.picard_jar} \
        RevertSam \
        INPUT={input.bam} \
        OUTPUT_BY_READGROUP=false \
        OUTPUT={output.bam} \
        VALIDATION_STRINGENCY=LENIENT \
        ATTRIBUTE_TO_CLEAR=FT \
        ATTRIBUTE_TO_CLEAR=CO \
        SORT_ORDER=queryname \
        RESTORE_ORIGINAL_QUALITIES=false) 2> {log}"

rule AlignAndMarkDuplicates:
    input:
        bam = "results/RevertSam/{tumor}/{tumor}.bam",
    output:
        mba_bam = "results/AlignAndMarkDuplicates/{tumor}/{tumor}_mba.bam",
        md_bam = "results/AlignAndMarkDuplicates/{tumor}/{tumor}_md.bam",
        bam = "results/AlignAndMarkDuplicates/{tumor}/{tumor}.bam",
        metrics = "results/AlignAndMarkDuplicates/{tumor}/{tumor}.metrics"
    params:
        bwa = config["bwa"],
        java = config["java"],
        picard_jar = config["picard_jar"],
        gatk = config["gatk_path"],
        reference_genome = config["mt_ref"]
    log:
        "logs/AlignAndMarkDuplicates/{tumor}.txt"
    shell:
        "print(os.path.basename({input.bam})"
